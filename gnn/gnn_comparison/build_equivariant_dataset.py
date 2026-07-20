"""Build the shared 3D PyG dataset used by all four equivariant models.

The repository data.csv is deterministically shuffled with seed 42 and split
into the same 4,392/1,098 train/test sets used in the reported experiments.
Molecular descriptors are fitted on the training split only.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import pickle
import random
import sys
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import Data


SCRIPT_DIR = Path(__file__).resolve().parent
GRAPH_TRANSFORMER_DIR = SCRIPT_DIR.parent / "graph_transformer"
sys.path.insert(0, str(GRAPH_TRANSFORMER_DIR))

from features import get_atom_features, get_bond_features  # noqa: E402
from ml_features import N_ML_FEATURES, compute_single, get_feature_names  # noqa: E402


RDLogger.DisableLog("rdApp.*")

SOLVENT_TO_ID = {
    "[2H]C(Cl)(Cl)Cl": 0,
    "[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]": 1,
    "[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]": 2,
    "[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]": 3,
    "[2H]C([2H])([2H])C#N": 4,
    "[2H]OC([2H])([2H])[2H]": 5,
    "[2H]C([2H])(Cl)Cl": 6,
    "[2H]C1([2H])OC([2H])([2H])C([2H])([2H])C1([2H])[2H]": 7,
    "[2H]c1c([2H])c([2H])c(C([2H])([2H])[2H])c([2H])c1[2H]": 8,
    "[2H]O[2H]": 9,
}


def read_rows(path: Path, limit: int | None = None) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {"Smiles", "solvent", "B_count", "ppm_values"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path} must contain columns: {sorted(required)}")
    return rows[:limit] if limit else rows


def parse_targets(row: dict[str, str], num_boron: int) -> torch.Tensor:
    values = [
        float(item.strip())
        for item in row["ppm_values"].replace(";", ",").split(",")
        if item.strip()
    ]
    if len(values) == 1 and num_boron > 1:
        values *= num_boron
    if len(values) != num_boron:
        raise ValueError(
            f"target count mismatch for {row['Smiles']}: "
            f"expected {num_boron}, got {len(values)}"
        )
    return torch.tensor(values, dtype=torch.float32)


def _embed(mol: Chem.Mol, seed: int, timeout: int, random_coords: bool):
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed % 2_147_483_647)
    params.maxIterations = 200
    params.timeout = int(timeout)
    params.useRandomCoords = random_coords
    conformer_id = AllChem.EmbedMolecule(mol_h, params)
    return mol_h, conformer_id


def generate_coordinates(
    mol: Chem.Mol, seed: int, timeout: int
) -> tuple[torch.Tensor, str]:
    """Generate one deterministic conformer without ever dropping a sample."""
    mol_h = None
    conformer_id = -1
    mode = "etkdg"
    try:
        mol_h, conformer_id = _embed(mol, seed, timeout, random_coords=False)
        if conformer_id < 0:
            mode = "random_fallback"
            mol_h, conformer_id = _embed(mol, seed, timeout, random_coords=True)
    except Exception:
        conformer_id = -1

    if conformer_id >= 0 and mol_h is not None:
        quality = mode + "_unoptimized"
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol_h):
                AllChem.MMFFOptimizeMolecule(mol_h, confId=conformer_id, maxIters=200)
                quality = mode + "_mmff"
            elif AllChem.UFFHasAllMoleculeParams(mol_h):
                AllChem.UFFOptimizeMolecule(mol_h, confId=conformer_id, maxIters=200)
                quality = mode + "_uff"
        except Exception:
            pass

        conformer = mol_h.GetConformer(conformer_id)
        # AddHs appends atoms, so the first mol.GetNumAtoms() coordinates retain
        # the exact atom order used by the graph and mask_b.
        positions = [
            conformer.GetAtomPosition(index) for index in range(mol.GetNumAtoms())
        ]
        pos = torch.tensor(
            [[point.x, point.y, point.z] for point in positions],
            dtype=torch.float32,
        )
        return pos, quality

    # Last-resort deterministic fallback: retain the sample and mark it clearly.
    mol_2d = Chem.Mol(mol)
    AllChem.Compute2DCoords(mol_2d)
    conformer = mol_2d.GetConformer()
    positions = [conformer.GetAtomPosition(i) for i in range(mol_2d.GetNumAtoms())]
    pos = torch.tensor(
        [[point.x, point.y, 0.0] for point in positions], dtype=torch.float32
    )
    return pos, "2d_fallback"


def build_spatial_edges(
    mol: Chem.Mol, pos: torch.Tensor, cutoff: float
) -> tuple[torch.Tensor, torch.Tensor]:
    """Build directed radius edges within each covalent component.

    Covalent edges are always retained. The final 11-dimensional edge feature
    contains the original 10 bond features plus an explicit is_bond flag.
    """
    bond_features: dict[tuple[int, int], torch.Tensor] = {}
    edge_pairs: set[tuple[int, int]] = set()
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        feature = get_bond_features(bond)
        bond_features[(i, j)] = feature
        bond_features[(j, i)] = feature
        edge_pairs.update({(i, j), (j, i)})

    for component in Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False):
        indices = torch.tensor(component, dtype=torch.long)
        component_pos = pos[indices]
        distances = torch.cdist(component_pos, component_pos)
        src_local, dst_local = torch.where(
            (distances <= cutoff) & (distances > 1.0e-8)
        )
        for src, dst in zip(src_local.tolist(), dst_local.tolist()):
            edge_pairs.add((int(indices[src]), int(indices[dst])))

    ordered_pairs = sorted(edge_pairs)
    if not ordered_pairs:
        return (
            torch.empty((2, 0), dtype=torch.long),
            torch.empty((0, 11), dtype=torch.float32),
        )

    edge_index = torch.tensor(ordered_pairs, dtype=torch.long).t().contiguous()
    edge_attr = []
    zero_bond = torch.zeros(10, dtype=torch.float32)
    for pair in ordered_pairs:
        is_bond = pair in bond_features
        bond_attr = bond_features[pair] if is_bond else zero_bond
        edge_attr.append(
            torch.cat([bond_attr, torch.tensor([float(is_bond)])], dim=0)
        )
    return edge_index, torch.stack(edge_attr)


def fit_ml_features(
    train_rows: list[dict[str, str]], test_rows: list[dict[str, str]]
) -> tuple[np.ndarray, np.ndarray, StandardScaler, np.ndarray]:
    train_raw = np.stack([compute_single(row["Smiles"]) for row in train_rows])
    test_raw = np.stack([compute_single(row["Smiles"]) for row in test_rows])
    medians = np.nanmedian(train_raw, axis=0)
    medians = np.where(np.isnan(medians), 0.0, medians)
    train_raw = np.where(np.isnan(train_raw), medians, train_raw)
    test_raw = np.where(np.isnan(test_raw), medians, test_raw)
    scaler = StandardScaler().fit(train_raw)
    return scaler.transform(train_raw), scaler.transform(test_raw), scaler, medians


def row_to_data(
    row: dict[str, str],
    ml_features: np.ndarray,
    row_index: int,
    split_name: str,
    cutoff: float,
    seed: int,
    embed_timeout: int,
) -> Data:
    mol = Chem.MolFromSmiles(row["Smiles"])
    if mol is None:
        raise ValueError(f"invalid SMILES at {split_name} row {row_index}")
    try:
        AllChem.ComputeGasteigerCharges(mol)
    except Exception:
        pass

    x = torch.stack([get_atom_features(atom) for atom in mol.GetAtoms()])
    mask_b = torch.tensor(
        [atom.GetSymbol() == "B" for atom in mol.GetAtoms()], dtype=torch.bool
    )
    num_boron = int(mask_b.sum())
    expected_boron = int(row["B_count"])
    if num_boron != expected_boron:
        raise ValueError(
            f"boron count mismatch at {split_name} row {row_index}: "
            f"CSV={expected_boron}, RDKit={num_boron}"
        )

    pos, geometry_quality = generate_coordinates(
        mol, seed=seed + row_index, timeout=embed_timeout
    )
    # Use the arithmetic mean of all modeled atom coordinates as the molecular
    # center. This changes only the global translation, never distances or edges.
    pos = pos - pos.mean(dim=0, keepdim=True)
    edge_index, edge_attr = build_spatial_edges(mol, pos, cutoff=cutoff)
    components = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
    component_id = torch.empty(mol.GetNumAtoms(), dtype=torch.long)
    for component_index, component in enumerate(components):
        component_id[list(component)] = component_index

    return Data(
        x=x,
        pos=pos,
        edge_index=edge_index,
        edge_attr=edge_attr,
        mask_b=mask_b,
        y_b=parse_targets(row, num_boron),
        solvent_id=torch.tensor(
            [SOLVENT_TO_ID.get(row["solvent"].strip(), 10)], dtype=torch.long
        ),
        ml_global_features=torch.tensor(ml_features, dtype=torch.float32).view(1, -1),
        component_id=component_id,
        smiles=row["Smiles"],
        row_index=torch.tensor([row_index], dtype=torch.long),
        geometry_quality=geometry_quality,
    )


def build_split(
    rows: list[dict[str, str]],
    ml_features: np.ndarray,
    split_name: str,
    cutoff: float,
    seed: int,
    embed_timeout: int,
    workers: int,
) -> tuple[list[Data], Counter]:
    dataset = []
    quality_counts: Counter = Counter()
    tasks = [
        (
            row,
            feature_row,
            index,
            split_name,
            cutoff,
            seed,
            embed_timeout,
        )
        for index, (row, feature_row) in enumerate(zip(rows, ml_features))
    ]
    if workers <= 1:
        for index, data in enumerate(map(_build_row_task, tasks)):
            dataset.append(data)
            quality_counts[data.geometry_quality] += 1
            if (index + 1) % 100 == 0 or index + 1 == len(rows):
                print(f"  {split_name}: {index + 1}/{len(rows)}")
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        max_buffer = max(workers * 8, workers)
        pending = {}
        ready = {}
        next_submit = 0
        next_emit = 0

        def submit_until_full():
            nonlocal next_submit
            while (
                next_submit < len(tasks)
                and len(pending) + len(ready) < max_buffer
            ):
                # Returning raw bytes avoids torch.multiprocessing creating one
                # shared-memory file descriptor per tensor for large datasets.
                future = executor.submit(
                    _build_row_task_serialized, tasks[next_submit]
                )
                pending[future] = next_submit
                next_submit += 1

        submit_until_full()
        try:
            while next_emit < len(tasks):
                if next_emit not in ready:
                    completed, _ = wait(pending, return_when=FIRST_COMPLETED)
                    for future in completed:
                        index = pending.pop(future)
                        ready[index] = pickle.loads(future.result())
                    submit_until_full()
                    continue

                data = ready.pop(next_emit)
                dataset.append(data)
                quality_counts[data.geometry_quality] += 1
                next_emit += 1
                if next_emit % 100 == 0 or next_emit == len(rows):
                    print(f"  {split_name}: {next_emit}/{len(rows)}")
                submit_until_full()
        finally:
            executor.shutdown(wait=True, cancel_futures=True)
    return dataset, quality_counts


def _build_row_task(task) -> Data:
    (
        row,
        feature_row,
        row_index,
        split_name,
        cutoff,
        seed,
        embed_timeout,
    ) = task
    return row_to_data(
        row=row,
        ml_features=feature_row,
        row_index=row_index,
        split_name=split_name,
        cutoff=cutoff,
        seed=seed,
        embed_timeout=embed_timeout,
    )


def _build_row_task_serialized(task) -> bytes:
    return pickle.dumps(_build_row_task(task), protocol=pickle.HIGHEST_PROTOCOL)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-csv", type=Path, default=GRAPH_TRANSFORMER_DIR / "data.csv"
    )
    parser.add_argument(
        "--train-size",
        type=int,
        default=4392,
        help="Number of shuffled rows assigned to the fixed training split.",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=SCRIPT_DIR / "processed_equivariant"
    )
    parser.add_argument("--cutoff", type=float, default=5.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--embed-timeout", type=int, default=5)
    parser.add_argument(
        "--workers",
        type=int,
        default=min(16, os.cpu_count() or 1),
        help="Number of parallel CPU processes used for conformer generation.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Use only the first N source rows, then make an 80/20 smoke split.",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    train_output = args.output_dir / "train_3d.pt"
    test_output = args.output_dir / "test_3d.pt"
    if not args.force and (train_output.exists() or test_output.exists()):
        raise FileExistsError(
            f"output already exists under {args.output_dir}; use --force to overwrite"
        )

    rows = read_rows(args.data_csv, args.limit)
    random.Random(args.seed).shuffle(rows)
    train_size = int(len(rows) * 0.8) if args.limit else args.train_size
    if not 1 <= train_size < len(rows):
        raise ValueError(
            f"train size must be between 1 and {len(rows) - 1}; got {train_size}"
        )
    train_rows = rows[:train_size]
    test_rows = rows[train_size:]
    print(f"Fixed split: train={len(train_rows)}, test={len(test_rows)}")
    train_ml, test_ml, scaler, medians = fit_ml_features(train_rows, test_rows)

    print("Building training conformers...")
    train_dataset, train_quality = build_split(
        train_rows,
        train_ml,
        "train",
        args.cutoff,
        args.seed,
        args.embed_timeout,
        args.workers,
    )
    print("Building test conformers...")
    test_dataset, test_quality = build_split(
        test_rows,
        test_ml,
        "test",
        args.cutoff,
        args.seed + 1_000_000,
        args.embed_timeout,
        args.workers,
    )

    torch.save(train_dataset, train_output)
    torch.save(test_dataset, test_output)
    with (args.output_dir / "ml_feature_scaler_3d.pkl").open("wb") as handle:
        pickle.dump({"scaler": scaler, "medians": medians}, handle)

    metadata = {
        "source_csv": str(args.data_csv.resolve()),
        "split_strategy": (
            f"random.Random({args.seed}).shuffle; first {train_size} rows for train"
        ),
        "train_samples": len(train_dataset),
        "test_samples": len(test_dataset),
        "cutoff_angstrom": args.cutoff,
        "seed": args.seed,
        "embed_timeout_seconds": args.embed_timeout,
        "workers": args.workers,
        "node_feature_dim": int(train_dataset[0].x.shape[1]),
        "edge_feature_dim": int(train_dataset[0].edge_attr.shape[1]),
        "ml_feature_dim": N_ML_FEATURES,
        "ml_feature_names": get_feature_names(),
        "train_geometry_quality": dict(train_quality),
        "test_geometry_quality": dict(test_quality),
        "coordinate_centering": "per-molecule arithmetic mean shifted to origin",
        "notes": "ML scaler fitted on training split only; no samples dropped.",
    }
    with (args.output_dir / "dataset_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)
    print(json.dumps(metadata, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
