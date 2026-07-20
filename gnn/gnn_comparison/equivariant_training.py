"""Shared fixed-split training/evaluation code for equivariant experiments."""

from __future__ import annotations

import csv
import json
import math
import random
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import torch
import torch.nn.functional as F
from sklearn.model_selection import KFold
from torch_geometric.loader import DataLoader


@dataclass
class TrainingConfig:
    seed: int = 42
    folds: int = 5
    epochs: int = 150
    batch_size: int = 16
    learning_rate: float = 3.0e-4
    weight_decay: float = 1.0e-5
    scheduler_patience: int = 10
    early_stopping_patience: int = 30
    num_workers: int = 0
    device: str = "cuda:0"


class Tee:
    def __init__(self, path: Path):
        self.terminal = sys.stdout
        self.file = path.open("w", encoding="utf-8")

    def write(self, message: str):
        self.terminal.write(message)
        self.file.write(message)
        self.file.flush()

    def flush(self):
        self.terminal.flush()
        self.file.flush()

    def close(self):
        self.file.close()


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def resolve_device(requested: str) -> torch.device:
    if requested.startswith("cuda") and not torch.cuda.is_available():
        print("CUDA is unavailable; falling back to CPU.")
        return torch.device("cpu")
    return torch.device(requested)


def load_fixed_datasets(data_dir: Path):
    train_path = data_dir / "train_3d.pt"
    test_path = data_dir / "test_3d.pt"
    if not train_path.exists() or not test_path.exists():
        raise FileNotFoundError(
            "Shared 3D datasets are missing. Run: "
            "python data/build_equivariant_dataset.py"
        )
    train_dataset = torch.load(train_path, map_location="cpu", weights_only=False)
    test_dataset = torch.load(test_path, map_location="cpu", weights_only=False)
    return train_dataset, test_dataset


def sorted_node_loss(
    predictions: torch.Tensor,
    targets: torch.Tensor,
    boron_batch: torch.Tensor,
) -> torch.Tensor:
    total = predictions.new_zeros(())
    count = 0
    for graph_index in torch.unique(boron_batch):
        mask = boron_batch == graph_index
        pred_sorted = torch.sort(predictions[mask]).values
        target_sorted = torch.sort(targets[mask]).values
        total = total + F.huber_loss(
            pred_sorted, target_sorted, delta=1.0, reduction="sum"
        )
        count += int(mask.sum())
    return total / max(count, 1)


def train_epoch(model, loader, optimizer, device: torch.device) -> float:
    model.train()
    total_loss = 0.0
    total_boron = 0
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad(set_to_none=True)
        predictions = model(batch)
        boron_batch = batch.batch[batch.mask_b]
        loss = sorted_node_loss(predictions, batch.y_b, boron_batch)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()
        count = int(batch.mask_b.sum())
        total_loss += float(loss) * count
        total_boron += count
    return total_loss / max(total_boron, 1)


def metrics(predictions: torch.Tensor, targets: torch.Tensor) -> dict[str, float]:
    errors = predictions - targets
    mae = torch.mean(torch.abs(errors)).item()
    rmse = torch.sqrt(torch.mean(errors * errors)).item()
    ss_res = torch.sum(errors * errors)
    centered = targets - torch.mean(targets)
    ss_tot = torch.sum(centered * centered)
    r2 = (1.0 - ss_res / ss_tot).item() if float(ss_tot) > 0 else 0.0
    return {"mae": mae, "rmse": rmse, "r2": r2}


@torch.no_grad()
def evaluate(model, loader, device: torch.device) -> dict[str, float]:
    model.eval()
    all_predictions = []
    all_targets = []
    for batch in loader:
        batch = batch.to(device)
        predictions = model(batch)
        boron_batch = batch.batch[batch.mask_b]
        for graph_index in torch.unique(boron_batch):
            mask = boron_batch == graph_index
            all_predictions.append(torch.sort(predictions[mask]).values.cpu())
            all_targets.append(torch.sort(batch.y_b[mask]).values.cpu())
    return metrics(torch.cat(all_predictions), torch.cat(all_targets))


def save_checkpoint(path: Path, model, model_kwargs: dict, fold: int) -> None:
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "model_kwargs": model_kwargs,
            "fold": fold,
        },
        path,
    )


def load_checkpoint(path: Path, model_class, device: torch.device):
    checkpoint = torch.load(path, map_location=device, weights_only=False)
    model = model_class(**checkpoint["model_kwargs"]).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    return model


@torch.no_grad()
def evaluate_ensemble(
    models,
    test_dataset,
    device: torch.device,
    prediction_path: Path,
) -> dict[str, float]:
    loader = DataLoader(test_dataset, batch_size=1, shuffle=False)
    flat_predictions = []
    flat_targets = []
    rows = []
    for batch in loader:
        batch = batch.to(device)
        fold_predictions = [torch.sort(model(batch)).values for model in models]
        prediction = torch.stack(fold_predictions).mean(dim=0).cpu()
        target = torch.sort(batch.y_b).values.cpu()
        flat_predictions.append(prediction)
        flat_targets.append(target)
        errors = torch.abs(prediction - target)
        rows.append(
            {
                "row_index": int(batch.row_index.view(-1)[0]),
                "smiles": batch.smiles[0],
                "targets_ppm": ";".join(f"{value:.8g}" for value in target.tolist()),
                "predictions_ppm": ";".join(
                    f"{value:.8g}" for value in prediction.tolist()
                ),
                "absolute_errors_ppm": ";".join(
                    f"{value:.8g}" for value in errors.tolist()
                ),
                "geometry_quality": batch.geometry_quality[0],
            }
        )

    with prediction_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return metrics(torch.cat(flat_predictions), torch.cat(flat_targets))


def run_experiment(
    model_name: str,
    model_class,
    model_kwargs: dict,
    experiment_dir: Path,
    data_dir: Path,
    config: TrainingConfig,
    output_dir: Path | None = None,
) -> dict:
    output_dir = output_dir or (experiment_dir / "outputs")
    missing_data_files = [
        path for path in (data_dir / "train_3d.pt", data_dir / "test_3d.pt")
        if not path.exists()
    ]
    if missing_data_files:
        raise FileNotFoundError(
            "missing fixed-split data files: "
            + ", ".join(str(path) for path in missing_data_files)
        )
    checkpoint_dir = output_dir / "trained_models"
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    try:
        # Atomic creation prevents two simultaneous jobs from sharing and
        # overwriting the same run directory.
        output_dir.mkdir(exist_ok=False)
    except FileExistsError as error:
        raise FileExistsError(
            f"run output already exists: {output_dir}. "
            "Choose a new --run-name to keep experiments isolated."
        ) from error
    checkpoint_dir.mkdir()
    tee = Tee(output_dir / "training.log")
    original_stdout = sys.stdout
    sys.stdout = tee
    try:
        seed_everything(config.seed)
        device = resolve_device(config.device)
        train_dataset, test_dataset = load_fixed_datasets(data_dir)
        print(f"Model: {model_name}")
        print(f"Device: {device}")
        print(f"Data directory: {data_dir.resolve()}")
        print(f"Fixed split: train={len(train_dataset)}, test={len(test_dataset)}")
        print(f"Configuration: {json.dumps(asdict(config), ensure_ascii=False)}")

        splitter = KFold(
            n_splits=config.folds, shuffle=True, random_state=config.seed
        )
        indices = np.arange(len(train_dataset))
        fold_metrics = []
        checkpoint_paths = []

        for fold, (train_indices, val_indices) in enumerate(
            splitter.split(indices), start=1
        ):
            print(f"\nFold {fold}/{config.folds}")
            train_subset = [train_dataset[index] for index in train_indices]
            val_subset = [train_dataset[index] for index in val_indices]
            train_loader = DataLoader(
                train_subset,
                batch_size=config.batch_size,
                shuffle=True,
                drop_last=False,
                num_workers=config.num_workers,
            )
            val_loader = DataLoader(
                val_subset,
                batch_size=config.batch_size,
                shuffle=False,
                drop_last=False,
                num_workers=config.num_workers,
            )
            model = model_class(**model_kwargs).to(device)
            if fold == 1:
                parameters = sum(p.numel() for p in model.parameters())
                print(f"Trainable parameters: {parameters:,}")
            optimizer = torch.optim.AdamW(
                model.parameters(),
                lr=config.learning_rate,
                weight_decay=config.weight_decay,
            )
            scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                optimizer,
                mode="min",
                factor=0.5,
                patience=config.scheduler_patience,
            )
            best = {"mae": math.inf, "rmse": math.inf, "r2": -math.inf}
            epochs_without_improvement = 0
            checkpoint_path = checkpoint_dir / f"model_fold_{fold}.pth"

            for epoch in range(1, config.epochs + 1):
                train_loss = train_epoch(model, train_loader, optimizer, device)
                val_result = evaluate(model, val_loader, device)
                scheduler.step(val_result["mae"])
                if val_result["mae"] < best["mae"]:
                    best = val_result
                    epochs_without_improvement = 0
                    save_checkpoint(
                        checkpoint_path, model, model_kwargs=model_kwargs, fold=fold
                    )
                else:
                    epochs_without_improvement += 1
                if epoch == 1 or epoch % 10 == 0:
                    print(
                        f"Epoch {epoch:03d} | loss={train_loss:.4f} | "
                        f"val_mae={val_result['mae']:.4f} | "
                        f"val_rmse={val_result['rmse']:.4f} | "
                        f"val_r2={val_result['r2']:.4f}"
                    )
                if epochs_without_improvement >= config.early_stopping_patience:
                    print(f"Early stopping at epoch {epoch}")
                    break

            print(f"Fold {fold} best: {best}")
            fold_metrics.append(best)
            checkpoint_paths.append(checkpoint_path)

        cv_summary = {
            key: {
                "mean": float(np.mean([result[key] for result in fold_metrics])),
                "std": float(np.std([result[key] for result in fold_metrics])),
            }
            for key in ("mae", "rmse", "r2")
        }
        models = [
            load_checkpoint(path, model_class=model_class, device=device)
            for path in checkpoint_paths
        ]
        test_metrics = evaluate_ensemble(
            models,
            test_dataset,
            device=device,
            prediction_path=output_dir / "test_predictions.csv",
        )
        result = {
            "model": model_name,
            "data_dir": str(data_dir.resolve()),
            "fixed_split": {
                "train_samples": len(train_dataset),
                "test_samples": len(test_dataset),
            },
            "training_config": asdict(config),
            "model_config": model_kwargs,
            "fold_metrics": fold_metrics,
            "cv_summary": cv_summary,
            "test_ensemble": test_metrics,
        }
        with (output_dir / "results.json").open("w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2, ensure_ascii=False)
        print(f"\nCV summary: {cv_summary}")
        print(f"Final fixed-test ensemble: {test_metrics}")
        return result
    finally:
        sys.stdout = original_stdout
        tee.close()
