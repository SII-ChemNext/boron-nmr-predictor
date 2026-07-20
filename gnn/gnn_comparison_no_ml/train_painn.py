"""Train the PaiNN ablation without prior ML features."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
WITH_ML_DIR = SCRIPT_DIR.parent / "gnn_comparison"
sys.path.insert(0, str(WITH_ML_DIR))
sys.path.insert(0, str(SCRIPT_DIR / "models"))

from equivariant_training import TrainingConfig, run_experiment  # noqa: E402
from model_painn import BoronPaiNNWithoutPriorML  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--epochs", type=int, default=150)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--learning-rate", type=float, default=3.0e-4)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--data-dir", type=Path, default=WITH_ML_DIR / "processed_equivariant")
    parser.add_argument("--run-name", default="split_current")
    return parser.parse_args()


def main():
    args = parse_args()
    if Path(args.run_name).name != args.run_name or args.run_name in {"", ".", ".."}:
        raise ValueError("--run-name must be one directory name without slashes")
    config = TrainingConfig(
        seed=args.seed, folds=args.folds, epochs=args.epochs,
        batch_size=args.batch_size, learning_rate=args.learning_rate,
        num_workers=args.num_workers, device=args.device,
    )
    model_config = {
        "node_in_dim": 58, "edge_in_dim": 11, "hidden_dim": 128,
        "num_layers": 4, "num_rbf": 32, "cutoff": 5.0,
        "dropout": 0.05, "num_solvents": 11, "solvent_dim": 32,
    }
    run_experiment(
        model_name="BoronPaiNN_WithoutPriorML",
        model_class=BoronPaiNNWithoutPriorML,
        model_kwargs=model_config,
        experiment_dir=SCRIPT_DIR,
        data_dir=args.data_dir.expanduser().resolve(),
        config=config,
        output_dir=SCRIPT_DIR / "painn_results" / args.run_name,
    )


if __name__ == "__main__":
    main()
