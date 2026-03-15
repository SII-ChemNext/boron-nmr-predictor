import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from torch_geometric.data import Batch
import numpy as np
import random
import optuna
from optuna.pruners import MedianPruner
from sklearn.model_selection import KFold
import sys
import os
from rdkit import Chem

# Import modules from the current directory
from models.model_v3 import BoronNMRNet_V3
from features import get_atom_features, get_bond_features

# =============================================================================
# Configuration
# =============================================================================
DATA_PATH = 'processed_boron_dataset_v3.pt'
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
SEED = 42

# Optuna configuration
N_TRIALS = 50           # number of trials
N_FOLDS = 5             # 5-fold cross validation
EPOCHS_PER_TRIAL = 30   # maximum epochs per trial
TIMEOUT = 28800         # 8-hour timeout (for all 50 trials)

# Fix random seeds
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)

# =============================================================================
# Utility functions
# =============================================================================
def get_dims():
    try:
        m = Chem.MolFromSmiles('CB')
        return get_atom_features(m.GetAtoms()[0]).shape[0], get_bond_features(m.GetBonds()[0]).shape[0]
    except:
        return 0, 0

def custom_collate_fn(batch):
    return Batch.from_data_list(batch)

def compute_sorted_loss(all_preds, all_targets, batch_indices):
    """Node-level loss aggregation"""
    unique_mols = torch.unique(batch_indices)
    total_loss = 0
    total_atoms = 0

    for mol_idx in unique_mols:
        mask = (batch_indices == mol_idx)
        p, t = all_preds[mask], all_targets[mask]
        p_sorted, _ = torch.sort(p)
        t_sorted, _ = torch.sort(t)
        mol_loss = F.huber_loss(p_sorted, t_sorted, delta=1.0, reduction='sum')
        total_loss += mol_loss
        total_atoms += len(p)

    return total_loss / total_atoms if total_atoms > 0 else total_loss

def train_one_epoch(model, loader, optimizer):
    model.train()
    total_loss = 0
    total_atoms = 0

    for batch in loader:
        batch = batch.to(DEVICE)
        optimizer.zero_grad()
        ml_feats = getattr(batch, 'ml_global_features', None)

        preds = model(
            batch.x,
            batch.edge_index,
            batch.edge_attr,
            batch.solvent_id.squeeze(),
            batch.mask_b,
            batch.batch,
            ml_global_features=ml_feats
        )

        b_batch_indices = batch.batch[batch.mask_b]
        loss = compute_sorted_loss(preds, batch.y_b, b_batch_indices)

        loss.backward()
        optimizer.step()

        num_b_atoms = batch.mask_b.sum().item()
        total_loss += loss.item() * num_b_atoms
        total_atoms += num_b_atoms

    return total_loss / total_atoms

def evaluate(model, loader):
    model.eval()
    all_preds = []
    all_targets = []

    with torch.no_grad():
        for batch in loader:
            batch = batch.to(DEVICE)
            ml_feats = getattr(batch, 'ml_global_features', None)

            preds = model(
                batch.x,
                batch.edge_index,
                batch.edge_attr,
                batch.solvent_id.squeeze(),
                batch.mask_b,
                batch.batch,
                ml_global_features=ml_feats
            )
            b_batch_indices = batch.batch[batch.mask_b]
            targets = batch.y_b

            unique_mols = torch.unique(b_batch_indices)
            for mol_idx in unique_mols:
                mask = (b_batch_indices == mol_idx)
                p_sorted, _ = torch.sort(preds[mask])
                t_sorted, _ = torch.sort(targets[mask])
                all_preds.append(p_sorted)
                all_targets.append(t_sorted)

    if not all_preds:
        return float('inf'), float('inf'), 0.0

    final_preds = torch.cat(all_preds)
    final_targets = torch.cat(all_targets)

    mae = F.l1_loss(final_preds, final_targets).item()
    mse = F.mse_loss(final_preds, final_targets).item()
    rmse = np.sqrt(mse)

    ss_res = torch.sum((final_targets - final_preds) ** 2)
    ss_tot = torch.sum((final_targets - torch.mean(final_targets)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot.item() != 0 else 0.0
    if isinstance(r2, torch.Tensor):
        r2 = r2.item()

    return mae, rmse, r2

# =============================================================================
# Optuna objective function
# =============================================================================
def objective(trial):
    """
    Optuna optimization objective function.
    Returns: average MAE from 5-fold cross-validation.
    """
    # 1. Hyperparameter search space
    batch_size = trial.suggest_categorical('batch_size', [8, 16, 32])
    lr = trial.suggest_float('lr', 1e-4, 1e-3, log=True)
    hidden_dim = trial.suggest_categorical('hidden_dim', [128, 256, 512])
    dropout = trial.suggest_float('dropout', 0.0, 0.2)
    solvent_dim = trial.suggest_categorical('solvent_dim', [32, 64, 128])
    ml_hidden_dim = trial.suggest_categorical('ml_hidden_dim', [32, 64, 128])

    # 2. Load data
    dataset = torch.load(DATA_PATH, weights_only=False)
    dev_dataset = dataset[:int(len(dataset)*0.8)]  # use only dev set for tuning

    # 3. Prepare K-Fold
    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    dev_dataset_list = list(dev_dataset)

    fold_maes = []
    node_dim, edge_dim = get_dims()

    for fold, (train_idx, val_idx) in enumerate(kf.split(dev_dataset_list)):
        train_fold = [dev_dataset_list[i] for i in train_idx]
        val_fold = [dev_dataset_list[i] for i in val_idx]

        train_loader = DataLoader(train_fold, batch_size=batch_size, shuffle=True,
                                  collate_fn=custom_collate_fn, drop_last=True)
        val_loader = DataLoader(val_fold, batch_size=batch_size, shuffle=False,
                                collate_fn=custom_collate_fn, drop_last=True)

        # 4. Initialize model (with ML global feature fusion)
        model = BoronNMRNet_V3(
            node_in_dim=node_dim,
            edge_in_dim=edge_dim,
            num_solvents=11,
            solvent_dim=solvent_dim,
            hidden_dim=hidden_dim,
            dropout=dropout,
            ml_feature_dim=20,
            ml_hidden_dim=ml_hidden_dim
        ).to(DEVICE)

        optimizer = torch.optim.Adam(model.parameters(), lr=lr)

        # 5. Training with early stopping
        patience_counter = 0
        patience = 5
        best_metrics = {'mae': float('inf'), 'rmse': float('inf'), 'r2': 0.0}

        for epoch in range(EPOCHS_PER_TRIAL):
            train_loss = train_one_epoch(model, train_loader, optimizer)
            val_mae, val_rmse, val_r2 = evaluate(model, val_loader)

            if val_mae < best_metrics['mae']:
                best_metrics = {'mae': val_mae, 'rmse': val_rmse, 'r2': val_r2}
                patience_counter = 0
            else:
                patience_counter += 1

            # Early stopping
            if patience_counter >= patience:
                break

        fold_maes.append(best_metrics['mae'])

        # Optuna pruning (report at fold level to avoid duplicate step warnings)
        trial.report(best_metrics['mae'], fold)
        if trial.should_prune():
            raise optuna.TrialPruned()

    # 6. Return average MAE
    avg_mae = np.mean(fold_maes)
    return avg_mae

# =============================================================================
# Main Program
# =============================================================================
if __name__ == "__main__":
    print("="*80)
    print("Hyperparameter Optimization (Optuna) - BoronNMRNet_V3")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Trials: {N_TRIALS}")
    print(f"  Folds: {N_FOLDS}")
    print(f"  Epochs per trial: {EPOCHS_PER_TRIAL}")
    print(f"  Device: {DEVICE}")
    print(f"  Timeout: {TIMEOUT}s ({TIMEOUT/3600:.1f}h)")
    print("="*80 + "\n")

    # Check data
    if not os.path.exists(DATA_PATH):
        print(f"Error: {DATA_PATH} not found")
        print("Please run csv_to_pyg_v2.py first")
        exit(1)

    # Create Optuna study
    study = optuna.create_study(
        direction='minimize',
        pruner=MedianPruner(n_startup_trials=5, n_warmup_steps=10),
        study_name='boron_nmr_v3'
    )

    # Start optimization
    print("Starting hyperparameter search...\n")
    study.optimize(objective, n_trials=N_TRIALS, timeout=TIMEOUT, show_progress_bar=True)

    # =============================================================================
    # Output results
    # =============================================================================
    print("\n" + "="*80)
    print("Optimization complete!")
    print("="*80)

    print("\nBest hyperparameters:")
    print("-" * 40)
    for key, value in study.best_params.items():
        print(f"  {key:15s}: {value}")

    print(f"\nBest validation MAE: {study.best_value:.4f}")

    # Statistics
    print("\n" + "="*80)
    print("Optimization statistics")
    print("="*80)
    print(f"  Completed trials: {len(study.trials)}")
    print(f"  Pruned trials: {len([t for t in study.trials if t.state == optuna.trial.TrialState.PRUNED])}")
    print(f"  Failed trials: {len([t for t in study.trials if t.state == optuna.trial.TrialState.FAIL])}")

    # Top 5 trials
    print("\nTop 5 hyperparameter combinations:")
    print("-" * 80)
    sorted_trials = sorted(study.trials, key=lambda t: t.value if t.value is not None else float('inf'))[:5]

    for i, trial in enumerate(sorted_trials, 1):
        print(f"\n#{i} Trial {trial.number} - MAE: {trial.value:.4f}")
        for key, value in trial.params.items():
            print(f"    {key:15s}: {value}")

    # Save results
    import json
    result = {
        'best_params': study.best_params,
        'best_value': study.best_value,
        'n_trials': len(study.trials),
        'top_5_trials': [
            {
                'trial_number': t.number,
                'value': t.value,
                'params': t.params
            }
            for t in sorted_trials
        ]
    }

    with open('tuning_results_v3.json', 'w') as f:
        json.dump(result, f, indent=2)

    print(f"\nResults saved to: tuning_results_v3.json")

    # Generate recommended configuration
    print("\n" + "="*80)
    print("Recommended configuration (for train_kfold_v3.py)")
    print("="*80)
    print("Set the following hyperparameters in train_kfold_v3.py:")
    print("-" * 40)
    print(f"BATCH_SIZE = {study.best_params['batch_size']}")
    print(f"LR = {study.best_params['lr']:.6f}")
    print(f"HIDDEN_DIM = {study.best_params['hidden_dim']}")
    print(f"DROPOUT = {study.best_params['dropout']:.4f}")
    print(f"SOLVENT_DIM = {study.best_params['solvent_dim']}")
    print(f"ML_HIDDEN_DIM = {study.best_params['ml_hidden_dim']}")
    print("="*80)

    # Visualize optimization history
    try:
        import matplotlib.pyplot as plt

        # Plot 1: Optimization history
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Subplot 1: Trial value history
        trial_values = [t.value for t in study.trials if t.value is not None]
        axes[0].plot(trial_values, 'o-', alpha=0.6)
        axes[0].axhline(y=study.best_value, color='r', linestyle='--',
                       label=f'Best: {study.best_value:.4f}')
        axes[0].set_xlabel('Trial Number')
        axes[0].set_ylabel('Validation MAE (ppm)')
        axes[0].set_title('Optimization History')
        axes[0].legend()
        axes[0].grid(alpha=0.3)

        # Subplot 2: Parameter importance (simplified - variance-based)
        param_importance = {}
        for param_name in study.best_params.keys():
            values = []
            for trial in study.trials:
                if trial.value is not None and param_name in trial.params:
                    values.append((trial.params[param_name], trial.value))

            if values:
                # Simple statistics: average MAE for different parameter values
                param_importance[param_name] = np.std([v[1] for v in values])

        if param_importance:
            params = list(param_importance.keys())
            importances = list(param_importance.values())
            axes[1].barh(params, importances)
            axes[1].set_xlabel('MAE Std Dev (importance proxy)')
            axes[1].set_title('Parameter Importance')

        plt.tight_layout()
        plt.savefig('tuning_history_v3.png', dpi=300, bbox_inches='tight')
        print("\nOptimization history plot saved to: tuning_history_v3.png")

    except Exception as e:
        print(f"\nWarning: Visualization failed ({e})")

    print("\nHyperparameter optimization complete! 🎉")
