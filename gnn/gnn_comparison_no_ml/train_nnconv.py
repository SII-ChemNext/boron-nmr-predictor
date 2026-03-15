import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from torch_geometric.data import Batch
import numpy as np
import random
import os
import sys
import datetime
from rdkit import Chem
from sklearn.model_selection import KFold

# Import modules from the current directory directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models.model_nnconv import BoronNMRNet_NNConv
from features import get_atom_features, get_bond_features

# =============================================================================
# 0. Logger class
# =============================================================================
class Logger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

# =============================================================================
# 1. Configuration
# =============================================================================
MODEL_NAME = "BoronNMRNet_NNConv"
LOG_FILENAME = "nnconv_results/logs/training_log.txt"

# Hyperparameters - Trial 5 (Best performing trial with MAE: 2.1272)
BATCH_SIZE = 16
LR = 0.0002779315932228636
HIDDEN_DIM = 256
DROPOUT = 0.012558398103042557
SOLVENT_DIM = 32

# Training strategy
EPOCHS = 150
PATIENCE = 10
K_FOLDS = 5

# System configuration
DATA_PATH = 'processed_boron_dataset_v2.pt'
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
SEED = 42
OUTPUT_DIR = "nnconv_results/"

# Fix random seeds
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)

# =============================================================================
# 2. Utility functions
# =============================================================================
def get_dims():
    try:
        m = Chem.MolFromSmiles('CB')
        return get_atom_features(m.GetAtoms()[0]).shape[0], get_bond_features(m.GetBonds()[0]).shape[0]
    except:
        return 0, 0

def custom_collate_fn(batch):
    """
    Custom collate function to handle solvent_id.
    """
    return Batch.from_data_list(batch)

def compute_sorted_loss(all_preds, all_targets, batch_indices):
    """
    Node-level loss aggregation: weighted by boron atom count to ensure
    consistency between the training objective and evaluation metrics.
    """
    unique_mols = torch.unique(batch_indices)
    total_loss = 0
    total_atoms = 0  # count total boron atoms, not molecules

    for mol_idx in unique_mols:
        mask = (batch_indices == mol_idx)
        p, t = all_preds[mask], all_targets[mask]
        p_sorted, _ = torch.sort(p)
        t_sorted, _ = torch.sort(t)

        # reduction='sum' gives the molecule-level total loss
        mol_loss = F.huber_loss(p_sorted, t_sorted, delta=1.0, reduction='sum')

        total_loss += mol_loss
        total_atoms += len(p)  # accumulate boron atom count

    # divide by total boron atoms -> node-level average loss
    return total_loss / total_atoms if total_atoms > 0 else total_loss

def train_one_epoch(model, loader, optimizer):
    model.train()
    total_loss = 0
    total_atoms = 0  # count atoms, not molecules

    for batch in loader:
        batch = batch.to(DEVICE)
        optimizer.zero_grad()

        # Call the updated forward method
        preds = model(
            batch.x,
            batch.edge_index,
            batch.edge_attr,
            batch.solvent_id.squeeze(),  # [batch_size]
            batch.mask_b,
            batch.batch
        )

        b_batch_indices = batch.batch[batch.mask_b]
        loss = compute_sorted_loss(preds, batch.y_b, b_batch_indices)

        loss.backward()
        optimizer.step()

        # Accumulate loss and atom count
        num_b_atoms = batch.mask_b.sum().item()
        total_loss += loss.item() * num_b_atoms
        total_atoms += num_b_atoms

    return total_loss / total_atoms  # node-level average loss

def evaluate(model, loader):
    model.eval()
    all_preds = []
    all_targets = []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(DEVICE)
            preds = model(
                batch.x,
                batch.edge_index,
                batch.edge_attr,
                batch.solvent_id.squeeze(),
                batch.mask_b,
                batch.batch
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
        return 0.0, 0.0, 0.0

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
# 3. Report generator
# =============================================================================
def print_experiment_header(train_len, test_len, total_len):
    print("="*80)
    print(f" EXPERIMENT REPORT: {MODEL_NAME}")
    print("="*80)

    print("\n1. Model Identity")
    print("-" * 40)
    print(f"   Model Name      : {MODEL_NAME}")
    print(f"   Script Run Time : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"   Description     : 5-Fold CV with Virtual Node + Learnable Solvent Embedding")

    print("\n2. Data Information")
    print("-" * 40)
    print(f"   Dataset File    : {DATA_PATH}")
    print(f"   Total Samples   : {total_len}")
    print(f"   Dev Set (80%)   : {train_len} (Used for 5-Fold CV)")
    print(f"   Test Set (20%)  : {test_len}  (Hold-out for Final Evaluation)")
    print(f"   Random Seed     : {SEED}")

    print("\n3. Configuration Details")
    print("-" * 40)
    print(f"   Learning Rate   : {LR}")
    print(f"   Batch Size      : {BATCH_SIZE}")
    print(f"   Hidden Dim      : {HIDDEN_DIM}")
    print(f"   Dropout         : {DROPOUT}")
    print(f"   Solvent Dim     : {SOLVENT_DIM}")
    print(f"   Max Epochs      : {EPOCHS}")
    print(f"   Patience        : {PATIENCE}")
    print(f"   Device          : {DEVICE}")
    print("="*80 + "\n")

# =============================================================================
# 4. Main Program
# =============================================================================
if __name__ == "__main__":
    # Redirect stdout
    sys.stdout = Logger(LOG_FILENAME)

    if not os.path.exists(DATA_PATH):
        print(f"Error: {DATA_PATH} not found")
        print("Please run csv_to_pyg_v2.py first to generate the dataset")
        exit()

    # Load and split data
    dataset = torch.load(DATA_PATH)
    random.shuffle(dataset)

    total_len = len(dataset)
    test_split_idx = int(total_len * 0.8)
    dev_dataset = dataset[:test_split_idx]
    final_test_dataset = dataset[test_split_idx:]

    # Print sections 1-3 (Header)
    print_experiment_header(len(dev_dataset), len(final_test_dataset), total_len)

    print("4. Training Process & Results")
    print("-" * 60)

    node_dim, edge_dim = get_dims()

    # Prepare 5-Fold
    kf = KFold(n_splits=K_FOLDS, shuffle=True, random_state=SEED)
    dev_dataset_list = list(dev_dataset)

    fold_metrics = {'mae': [], 'rmse': [], 'r2': []}

    for fold, (train_idx, val_idx) in enumerate(kf.split(dev_dataset_list)):
        print(f"\n>>> Fold {fold + 1} / {K_FOLDS}")

        train_fold = [dev_dataset_list[i] for i in train_idx]
        val_fold = [dev_dataset_list[i] for i in val_idx]

        train_loader = DataLoader(train_fold, batch_size=BATCH_SIZE, shuffle=True,
                                  collate_fn=custom_collate_fn, drop_last=True)
        val_loader = DataLoader(val_fold, batch_size=BATCH_SIZE, shuffle=False,
                                collate_fn=custom_collate_fn, drop_last=True)

        # Initialize model
        model = BoronNMRNet_NNConv(
            node_in_dim=node_dim,
            edge_in_dim=edge_dim,
            num_solvents=11,      # 10 solvents + 1 unknown
            solvent_dim=SOLVENT_DIM,
            hidden_dim=HIDDEN_DIM,
            dropout=DROPOUT
        ).to(DEVICE)

        optimizer = torch.optim.Adam(model.parameters(), lr=LR)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='min', factor=0.5, patience=PATIENCE
        )

        best_res = {'mae': float('inf'), 'rmse': 0, 'r2': 0}

        for epoch in range(1, EPOCHS + 1):
            train_loss = train_one_epoch(model, train_loader, optimizer)
            val_mae, val_rmse, val_r2 = evaluate(model, val_loader)
            scheduler.step(val_mae)

            if val_mae < best_res['mae']:
                best_res = {'mae': val_mae, 'rmse': val_rmse, 'r2': val_r2}
                torch.save(model.state_dict(), f"{OUTPUT_DIR}model_fold_{fold+1}.pth")

            if epoch % 10 == 0:
                print(f"   Epoch {epoch:03d} | Loss: {train_loss:.4f} | Val MAE: {val_mae:.4f}")

        print(f"   [Fold {fold+1} BEST] MAE: {best_res['mae']:.4f} | RMSE: {best_res['rmse']:.4f} | R2: {best_res['r2']:.4f}")
        for k in fold_metrics:
            fold_metrics[k].append(best_res[k])

    print("\n" + "="*60)
    print(f"Cross-Validation Summary (Average of {K_FOLDS} folds):")
    print(f"   Avg MAE  : {np.mean(fold_metrics['mae']):.4f} ± {np.std(fold_metrics['mae']):.4f}")
    print(f"   Avg RMSE : {np.mean(fold_metrics['rmse']):.4f} ± {np.std(fold_metrics['rmse']):.4f}")
    print(f"   Avg R2   : {np.mean(fold_metrics['r2']):.4f} ± {np.std(fold_metrics['r2']):.4f}")
    print("="*60)

    # Ensemble
    print("\n>>> Final Ensemble Evaluation (on Hold-out Test Set)...")

    test_loader = DataLoader(final_test_dataset, batch_size=1, shuffle=False,
                             collate_fn=custom_collate_fn)
    models = []

    for fold in range(K_FOLDS):
        m = BoronNMRNet_NNConv(
            node_in_dim=node_dim,
            edge_in_dim=edge_dim,
            num_solvents=11,
            solvent_dim=SOLVENT_DIM,
            hidden_dim=HIDDEN_DIM,
            dropout=DROPOUT
        ).to(DEVICE)

        try:
            m.load_state_dict(torch.load(f"{OUTPUT_DIR}model_fold_{fold+1}.pth",
                                         map_location=DEVICE, weights_only=True))
        except:
            m.load_state_dict(torch.load(f"{OUTPUT_DIR}model_fold_{fold+1}.pth",
                                         map_location=DEVICE))
        m.eval()
        models.append(m)

    all_preds, all_truths = [], []
    with torch.no_grad():
        for batch in test_loader:
            batch = batch.to(DEVICE)
            targets = batch.y_b
            t_sorted, _ = torch.sort(targets)

            fold_preds = []
            for m in models:
                p = m(batch.x, batch.edge_index, batch.edge_attr,
                      batch.solvent_id.squeeze(), batch.mask_b, batch.batch)
                p_sorted, _ = torch.sort(p)
                fold_preds.append(p_sorted)

            all_preds.append(torch.stack(fold_preds).mean(dim=0))
            all_truths.append(t_sorted)

    final_preds = torch.cat(all_preds)
    final_truths = torch.cat(all_truths)

    final_mae = F.l1_loss(final_preds, final_truths).item()
    final_rmse = np.sqrt(F.mse_loss(final_preds, final_truths).item())

    ss_res = torch.sum((final_truths - final_preds) ** 2)
    ss_tot = torch.sum((final_truths - torch.mean(final_truths)) ** 2)
    final_r2 = 1 - (ss_res / ss_tot) if ss_tot.item() != 0 else 0.0

    print(f"\n[FINAL ENSEMBLE RESULT]")
    print(f"   MAE  : {final_mae:.4f}")
    print(f"   RMSE : {final_rmse:.4f}")
    print(f"   R2   : {final_r2:.4f}")
    print("="*80)
