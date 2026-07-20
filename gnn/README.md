# GNN Models for ¹¹B NMR Chemical Shift Prediction

Graph neural network models for predicting ¹¹B NMR chemical shifts, including
2D message-passing baselines and four 3D equivariant architectures. Both the
full models and the ablation without the 20 prior ML features are provided.

---

## Directory Structure

```
gnn/
├── graph_transformer/           # Main model: Graph Transformer + ML feature fusion
│   ├── models/model_v3.py       # Model architecture (BoronNMRNet_V3)
│   ├── trained_models/          # 5-fold cross-validated model weights
│   ├── inference/predict.py     # Single-molecule prediction script
│   ├── train_kfold_v3.py        # 5-fold CV training
│   ├── tune_cv_v3.py            # Hyperparameter tuning (Optuna)
│   ├── features.py              # Atom/bond feature extraction
│   ├── ml_features.py           # ML global features (SHAP Top-20)
│   ├── ml_feature_scaler.pkl    # Pre-fitted feature scaler
│   ├── build_dataset_v3.py      # Dataset construction
│   ├── csv_to_pyg.py            # CSV to PyG dataset conversion
│   └── data.csv                 # Raw dataset
│
├── graph_transformer_no_ml/     # Ablation: Graph Transformer without ML features
│   ├── models/model_v2.py
│   └── train_kfold_v2.py
│
├── gnn_comparison/              # Benchmarks with prior ML features
│   ├── models/                  # 2D models + EGNN, PaiNN, MACE, EquiformerV2
│   ├── build_equivariant_dataset.py
│   ├── equivariant_training.py
│   ├── train_*.py               # One training entry point per architecture
│   └── GNN_Comparison_Results.md
│
├── gnn_comparison_no_ml/        # Same benchmarks without prior ML features
│   ├── models/                  # Includes four no-ML equivariant wrappers
│   ├── train_*.py
│   └── GNN_Results_Summary.md
│
└── requirements.txt
```

---

## Installation

```bash
# Install PyTorch (CUDA 12.1 recommended)
pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121

# Install remaining dependencies
pip install -r requirements.txt
```

---

## Usage

### Prepare the dataset

`processed_boron_dataset_v3.pt` is not included in this repository due to file size. Generate it from the raw data before training:

```bash
cd graph_transformer

# Step 1: Convert raw CSV to a base PyG dataset
python csv_to_pyg.py

# Step 2: Append ML global features (SHAP Top-20) to get the final dataset
python build_dataset_v3.py
```

> `features.py` and `ml_features.py` are feature definition modules used internally by the above scripts and the training pipeline — they are not meant to be run directly.

The four equivariant models require a separate 3D dataset. This command reads
the same `graph_transformer/data.csv`, reproduces the seed-42 4392/1098 split,
generates one conformer per molecule, and writes ignored `.pt` files under
`gnn_comparison/processed_equivariant/`:

```bash
cd gnn_comparison
python build_equivariant_dataset.py
```

The builder uses RDKit ETKDGv3 with deterministic seeds, MMFF94 optimization
and UFF fallback. It centers each molecule and constructs directed spatial
edges with a 5.0 Å cutoff. The training split alone is used to fit the scaler
for the 20 prior ML features.

---

### Predict ¹¹B NMR chemical shifts for a new molecule

Edit `MOLECULE_SMILES` and `SOLVENT_SMILES` in `graph_transformer/inference/predict.py`, then run:

```bash
cd graph_transformer
python inference/predict.py
```

### Retrain the model

```bash
cd graph_transformer
python train_kfold_v3.py
```

### Hyperparameter tuning

Runs 50 Optuna trials (8-hour timeout) with 5-fold CV per trial. The best hyperparameters found are already applied in `train_kfold_v3.py`.

```bash
cd graph_transformer
python tune_cv_v3.py
```

### Benchmark other GNN architectures

```bash
cd gnn_comparison
python train_gatv2.py   # or train_gcn.py / train_gine.py / train_nnconv.py / train_pna.py
```

After preparing the 3D dataset, train an equivariant model with prior ML
features on one GPU:

```bash
cd gnn_comparison
python train_egnn.py --device cuda:0
python train_painn.py --device cuda:0
python train_mace.py --device cuda:0
python train_equiformer_v2.py --device cuda:0
```

Run the corresponding no-prior-ML ablations from the sibling directory. These
classes do not access `ml_global_features`, although the shared graph files
retain that tensor so the split and geometry are identical:

```bash
cd gnn_comparison_no_ml
python train_egnn.py --device cuda:0
python train_painn.py --device cuda:0
python train_mace.py --device cuda:0
python train_equiformer_v2.py --device cuda:0
```

Every training entry point performs five folds sequentially and rejects an
existing output run directory. Use a distinct `--run-name` for each rerun to
avoid overwriting results.
