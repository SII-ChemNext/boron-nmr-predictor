# GNN Models for ¹¹B NMR Chemical Shift Prediction

Graph Neural Network models for predicting ¹¹B NMR chemical shifts, based on a Graph Transformer architecture with learnable solvent embeddings and ML feature fusion.

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
├── gnn_comparison/              # Benchmark: GCN, GATv2, GINE, NNConv, PNA (with ML features)
│   ├── models/
│   ├── train_gcn/gatv2/gine/nnconv/pna.py
│   └── GNN_Comparison_Results.md
│
├── gnn_comparison_no_ml/        # Benchmark: same architectures without ML features
│   ├── models/
│   ├── train_gcn/gatv2/gine/nnconv/pna.py
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
