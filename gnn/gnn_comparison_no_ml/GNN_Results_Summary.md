# GNN Variant Model Experiment Results Summary (Without ML Features)

**Task**: Boron NMR Chemical Shift Prediction
**Framework**: 5-Fold Cross-Validation + Hold-out Test Set Ensemble Evaluation

---

## Experiment Configuration

| Parameter | Value |
|------|----|
| Raw Dataset | `../graph_transformer/data.csv` |
| Total Samples | 5490 |
| Dev Set (80%, used for 5-Fold CV) | 4392 |
| Test Set (20%, Hold-out) | 1098 |
| Random Seed | 42 |

The 2D models use the virtual-node dataset without prior ML features. The four
3D models reuse `../gnn_comparison/processed_equivariant`; the 20-dimensional
tensor may be present in each graph for data sharing, but the no-ML model classes
never read it and contain no descriptor encoder or descriptor-fusion branch.

### Model-family configuration

| Family | Graph / geometry | Main training settings |
|------|-------------------|------------------------|
| GCN, GATv2, GINE, NNConv, PNA | 2D covalent graph + virtual node | batch 16, lr 0.00027793, max 150 epochs, patience 10 |
| EGNN, PaiNN, MACE, EquiformerV2 | one optimized 3D conformer, 5.0 Å radius graph, 32 Gaussian RBFs | batch 16, lr 0.0003, max 150 epochs, patience 30 |

All models retain the learnable 32-dimensional solvent embedding.

---

## Results Summary

### 5-Fold Cross-Validation Mean

| Model | Avg MAE ↓ | Avg RMSE ↓ | Avg R² ↑ |
|------|-----------|------------|----------|
| GCN | 2.7023 ± 0.1899 | 6.8933 ± 0.6157 | 0.9111 ± 0.0141 |
| GATv2 | 2.5163 ± 0.2306 | 6.5734 ± 0.6444 | 0.9187 ± 0.0160 |
| GINE | 2.3662 ± 0.2469 | 6.2665 ± 0.8559 | 0.9255 ± 0.0190 |
| NNConv | 2.6286 ± 0.1784 | 6.6858 ± 0.2777 | 0.9170 ± 0.0020 |
| PNA | 2.5502 ± 0.2735 | 6.7460 ± 0.8191 | 0.9148 ± 0.0162 |
| EGNN | 2.8833 ± 0.3425 | 6.9052 ± 1.3767 | 0.9081 ± 0.0293 |
| PaiNN | 2.5667 ± 0.2665 | 6.6412 ± 1.0029 | 0.9161 ± 0.0212 |
| MACE | 2.9686 ± 0.2642 | 7.1242 ± 0.7178 | 0.9049 ± 0.0146 |
| EquiformerV2 | 2.6035 ± 0.2426 | 7.0039 ± 0.7251 | 0.9077 ± 0.0157 |

### Hold-out Test Set Ensemble Results

| Model | MAE ↓ | RMSE ↓ | R² ↑ |
|------|-------|--------|------|
| GCN | 2.2908 | 5.8864 | 0.9316 |
| GATv2 | 2.0423 | 5.3603 | 0.9433 |
| **GINE** | **2.0419** | 5.3744 | 0.9430 |
| NNConv | 2.1876 | 5.5032 | 0.9402 |
| **PNA** | 2.0896 | **5.3211** | **0.9441** |
| EGNN | 2.3732 | 5.7919 | 0.9338 |
| PaiNN | 2.1019 | 5.6313 | 0.9374 |
| MACE | 2.3714 | 5.9641 | 0.9298 |
| EquiformerV2 | 2.0823 | 5.6649 | 0.9367 |

> **Best Model**: PNA achieves the best RMSE (5.3211) and highest R² (0.9441) on the Hold-out test set

### 3D architecture details

| Model | Main architecture hyperparameters | Trainable parameters |
|------|-----------------------------------|---------------------:|
| EGNN | 4 layers, hidden dim 128, coordinate updates disabled | 533,473 |
| PaiNN | 4 scalar-vector blocks, hidden dim 128 | 1,095,649 |
| MACE | 2 interactions, 32 channels, `lmax=2`, correlation 3 | 4,528,641 |
| EquiformerV2 | 4 blocks, 4 heads, 32 sphere channels, `lmax=2` | 6,232,721 |
