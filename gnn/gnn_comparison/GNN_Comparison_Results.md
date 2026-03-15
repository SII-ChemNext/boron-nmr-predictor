# GNN Variant Model Comparison Experiment Results

**Task**: Boron NMR Chemical Shift Prediction
**Framework**: 5-Fold Cross-Validation + Hold-out Test Set Ensemble Evaluation

---

## Experiment Configuration

| Item | Value |
|------|----|
| Dataset File | `processed_boron_dataset_v3.pt` |
| Total Samples | 5490 |
| Dev Set (80%, used for 5-Fold CV) | 4392 |
| Test Set (20%, Hold-out) | 1098 |
| Random Seed | 42 |
| Learning Rate | 0.0002779315932228636 |
| Batch Size | 16 |
| Hidden Dim | 256 |
| Dropout | 0.012558398103042557 |
| Solvent Dim | 32 |
| ML Feature Dim | 20 |
| ML Hidden Dim | 64 |
| Max Epochs | 150 |
| Early Stopping Patience | 10 |

> All models adopt the Virtual Node + Solvent Embedding + ML Feature Fusion architecture

---

## Results Summary

### 5-Fold Cross-Validation Mean

| Model | Avg MAE ↓ | Avg RMSE ↓ | Avg R² ↑ |
|------|-----------|------------|----------|
| GCN | 2.6423 ± 0.1848 | 6.6615 ± 0.6297 | 0.9168 ± 0.0152 |
| GATv2 | 2.4755 ± 0.2402 | 6.4710 ± 0.6549 | 0.9214 ± 0.0148 |
| GINE | 2.4085 ± 0.2294 | 6.3268 ± 0.7939 | 0.9246 ± 0.0165 |
| NNConv | 2.5466 ± 0.2547 | 6.2933 ± 0.6652 | 0.9256 ± 0.0155 |
| PNA | 2.5442 ± 0.2810 | 6.8300 ± 1.0968 | 0.9117 ± 0.0232 |

### Hold-out Test Set Ensemble Results

| Model | MAE ↓ | RMSE ↓ | R² ↑ |
|------|-------|--------|------|
| GCN | 2.2311 | 5.7127 | 0.9356 |
| **GATv2** | **1.9711** | **5.2677** | **0.9452** |
| GINE | 2.0739 | 5.5433 | 0.9394 |
| NNConv | 2.1901 | 5.5323 | 0.9396 |
| PNA | 2.0827 | 5.3440 | 0.9436 |

