# GNN Variant Model Experiment Results Summary (Without ML Features)

**Task**: Boron NMR Chemical Shift Prediction
**Framework**: 5-Fold Cross-Validation + Hold-out Test Set Ensemble Evaluation (Virtual Node + Learnable Solvent Embedding)

---

## Experiment Configuration

| Parameter | Value |
|------|----|
| Dataset File | `processed_boron_dataset_v2.pt` |
| Total Samples | 5490 |
| Dev Set (80%, used for 5-Fold CV) | 4392 |
| Test Set (20%, Hold-out) | 1098 |
| Random Seed | 42 |
| Learning Rate | 0.0002779315932228636 |
| Batch Size | 16 |
| Hidden Dim | 256 |
| Dropout | 0.012558398103042557 |
| Solvent Dim | 32 |
| Max Epochs | 150 |
| Early Stopping Patience | 10 |

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

### Hold-out Test Set Ensemble Results

| Model | MAE ↓ | RMSE ↓ | R² ↑ |
|------|-------|--------|------|
| GCN | 2.2908 | 5.8864 | 0.9316 |
| GATv2 | 2.0423 | 5.3603 | 0.9433 |
| **GINE** | **2.0419** | 5.3744 | 0.9430 |
| NNConv | 2.1876 | 5.5032 | 0.9402 |
| **PNA** | 2.0896 | **5.3211** | **0.9441** |

> **Best Model**: PNA achieves the best RMSE (5.3211) and highest R² (0.9441) on the Hold-out test set
