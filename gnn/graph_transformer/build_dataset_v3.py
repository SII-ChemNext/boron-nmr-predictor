"""
Build Dataset V3: Add ML Global Features on Top of V2 Dataset
===================================================
Load processed_boron_dataset_v2.pt,
compute Top-20 SHAP ML features (median imputation + standardization) for each molecule,
and save as processed_boron_dataset_v3.pt.

Usage:
    python build_dataset_v3.py
"""

import torch
import numpy as np
from ml_features import compute_batch, N_ML_FEATURES, get_feature_names

# =============================================================================
# Configuration
# =============================================================================
INPUT_PATH = 'processed_boron_dataset_v2.pt'
OUTPUT_PATH = 'processed_boron_dataset_v3.pt'
SCALER_SAVE_DIR = '.'  # Directory to save scaler and medians

# =============================================================================
# Main Program
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("Build Dataset V3: Adding Top-20 ML Global Features")
    print("=" * 70)

    # 1. Load V2 dataset
    print(f"\n[1/4] Loading dataset: {INPUT_PATH}")
    dataset = torch.load(INPUT_PATH, weights_only=False)
    print(f"  Number of samples: {len(dataset)}")

    # 2. Extract SMILES
    print(f"\n[2/4] Extracting SMILES...")
    smiles_list = []
    missing_smiles_indices = []

    for i, data in enumerate(dataset):
        if hasattr(data, 'smiles') and data.smiles:
            smiles_list.append(data.smiles)
        else:
            # If v2 data lacks a smiles attribute, record it and use a placeholder
            missing_smiles_indices.append(i)
            smiles_list.append('')

    if missing_smiles_indices:
        print(f"  [WARNING] {len(missing_smiles_indices)} samples are missing the SMILES attribute!")
        print(f"  ML features for these samples will be filled with median values.")
    else:
        print(f"  Successfully extracted {len(smiles_list)} SMILES")

    # 3. Compute ML features
    print(f"\n[3/4] Computing Top-20 ML features...")
    scaled_features, scaler, medians = compute_batch(
        smiles_list, save_dir=SCALER_SAVE_DIR
    )

    print(f"\n  Feature matrix shape: {scaled_features.shape}")
    print(f"  Feature list:")
    for idx, name in enumerate(get_feature_names()):
        print(f"    [{idx:2d}] {name}")

    # 4. Add to Data objects and save
    print(f"\n[4/4] Adding ml_global_features to each Data object...")
    for i, data in enumerate(dataset):
        # Store as a 2D tensor of shape (1, 20)
        # PyG Batch will concatenate along dim=0 -> (batch_size, 20)
        data.ml_global_features = torch.tensor(
            scaled_features[i], dtype=torch.float
        ).unsqueeze(0)  # shape: (1, 20)

    # Verification
    sample = dataset[0]
    print(f"\n  Verification (first sample):")
    print(f"    x shape:                  {sample.x.shape}")
    print(f"    edge_index shape:         {sample.edge_index.shape}")
    print(f"    ml_global_features shape: {sample.ml_global_features.shape}")
    print(f"    ml_global_features:       {sample.ml_global_features}")

    # Save
    torch.save(dataset, OUTPUT_PATH)
    print(f"\n{'=' * 70}")
    print(f"V3 dataset saved to: {OUTPUT_PATH}")
    print(f"Scaler saved to:     {SCALER_SAVE_DIR}/ml_feature_scaler.pkl")
    print(f"{'=' * 70}")

    # File size statistics
    import os
    v2_size = os.path.getsize(INPUT_PATH) / (1024 * 1024)
    v3_size = os.path.getsize(OUTPUT_PATH) / (1024 * 1024)
    print(f"\n  V2 file size: {v2_size:.1f} MB")
    print(f"  V3 file size: {v3_size:.1f} MB")
    print(f"  Added: {v3_size - v2_size:.1f} MB")
