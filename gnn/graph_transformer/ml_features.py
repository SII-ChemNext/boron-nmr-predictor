"""
ML Global Features for GNN Fusion
==================================
Compute Top-20 SHAP features (from XGBoost SHAP analysis) for each molecule.
Uses the same computation method as the original Chemia ML training pipeline
(utils/mol_fp_features.py: Descriptors._descList + rdFingerprintGenerator).

Features:
  - 16 RDKit molecular descriptors (compound-level)
  - 4 Morgan fingerprint bits (radius=2, 1024 bits)

NaN handling: median imputation + StandardScaler normalization
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.preprocessing import StandardScaler
import pickle
import os

# ==========================================
# Morgan fingerprint API (consistent with Chemia)
# ==========================================
try:
    from rdkit.Chem import rdFingerprintGenerator
    HAS_NEW_FP_GENERATOR = True
except ImportError:
    from rdkit.Chem import AllChem
    HAS_NEW_FP_GENERATOR = False

# ==========================================
# Top-20 SHAP Features Definition
# (from feature_importance_ranking.csv)
# ==========================================

# 16 RDKit descriptors (ordered by SHAP importance)
TOP_RDKIT_DESCRIPTORS = [
    'BCUT2D_MRLOW',               # #1  SHAP=4.376
    'SlogP_VSA7',                 # #4  SHAP=1.294
    'MaxPartialCharge',           # #5  SHAP=1.228
    'BCUT2D_MWLOW',               # #7  SHAP=0.934
    'MinPartialCharge',           # #8  SHAP=0.664
    'MaxAbsPartialCharge',        # #9  SHAP=0.661
    'NOCount',                    # #10 SHAP=0.493
    'BCUT2D_LOGPHI',              # #11 SHAP=0.448
    'MinEStateIndex',             # #12 SHAP=0.437
    'MinAbsPartialCharge',        # #13 SHAP=0.411
    'fr_quatN',                   # #14 SHAP=0.383
    'VSA_EState8',                # #15 SHAP=0.367
    'NumAliphaticHeterocycles',   # #16 SHAP=0.340
    'PEOE_VSA3',                  # #17 SHAP=0.335
    'NumAliphaticRings',          # #19 SHAP=0.270
    'MolLogP',                    # #20 SHAP=0.270
]

# 4 Morgan fingerprint bit indices (0-indexed, radius=2, 1024 bits)
TOP_MORGAN_BITS = [690, 271, 390, 656]

# Total
N_ML_FEATURES = len(TOP_RDKIT_DESCRIPTORS) + len(TOP_MORGAN_BITS)  # 20

# Build descriptor function lookup from Descriptors._descList
# (consistent with Chemia mol_fp_features.py)
_RDKIT_DESC_FUNCS = {}
for _name, _func in Descriptors._descList:
    if _name in TOP_RDKIT_DESCRIPTORS:
        _RDKIT_DESC_FUNCS[_name] = _func

# Check that all descriptors were found
_missing = set(TOP_RDKIT_DESCRIPTORS) - set(_RDKIT_DESC_FUNCS.keys())
if _missing:
    print(f"[WARNING] The following RDKit descriptors were not found in Descriptors._descList: {_missing}")


def compute_single(smiles: str) -> np.ndarray:
    """
    Compute Top-20 ML features for a single molecule.

    Args:
        smiles: SMILES string

    Returns:
        np.ndarray of shape (20,) — raw values, may contain NaN
    """
    features = np.full(N_ML_FEATURES, np.nan)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return features

    # --- Part 1: 16 RDKit descriptors ---
    for i, desc_name in enumerate(TOP_RDKIT_DESCRIPTORS):
        func = _RDKIT_DESC_FUNCS.get(desc_name)
        if func is None:
            continue
        try:
            val = func(mol)
            if val is not None:
                val = float(val)
                if not (np.isnan(val) or np.isinf(val)):
                    features[i] = val
        except Exception:
            pass

    # --- Part 2: 4 Morgan fingerprint bits ---
    offset = len(TOP_RDKIT_DESCRIPTORS)
    try:
        if HAS_NEW_FP_GENERATOR:
            gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
            fp = gen.GetFingerprint(mol)
            fp_array = np.array(fp)
        else:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            fp_array = np.array(fp)

        for j, bit_idx in enumerate(TOP_MORGAN_BITS):
            features[offset + j] = float(fp_array[bit_idx])
    except Exception:
        pass

    return features


def compute_batch(smiles_list: list, save_dir: str = None):
    """
    Batch computation, imputation, and standardization of Top-20 ML features.

    Args:
        smiles_list: list of SMILES strings
        save_dir: if provided, save scaler and medians to this directory (for inference)

    Returns:
        scaled_features: np.ndarray (n, 20) — standardized feature matrix
        scaler: fitted StandardScaler
        medians: np.ndarray (20,) — medians used for imputation
    """
    n = len(smiles_list)
    raw = np.zeros((n, N_ML_FEATURES))

    print(f"  Computing {N_ML_FEATURES} ML features for {n} molecules...")
    for i, smi in enumerate(smiles_list):
        raw[i] = compute_single(smi)
        if (i + 1) % 1000 == 0:
            print(f"    Processed {i+1}/{n}")

    # --- NaN statistics ---
    nan_count = np.isnan(raw).sum(axis=0)
    total_nan = int(np.isnan(raw).sum())
    print(f"  NaN statistics: {total_nan} total NaN values in {n} x {N_ML_FEATURES} matrix")
    feature_names = get_feature_names()
    for idx in range(N_ML_FEATURES):
        if nan_count[idx] > 0:
            print(f"    {feature_names[idx]}: {int(nan_count[idx])} NaN ({nan_count[idx]/n*100:.1f}%)")

    # --- Median imputation ---
    medians = np.nanmedian(raw, axis=0)
    # If a feature is entirely NaN, use 0.0
    medians = np.where(np.isnan(medians), 0.0, medians)

    for j in range(N_ML_FEATURES):
        nan_mask = np.isnan(raw[:, j])
        raw[nan_mask, j] = medians[j]

    # --- StandardScaler normalization ---
    scaler = StandardScaler()
    scaled = scaler.fit_transform(raw)

    print(f"  Scaling complete.")
    print(f"    Mean range: [{scaled.mean(axis=0).min():.6f}, {scaled.mean(axis=0).max():.6f}]")
    print(f"    Std  range: [{scaled.std(axis=0).min():.4f}, {scaled.std(axis=0).max():.4f}]")

    # --- Save scaler and medians (needed for inference) ---
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, 'ml_feature_scaler.pkl')
        with open(save_path, 'wb') as f:
            pickle.dump({'scaler': scaler, 'medians': medians}, f)
        print(f"  Scaler & medians saved to: {save_path}")

    return scaled, scaler, medians


def get_feature_names():
    """Get the ordered list of feature names."""
    names = list(TOP_RDKIT_DESCRIPTORS)
    names += [f"morgan_bit_{idx}" for idx in TOP_MORGAN_BITS]
    return names


if __name__ == "__main__":
    # Simple test
    test_smiles = [
        'CC1(C)OB(C)OC1(C)C',          # pinacol boronate
        'c1ccc(B(O)O)cc1',              # phenylboronic acid
        'CB(C)C',                        # trimethylborane
        'INVALID_SMILES',                # test invalid SMILES
    ]

    print("=" * 60)
    print("ML Features Test")
    print("=" * 60)

    # Test single molecule
    for smi in test_smiles:
        feats = compute_single(smi)
        nan_count = np.isnan(feats).sum()
        print(f"\n  {smi[:40]:40s} | NaN: {nan_count}/{N_ML_FEATURES}")

    # Test batch
    print("\n" + "=" * 60)
    scaled, scaler, medians = compute_batch(test_smiles)
    print(f"\nScaled features shape: {scaled.shape}")
    print(f"Feature names: {get_feature_names()}")
