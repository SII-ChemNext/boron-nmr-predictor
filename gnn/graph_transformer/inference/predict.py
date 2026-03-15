"""
11B NMR Chemical Shift Prediction Script
Ensemble prediction using trained 5-fold cross-validation models.

Usage:
    python predict.py

Before running, modify:
    1. MOLECULE_SMILES: the SMILES of the molecule to predict
    2. SOLVENT_SMILES: choose one from the predefined solvent list below
"""

import torch
import torch.nn.functional as F
from torch_geometric.data import Data
import numpy as np
import pickle
import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from datetime import datetime

# Add parent directory to path to import modules
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_DIR)

from models.model_v3 import BoronNMRNet_V3
from features import get_atom_features, get_bond_features
from ml_features import compute_single, N_ML_FEATURES

# =============================================================================
# User configuration area - modify the molecule and solvent to predict here
# =============================================================================

# [Modify here] Enter the SMILES of the molecule to predict
MOLECULE_SMILES = "CC1(C)C(C)(C)OB(O1)[C@H]2COCC[C@@](C3=CC=C(OC)C=C3)2B4OC(C)(C)C(C)(C)O4"  # Example: molecule with 2 boron atoms

# [Modify here] Select solvent SMILES (copy from the solvent list below)
SOLVENT_SMILES = "[2H]C(Cl)(Cl)Cl"  # Example: CDCl3

# =============================================================================
# Solvent list - copy the desired solvent SMILES to SOLVENT_SMILES above
# =============================================================================
"""
Available solvents and their SMILES strings (copy and use directly):

1. CDCl3 (deuterated chloroform) - most common NMR solvent
   [2H]C(Cl)(Cl)Cl

2. C6D6 (deuterated benzene)
   [2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]

3. DMSO-d6 (deuterated dimethyl sulfoxide) - strongly polar protic solvent
   [2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]

4. Acetone-d6 (deuterated acetone)
   [2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]

5. CD3CN (deuterated acetonitrile)
   [2H]C([2H])([2H])C#N

6. CD3OD (deuterated methanol) - protic solvent
   [2H]OC([2H])([2H])[2H]

7. CD2Cl2 (deuterated dichloromethane)
   [2H]C([2H])(Cl)Cl

8. THF-d8 (deuterated tetrahydrofuran)
   [2H]C1([2H])OC([2H])([2H])C([2H])([2H])C1([2H])[2H]

9. Toluene-d8 (deuterated toluene)
   [2H]c1c([2H])c([2H])c(C([2H])([2H])[2H])c([2H])c1[2H]

10. D2O (heavy water) - polar protic solvent
    [2H]O[2H]
"""

# =============================================================================
# Model configuration (do not modify)
# =============================================================================
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
HIDDEN_DIM = 256
DROPOUT = 0.012558398103042557
SOLVENT_DIM = 32
ML_FEATURE_DIM = 20
ML_HIDDEN_DIM = 64
NUM_SOLVENTS = 11
K_FOLDS = 5

# Model file paths
MODEL_DIR = os.path.join(PROJECT_DIR, 'trained_models')
MODEL_PATHS = [os.path.join(MODEL_DIR, f'model_v3_fold_{i+1}.pth') for i in range(K_FOLDS)]
SCALER_PATH = os.path.join(PROJECT_DIR, 'ml_feature_scaler.pkl')

# Solvent ID mapping
SOLVENT_TO_ID = {
    "[2H]C(Cl)(Cl)Cl": 0,                   # CDCl3
    "[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]": 1,  # C6D6
    "[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]": 2,   # DMSO-d6
    "[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]": 3,   # Acetone-d6
    "[2H]C([2H])([2H])C#N": 4,              # CD3CN
    "[2H]OC([2H])([2H])[2H]": 5,            # CD3OD
    "[2H]C([2H])(Cl)Cl": 6,                 # CD2Cl2
    "[2H]C1([2H])OC([2H])([2H])C([2H])([2H])C1([2H])[2H]": 7,  # THF-d8
    "[2H]c1c([2H])c([2H])c(C([2H])([2H])[2H])c([2H])c1[2H]": 8, # Toluene-d8
    "[2H]O[2H]": 9                          # D2O
}

SOLVENT_NAMES = {
    0: "CDCl3", 1: "C6D6", 2: "DMSO-d6", 3: "Acetone-d6",
    4: "CD3CN", 5: "CD3OD", 6: "CD2Cl2", 7: "THF-d8",
    8: "Toluene-d8", 9: "D2O", 10: "Unknown"
}

# =============================================================================
# ML global feature computation
# =============================================================================

def load_scaler(scaler_path):
    """Load the StandardScaler and medians saved during training."""
    with open(scaler_path, 'rb') as f:
        obj = pickle.load(f)
    return obj['scaler'], obj['medians']

def compute_ml_features_single(smiles, scaler, medians):
    """
    Compute ML global features for a single molecule (consistent with training).
    Returns: torch.Tensor of shape (1, 20)
    """
    raw = compute_single(smiles)
    nan_mask = np.isnan(raw)
    raw[nan_mask] = medians[nan_mask]
    scaled = scaler.transform(raw.reshape(1, -1))
    return torch.tensor(scaled, dtype=torch.float32)

# =============================================================================
# Utility functions
# =============================================================================

def smiles_to_graph(smiles, solvent_smiles, scaler, medians):
    """
    Convert a SMILES string to a PyG Data object.

    Important: do not add explicit hydrogen atoms — consistent with training data processing.

    Returns:
        data: PyG Data object
        num_boron: number of boron atoms
        solvent_id: solvent ID
        canonical_smiles: canonicalized SMILES (for visualization)
    """
    # Parse molecule and canonicalize SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES: {smiles}")

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)

    # Re-parse canonicalized SMILES
    mol = Chem.MolFromSmiles(canonical_smiles)
    if mol is None:
        raise ValueError(f"Canonicalization failed: {smiles} -> {canonical_smiles}")

    # Compute Gasteiger charges
    try:
        AllChem.ComputeGasteigerCharges(mol)
    except:
        print("Warning: Failed to compute Gasteiger charges, using default value 0")

    # Extract node features
    node_features = []
    boron_mask = []

    for atom in mol.GetAtoms():
        node_features.append(get_atom_features(atom))
        boron_mask.append(atom.GetSymbol() == 'B')

    x = torch.stack(node_features)
    mask_b = torch.tensor(boron_mask, dtype=torch.bool)

    num_boron = mask_b.sum().item()
    if num_boron == 0:
        raise ValueError(f"No boron (B) atoms detected in molecule! SMILES: {smiles}")

    # Extract edge features
    edge_indices = []
    edge_features = []

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        edge_feat = get_bond_features(bond)

        edge_indices.append([i, j])
        edge_indices.append([j, i])
        edge_features.append(edge_feat)
        edge_features.append(edge_feat)

    if len(edge_indices) == 0:
        raise ValueError(f"Molecule has no chemical bonds! SMILES: {smiles}")

    edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
    edge_attr = torch.stack(edge_features)

    # Get solvent ID
    solvent_id = SOLVENT_TO_ID.get(solvent_smiles, 10)
    solvent_id_tensor = torch.tensor([solvent_id], dtype=torch.long)

    # Compute ML global features
    ml_global_features = compute_ml_features_single(canonical_smiles, scaler, medians)

    # Create PyG Data object
    data = Data(
        x=x,
        edge_index=edge_index,
        edge_attr=edge_attr,
        mask_b=mask_b,
        solvent_id=solvent_id_tensor,
        ml_global_features=ml_global_features
    )

    return data, num_boron, solvent_id, canonical_smiles

def load_ensemble_models(node_dim, edge_dim):
    """Load models from all 5 folds."""
    models = []

    for fold_idx, model_path in enumerate(MODEL_PATHS):
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Model file not found: {model_path}")

        model = BoronNMRNet_V3(
            node_in_dim=node_dim,
            edge_in_dim=edge_dim,
            num_solvents=NUM_SOLVENTS,
            solvent_dim=SOLVENT_DIM,
            hidden_dim=HIDDEN_DIM,
            dropout=DROPOUT,
            ml_feature_dim=ML_FEATURE_DIM,
            ml_hidden_dim=ML_HIDDEN_DIM
        ).to(DEVICE)

        try:
            model.load_state_dict(torch.load(model_path, map_location=DEVICE, weights_only=True))
        except:
            model.load_state_dict(torch.load(model_path, map_location=DEVICE))

        model.eval()
        models.append(model)
        print(f"  [√] Fold {fold_idx + 1} model loaded successfully")

    return models

def predict_ensemble(models, data):
    """
    Perform ensemble prediction.

    Returns:
        predictions: predicted values for each boron atom (in atom index order)
        boron_indices: list of boron atom indices in the molecule
    """
    all_fold_preds = []

    data = data.to(DEVICE)
    ml_feats = getattr(data, 'ml_global_features', None)

    boron_indices = torch.where(data.mask_b)[0].cpu().tolist()

    with torch.no_grad():
        for model in models:
            preds = model(
                data.x,
                data.edge_index,
                data.edge_attr,
                data.solvent_id.squeeze(),
                data.mask_b,
                torch.zeros(data.x.size(0), dtype=torch.long, device=DEVICE),
                ml_global_features=ml_feats
            )
            # Do not sort — preserve atom index order
            all_fold_preds.append(preds)

    ensemble_pred = torch.stack(all_fold_preds).mean(dim=0)
    return ensemble_pred.cpu().numpy(), boron_indices

def visualize_molecule_with_predictions(smiles, predictions, solvent_name, boron_indices):
    """
    Visualize the molecular structure with boron atoms labeled by their predicted chemical shifts.
    """
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Warning: Cannot visualize molecule")
        return None

    if len(boron_indices) == 0:
        print("Warning: No boron atoms in molecule")
        return None

    boron_annotations = {}
    for i, atom_idx in enumerate(boron_indices):
        boron_annotations[atom_idx] = (i + 1, predictions[i])

    rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(1400, 1000)
    drawer.drawOptions().addAtomIndices = False

    for atom in mol.GetAtoms():
        if atom.GetIdx() in boron_annotations:
            boron_num, pred_value = boron_annotations[atom.GetIdx()]
            atom.SetProp("atomLabel", f"B{boron_num}")

    highlight_atoms = boron_indices
    highlight_colors = {idx: (1.0, 0.5, 0.5) for idx in boron_indices}
    highlight_radii = {idx: 0.4 for idx in boron_indices}

    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms,
        highlightBonds=[],
        highlightAtomColors=highlight_colors,
        highlightBondColors={},
        highlightAtomRadii=highlight_radii,
        confId=-1,
        legend=f""
    )

    drawer.FinishDrawing()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = os.path.join(SCRIPT_DIR, f"molecule_prediction_{timestamp}.png")

    with open(output_path, 'wb') as f:
        f.write(drawer.GetDrawingText())

    try:
        from PIL import Image, ImageDraw, ImageFont

        img = Image.open(output_path)
        draw = ImageDraw.Draw(img)

        try:
            font_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 28)
            font_label = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 24)
        except:
            font_title = ImageFont.load_default()
            font_label = ImageFont.load_default()

        y_start = 30
        title_text = f"11B NMR Chemical Shift Predictions (Solvent: {solvent_name})"
        draw.text((50, y_start), title_text, fill=(0, 0, 150), font=font_title)

        y_position = y_start + 60
        for i, atom_idx in enumerate(boron_indices):
            boron_num, pred_value = boron_annotations[atom_idx]
            label_text = f"B{boron_num}: {pred_value:.2f} ppm"

            box_size = 20
            box_left = 60
            box_top = y_position
            draw.rectangle(
                [(box_left, box_top), (box_left + box_size, box_top + box_size)],
                fill=(255, 128, 128),
                outline=(200, 0, 0),
                width=3
            )

            draw.text((box_left + box_size + 15, y_position - 5),
                     label_text, fill=(0, 0, 0), font=font_label)
            y_position += 45

        img.save(output_path)

    except Exception as e:
        print(f"   [!] PIL processing failed: {e}")
        pass

    return output_path

# =============================================================================
# Main Program
# =============================================================================

def main():
    print("="*80)
    print(" 11B NMR Chemical Shift Prediction System (Ensemble)")
    print("="*80)
    print(f"Device: {DEVICE}")
    print(f"Number of models: {K_FOLDS} folds")
    print()

    # 0. Load scaler
    print("0. Loading ML feature scaler...")
    if not os.path.exists(SCALER_PATH):
        print(f"   [x] Scaler file not found: {SCALER_PATH}")
        return
    scaler, medians = load_scaler(SCALER_PATH)
    print(f"   [√] Scaler loaded successfully")

    # 1. Input validation
    print("\n1. Validating inputs...")
    print(f"   Molecule SMILES: {MOLECULE_SMILES}")
    print(f"   Solvent SMILES: {SOLVENT_SMILES}")

    # 2. Convert to graph data
    print("\n2. Parsing molecular structure...")
    try:
        test_mol = Chem.MolFromSmiles('CB')
        node_dim = get_atom_features(test_mol.GetAtoms()[0]).shape[0]
        edge_dim = get_bond_features(test_mol.GetBonds()[0]).shape[0]
        print(f"   Node feature dimension: {node_dim}")
        print(f"   Edge feature dimension: {edge_dim}")

        data, num_boron, solvent_id, canonical_smiles = smiles_to_graph(
            MOLECULE_SMILES, SOLVENT_SMILES, scaler, medians
        )
        print(f"   [√] Molecule parsed successfully")
        print(f"   - Number of atoms: {data.x.size(0)}")
        print(f"   - Number of boron atoms: {num_boron}")
        print(f"   - Number of bonds: {data.edge_index.size(1) // 2}")
        print(f"   - Solvent: {SOLVENT_NAMES[solvent_id]}")
        if canonical_smiles != MOLECULE_SMILES:
            print(f"   - Canonical SMILES: {canonical_smiles}")

    except Exception as e:
        print(f"   [x] Error: {e}")
        return

    # 3. Load models
    print("\n3. Loading trained models...")
    try:
        models = load_ensemble_models(node_dim, edge_dim)
        print(f"   [√] All models loaded successfully")
    except Exception as e:
        print(f"   [x] Model loading failed: {e}")
        return

    # 4. Prediction
    print("\n4. Running prediction...")
    try:
        predictions, boron_indices = predict_ensemble(models, data)
        print(f"   [√] Prediction complete")
    except Exception as e:
        print(f"   [x] Prediction failed: {e}")
        return

    # 5. Visualize molecular structure
    print("\n5. Generating molecular structure image...")
    try:
        img_path = visualize_molecule_with_predictions(
            canonical_smiles,
            predictions,
            SOLVENT_NAMES[solvent_id],
            boron_indices
        )
        if img_path:
            print(f"   [√] Molecular structure image saved: {img_path}")
        else:
            print(f"   [!] Unable to generate molecular structure image")
    except Exception as e:
        print(f"   [!] Visualization failed: {e}")

    # 6. Output results
    print("\n" + "="*80)
    print(" Prediction Results")
    print("="*80)
    print(f"\nInput molecule: {MOLECULE_SMILES}")
    print(f"Canonical SMILES: {canonical_smiles}")
    print(f"Solvent: {SOLVENT_NAMES[solvent_id]}")
    print(f"\nPredicted 11B NMR chemical shifts (ppm):")
    print("-" * 40)

    for i, (atom_idx, pred) in enumerate(zip(boron_indices, predictions), 1):
        print(f"  Boron atom B{i} (atom index {atom_idx}): {pred:.2f} ppm")

    print("\n" + "="*80)
    print("Prediction complete!")
    print("="*80)

if __name__ == "__main__":
    main()
