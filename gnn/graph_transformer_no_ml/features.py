import torch
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdchem
from torch_geometric.data import Data

# ==========================================
# 1. Global Configuration: Atom List and Fingerprint Parameters
# ==========================================

# Allowed atom list updated according to dataset statistics
# Includes: common organic non-metals + halogens + metals/metalloids present in the data
ALLOWED_ATOMS = [
    # --- Core organic elements & halogens ---
    'B', 'C', 'N', 'O', 'H', 'F', 'Cl', 'Br', 'I', 'P', 'S', 'Si',
    # --- Other elements from the statistics plot (alphabetical order) ---
    'Al', 'As', 'Au', 'Bi', 'Cd', 'Co', 'Cr', 'Cs', 'Cu',
    'Fe', 'Ga', 'Ge', 'In', 'K', 'Li', 'Mg', 'Mn', 'Na',
    'Pb', 'Re', 'Sb', 'Se', 'Sn', 'Te', 'W'
]

# Solvent feature extraction parameters (Morgan Fingerprint)
SOLVENT_FP_SIZE = 1024
SOLVENT_RADIUS = 2  # Equivalent to ECFP4

# ==========================================
# 2. Utility Functions
# ==========================================

def one_hot_encoding(x, permitted_list):
    """
    One-hot encoding helper function.
    Converts input x into a [0, 1, 0, ...] vector.
    If x is not in the list, it is mapped to the last element (Unknown).
    """
    if x not in permitted_list:
        x = permitted_list[-1]
    return [int(x == possible) for possible in permitted_list]

# ==========================================
# 3. Core Feature Extraction Functions (with detailed comments)
# ==========================================

def get_atom_features(atom):
    """
    Extract the feature vector of a single atom.

    Input: rdkit.Chem.rdchem.Atom object
    Output: torch.Tensor (1D vector)
    """
    features = []

    # ---------------------------------------------------------
    # Feature 1: Atom Type (Atom Symbol)
    # [Chemical meaning]: Determines the core electronic structure, electronegativity,
    # and electron configuration. This is the most fundamental feature,
    # distinguishing boron (B), carbon (C), and oxygen (O).
    # ---------------------------------------------------------
    features += one_hot_encoding(atom.GetSymbol(), ALLOWED_ATOMS)

    # ---------------------------------------------------------
    # Feature 2: Atom Degree (Connectivity)
    # [Range]: 0 to 5+ (One-Hot)
    # [Chemical meaning]: Describes steric hindrance.
    # For example, primary, secondary, and tertiary carbons have different degrees,
    # leading to very different electronic environments.
    # ---------------------------------------------------------
    features += one_hot_encoding(atom.GetDegree(), [0, 1, 2, 3, 4, 5])

    # ---------------------------------------------------------
    # Feature 3: Total Number of Attached Hydrogens (Total Num Hs)
    # [Range]: 0 to 4+ (One-Hot)
    # [Chemical meaning]: Critically important for NMR.
    # Hydrogen atoms influence the shielding constant of the central atom
    # through hyperconjugation or direct electronic effects.
    # ---------------------------------------------------------
    features += one_hot_encoding(atom.GetTotalNumHs(), [0, 1, 2, 3, 4])

    # ---------------------------------------------------------
    # Feature 4: Hybridization
    # [Range]: sp, sp2, sp3, sp3d, sp3d2 (One-Hot)
    # [Chemical meaning]: Determines the s-character of the orbital.
    # Higher s-character means electrons are closer to the nucleus,
    # leading to stronger deshielding.
    # For example: sp2 boron (planar) and sp3 boron (tetrahedral) typically
    # differ by tens of ppm in chemical shift.
    # ---------------------------------------------------------
    features += one_hot_encoding(
        atom.GetHybridization(),
        [
            rdchem.HybridizationType.SP,
            rdchem.HybridizationType.SP2,
            rdchem.HybridizationType.SP3,
            rdchem.HybridizationType.SP3D,
            rdchem.HybridizationType.SP3D2
        ]
    )

    # ---------------------------------------------------------
    # Feature 5: Aromaticity (Is Aromatic)
    # [Type]: Boolean (0 or 1)
    # [Chemical meaning]: Aromatic rings generate a "ring current effect".
    # If a boron atom is above or bonded to a benzene ring, its shift
    # is significantly affected by anisotropic effects.
    # ---------------------------------------------------------
    features.append(1 if atom.GetIsAromatic() else 0)

    # ---------------------------------------------------------
    # Feature 6: Formal Charge
    # [Type]: Integer value (e.g., -1, 0, +1)
    # [Chemical meaning]: Directly affects electron cloud density.
    # Positive charge (electron-poor) -> deshielding -> downfield shift (larger value)
    # Negative charge (electron-rich) -> shielding -> upfield shift (smaller value)
    # ---------------------------------------------------------
    features.append(atom.GetFormalCharge())

    # ---------------------------------------------------------
    # Feature 7: Is In Ring
    # [Type]: Boolean
    # [Chemical meaning]: Ring structure constrains bond angles, introduces
    # ring strain, and alters hybridization orbital composition.
    # ---------------------------------------------------------
    features.append(1 if atom.IsInRing() else 0)

    # ---------------------------------------------------------
    # Feature 8: Chirality Tag
    # [Type]: Boolean (whether there is a chiral center)
    # [Chemical meaning]: Distinguishes R/S configuration.
    # Although 2D graphs have limited stereochemical capture,
    # chiral environments in diastereomers lead to different shifts.
    # ---------------------------------------------------------
    features.append(1 if atom.GetChiralTag() != rdchem.ChiralType.CHI_UNSPECIFIED else 0)


    # ---------------------------------------------------------
    # [Added] Feature 9: Gasteiger Partial Charge (electron density)
    # ---------------------------------------------------------
    try:
        # Attempt to retrieve the pre-computed charge (must be calculated on the molecule beforehand)
        # Default to 0.0 if retrieval fails
        charge = float(atom.GetProp('_GasteigerCharge'))

        # Normalization note: charges are typically between -0.5 and +0.5
        # Can be used directly, or checked for extreme values (NaN/Inf)
        if np.isnan(charge) or np.isinf(charge):
            charge = 0.0
    except:
        charge = 0.0

    features.append(charge)

    return torch.tensor(features, dtype=torch.float)

def get_bond_features(bond):
    """
    Extract the feature vector of a chemical bond.
    """
    bond_type = bond.GetBondType()

    # ---------------------------------------------------------
    # Feature 1: Bond Type (Bond Order)
    # [Range]: single, double, triple, aromatic (One-Hot)
    # [Chemical meaning]: Determines the tightness of bonding and
    # the degree of pi electron participation.
    # ---------------------------------------------------------
    features = [
        int(bond_type == rdchem.BondType.SINGLE),
        int(bond_type == rdchem.BondType.DOUBLE),
        int(bond_type == rdchem.BondType.TRIPLE),
        int(bond_type == rdchem.BondType.AROMATIC)
    ]

    # ---------------------------------------------------------
    # Feature 2: Is Conjugated
    # [Type]: Boolean
    # [Chemical meaning]: Conjugated systems allow electron delocalization.
    # This means distant electron-withdrawing groups can directly affect
    # the electron density of boron through the conjugated chain.
    # ---------------------------------------------------------
    features.append(int(bond.GetIsConjugated()))

    # ---------------------------------------------------------
    # Feature 3: Is In Ring
    # [Type]: Boolean
    # ---------------------------------------------------------
    features.append(int(bond.IsInRing()))

    # ---------------------------------------------------------
    # Feature 4: Stereo Configuration
    # [Type]: One-Hot (none, Z/Cis, E/Trans, any)
    # [Chemical meaning]: Cis/Trans isomerism has a large effect on the
    # spatial environment. For example, cis structures may cause
    # steric hindrance and change chemical shifts.
    # ---------------------------------------------------------
    stereo = bond.GetStereo()
    features.append(1 if stereo == rdchem.BondStereo.STEREONONE else 0)
    features.append(1 if stereo == rdchem.BondStereo.STEREOANY else 0)
    features.append(1 if stereo == rdchem.BondStereo.STEREOZ else 0)
    features.append(1 if stereo == rdchem.BondStereo.STEREOE else 0)

    return torch.tensor(features, dtype=torch.float)

def get_solvent_features(solvent_smiles):
    """
    Extract solvent features: Morgan Fingerprint (ECFP)
    """
    mol = Chem.MolFromSmiles(solvent_smiles)
    if mol is None:
        # If parsing fails, return a zero vector (treated as unknown solvent)
        return torch.zeros(SOLVENT_FP_SIZE, dtype=torch.float)

    # Generate fingerprint (Radius=2, 1024 bits)
    # This captures functional group information of the solvent
    # (e.g., whether it has hydroxyl, benzene ring, chlorine, etc.)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=SOLVENT_RADIUS, nBits=SOLVENT_FP_SIZE)

    array = np.zeros((0,), dtype=np.float32)
    AllChem.DataStructs.ConvertToNumpyArray(fp, array)

    return torch.tensor(array, dtype=torch.float)

# ==========================================
# 4. Main Processing Entry Point (Process Entry)
# ==========================================

def process_entry(row):
    """
    Process a single data row, converting a Pandas Row into a PyG Data object.
    """
    mol_smiles = row['Smiles']          # Molecule SMILES
    solv_smiles = row['solvent_smiles'] # Solvent SMILES
    ppm_str = str(row['ppm_values'])    # Chemical shift value string

    # --- A. Build the molecular graph ---
    mol = Chem.MolFromSmiles(mol_smiles)
    if mol is None: return None

    atom_feats = []
    boron_mask = []

    for atom in mol.GetAtoms():
        # Extract detailed atom features
        atom_feats.append(get_atom_features(atom))
        # Mark whether this atom is boron
        boron_mask.append(atom.GetSymbol() == 'B')

    x = torch.stack(atom_feats)
    mask_b = torch.tensor(boron_mask, dtype=torch.bool)

    # Validation: must contain at least one boron atom
    num_b = mask_b.sum().item()
    if num_b == 0: return None

    # Build edges (Edge Index & Features)
    edge_indices = []
    edge_attrs = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        # Undirected graph: add both (i, j) and (j, i)
        edge_indices.append([i, j])
        edge_indices.append([j, i])

        # Extract detailed bond features
        b_feat = get_bond_features(bond)
        edge_attrs.append(b_feat)
        edge_attrs.append(b_feat)

    if not edge_indices:
        # Handle single-atom molecules (rare)
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr_dim = len(get_bond_features(mol.GetBonds()[0])) if mol.GetNumBonds() > 0 else 10  # estimated dimension
        edge_attr = torch.empty((0, edge_attr_dim), dtype=torch.float)
    else:
        edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
        edge_attr = torch.stack(edge_attrs)

    # --- B. Extract solvent features ---
    # Output shape: [1, 1024]
    solvent_feat = get_solvent_features(solv_smiles).unsqueeze(0)

    # --- C. Parse labels ---
    try:
        shifts = [float(s) for s in ppm_str.replace(';', ',').split(',') if s.strip()]
    except:
        return None

    # Handle "multiple borons, single value" (symmetric) and "multiple borons, multiple values" (asymmetric)
    if len(shifts) == 1 and num_b > 1:
        shifts = shifts * num_b  # broadcast
    elif len(shifts) != num_b:
        # Mismatch, skip
        return None

    y_b = torch.tensor(shifts, dtype=torch.float)

    # --- D. Pack Data object ---
    data = Data(
        x=x,                 # Node features [Num_Nodes, Feature_Dim]
        edge_index=edge_index, # Connectivity [2, Num_Edges]
        edge_attr=edge_attr,   # Edge features [Num_Edges, Edge_Feature_Dim]
        mask_b=mask_b,       # Boron atom mask [Num_Nodes]
        y_b=y_b,             # Ground-truth chemical shifts [Num_Borons]
        solvent_x=solvent_feat # Solvent fingerprint [1, 1024]
    )

    return data

if __name__ == "__main__":
    # Assume CSV has been read into df
    # df = pd.read_csv('your_data.csv')
    # dataset = [process_entry(row) for _, row in df.iterrows() if process_entry(row) is not None]

    # Simulation test
    test_row = {
        'Smiles': 'CN1CC(=O)OB(C(B2OC(C)(C)C(C)(C)O2)=C(B2OC(C)(C)C(C)(C)O2)B2OC(C)(C)C(C)(C)O2)OC(=O)C1',
        'solvent_smiles': 'ClC(Cl)Cl', # chloroform
        'ppm_values': '28.5'
    }
    data = process_entry(test_row)

    print("X Shape:", data.x.shape)          # [Num_Atoms, 50+]
    print("Solvent Shape:", data.solvent_x.shape) # [1, 1024]
    print("Target:", data.y_b)
    # print("edge_index:", data.edge_index)
    # print("edge_attr:", data.edge_attr)
