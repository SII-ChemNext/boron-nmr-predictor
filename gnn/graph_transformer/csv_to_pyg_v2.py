import torch
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem  # Required for Gasteiger charge calculation
from torch_geometric.data import Data
from tqdm import tqdm

from features import get_atom_features, get_bond_features

# =============================================================================
# Solvent mapping table - map solvent SMILES to IDs
# =============================================================================
SOLVENT_TO_ID = {
    '[2H]C(Cl)(Cl)Cl': 0,                                           # CDCl3
    '[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]': 1,                 # C6D6
    '[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]': 2,                  # DMSO-d6
    '[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]': 3,                  # Acetone-d6
    '[2H]C([2H])([2H])C#N': 4,                                     # CD3CN
    '[2H]OC([2H])([2H])[2H]': 5,                                   # CD3OD
    '[2H]C([2H])(Cl)Cl': 6,                                        # CD2Cl2
    '[2H]C1([2H])OC([2H])([2H])C([2H])([2H])C1([2H])[2H]': 7,      # THF-d8
    '[2H]c1c([2H])c([2H])c(C([2H])([2H])[2H])c([2H])c1[2H]': 8,    # Toluene-d8
    '[2H]O[2H]': 9,                                                 # D2O
}

# =============================================================================
# Configuration
# =============================================================================
CSV_PATH = 'data.csv'                  # Input CSV file name
OUTPUT_PATH = 'processed_boron_dataset_v2.pt' # Output .pt file name
ERROR_LOG_PATH = 'skipped_errors_v2.csv'      # Error log file name

# CSV column mapping
COL_SMILES = 'Smiles'
COL_SOLVENT = 'solvent'
COL_B_COUNT = 'B_count'
COL_PPM = 'ppm_values'

# =============================================================================
# Core logic: process a single row
# =============================================================================
def process_row(row):
    """
    Process a single row and return (Data, error_message)
    """
    smiles_str = row[COL_SMILES]
    solvent_str = row[COL_SOLVENT]
    ppm_str = str(row[COL_PPM])
    
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        return None, "RDKit Parsing Error (Invalid SMILES)"

    # =================================================================
    # [Added step] Compute Gasteiger partial charges (electron density)
    # This modifies the mol object by adding the '_GasteigerCharge' property to each atom
    # get_atom_features in features.py reads this property
    # =================================================================
    try:
        AllChem.ComputeGasteigerCharges(mol)
    except Exception as e:
        # In rare cases, the calculation may fail (for example, with unusual valence states), but this should not stop processing
        # If it fails, the property is absent and get_atom_features defaults to 0
        pass

    # 2. Extract node features
    atom_features_list = []
    boron_mask = []

    for atom in mol.GetAtoms():
        # Call features.py to extract features (including charge information)
        atom_features_list.append(get_atom_features(atom))
        boron_mask.append(atom.GetSymbol() == 'B')

    x = torch.stack(atom_features_list)
    mask_b = torch.tensor(boron_mask, dtype=torch.bool)
    
    # 3. Validate the number of boron atoms
    num_b_found = mask_b.sum().item()
    csv_b_count = int(row[COL_B_COUNT])
    
    if num_b_found == 0:
        return None, "No Boron Atoms Found"
        
    if num_b_found != csv_b_count:
        return None, f"Count Mismatch: CSV says {csv_b_count}, RDKit found {num_b_found}"

    # 4. Build edges (chemical bonds)
    edge_indices = []
    edge_attrs = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        
        # Represent each undirected bond with two directed edges
        edge_indices.append([i, j])
        edge_indices.append([j, i])
        
        b_feat = get_bond_features(bond)
        edge_attrs.append(b_feat)
        edge_attrs.append(b_feat)

    if len(edge_indices) == 0:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0, edge_attrs[0].size(0) if edge_attrs else 0), dtype=torch.float)
    else:
        edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
        edge_attr = torch.stack(edge_attrs)

    # 5. Map solvent SMILES to an ID; 10 represents an unknown solvent
    solvent_id = SOLVENT_TO_ID.get(solvent_str.strip(), 10)

    # 6. Process labels (PPM values)
    try:
        shifts = [float(s) for s in ppm_str.replace(';', ',').split(',') if s.strip()]
    except ValueError:
        return None, f"PPM Format Error: {ppm_str}"

    # Label alignment and broadcasting logic
    if len(shifts) == num_b_found:
        final_shifts = shifts
    elif len(shifts) == 1 and num_b_found > 1:
        # Broadcast one value to multiple equivalent B atoms
        final_shifts = shifts * num_b_found
    else:
        return None, f"Label Mismatch: Need {num_b_found} values, got {len(shifts)} ({ppm_str})"

    y_b = torch.tensor(final_shifts, dtype=torch.float)

    # 7. Pack the Data object
    data = Data(
        x=x, 
        edge_index=edge_index, 
        edge_attr=edge_attr,
        mask_b=mask_b, 
        y_b=y_b,
        solvent_id=torch.tensor([solvent_id], dtype=torch.long)
    )
    # Store raw information for debugging
    data.smiles = smiles_str 
    
    return data, None

# =============================================================================
# Main program
# =============================================================================
if __name__ == "__main__":
    print(f"1. Reading CSV file: {CSV_PATH}")
    try:
        df = pd.read_csv(CSV_PATH)
    except FileNotFoundError:
        print("Error: CSV file not found.")
        exit()
    
    data_list = []
    skipped_records = [] 
    
    print("\n2. Starting processing (computing charges + building graphs)...")
    for index, row in tqdm(df.iterrows(), total=len(df)):
        try:
            data, error_msg = process_row(row)
            
            if data is not None:
                data_list.append(data)
            else:
                skipped_records.append({
                    'Index': index,
                    'Smiles': row[COL_SMILES],
                    'Error': error_msg,
                    'Raw_PPM': row[COL_PPM]
                })
                
        except Exception as e:
            skipped_records.append({
                'Index': index,
                'Smiles': row[COL_SMILES],
                'Error': f"Unexpected Exception: {str(e)}",
                'Raw_PPM': row[COL_PPM]
            })

    print(f"\n3. Processing complete!")
    print(f"   - Successfully converted: {len(data_list)}")
    print(f"   - Skipped (invalid): {len(skipped_records)}")

    # Save valid data
    if len(data_list) > 0:
        torch.save(data_list, OUTPUT_PATH)
        print(f"   - Dataset saved to: {OUTPUT_PATH}")
        print("     (Note: this dataset uses solvent IDs rather than Morgan fingerprints)")

    # Save error log
    if len(skipped_records) > 0:
        error_df = pd.DataFrame(skipped_records)
        error_df.to_csv(ERROR_LOG_PATH, index=False)
        print(f"   - Error log saved to: {ERROR_LOG_PATH}")
