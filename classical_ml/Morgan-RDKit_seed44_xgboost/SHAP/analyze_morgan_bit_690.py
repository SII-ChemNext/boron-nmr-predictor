#!/usr/bin/env python3
"""
Morgan Bit 690 Substructure Analysis
======================================
Generate highlighted structure SVGs for molecules with and without Morgan bit 690.
"""

import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt

# Add module path
sys.path.insert(0, str(Path(__file__).parent))
from modules.utils import load_config, FeatureMapper, DataLoader


def get_morgan_fingerprint_bits(smiles: str, radius: int = 2, n_bits: int = 1024):
    """
    Get Morgan fingerprint bit information for a molecule.

    Returns:
        dict: {bit_index: [atom_indices]}
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Use bitInfo to retrieve atoms corresponding to each bit
    bit_info = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius, nBits=n_bits, bitInfo=bit_info
    )

    return bit_info, fp


def draw_molecule_with_highlights(smiles: str, bit_info: dict, target_bit: int,
                                   mol_size=(400, 400)):
    """
    Draw an SVG of a molecule with the substructure for a specific bit highlighted.

    Parameters:
        smiles: SMILES string
        bit_info: Morgan fingerprint bit information
        target_bit: Bit index to highlight
        mol_size: Image dimensions

    Returns:
        svg_string: SVG string
        highlight_atoms: List of highlighted atom indices
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, []

    # Collect atoms/bonds for target_bit
    highlight_atoms = []
    highlight_bonds = []

    if target_bit in bit_info:
        for atom_idx, radius in bit_info[target_bit]:
            if radius == 0:
                # Radius=0: center atom + direct neighbors
                highlight_atoms.append(atom_idx)
                atom = mol.GetAtomWithIdx(atom_idx)

                # Add all neighbor atoms
                for neighbor in atom.GetNeighbors():
                    highlight_atoms.append(neighbor.GetIdx())

                # Add all bonds connected to the center atom
                for bond in mol.GetBonds():
                    if bond.GetBeginAtomIdx() == atom_idx or bond.GetEndAtomIdx() == atom_idx:
                        highlight_bonds.append(bond.GetIdx())
            else:
                # Radius > 0: use standard environment method
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)
                highlight_bonds.extend(env)

                atoms_in_env = set([atom_idx])
                for bond_idx in env:
                    bond = mol.GetBondWithIdx(bond_idx)
                    atoms_in_env.add(bond.GetBeginAtomIdx())
                    atoms_in_env.add(bond.GetEndAtomIdx())
                highlight_atoms.extend(atoms_in_env)

    highlight_atoms = list(set(highlight_atoms))
    highlight_bonds = list(set(highlight_bonds))

    # Set up drawer
    drawer = rdMolDraw2D.MolDraw2DSVG(mol_size[0], mol_size[1])
    drawer.drawOptions().addAtomIndices = False
    drawer.drawOptions().bondLineWidth = 2

    # Highlight colors — vivid red for Morgan bit 690 substructure
    atom_colors = {}
    bond_colors = {}

    for atom_idx in highlight_atoms:
        atom_colors[atom_idx] = (1.0, 0.4, 0.4)  # bright red

    for bond_idx in highlight_bonds:
        bond_colors[bond_idx] = (0.9, 0.1, 0.1)  # dark red

    # Draw molecule
    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms,
        highlightBonds=highlight_bonds,
        highlightAtomColors=atom_colors,
        highlightBondColors=bond_colors
    )
    drawer.FinishDrawing()

    svg = drawer.GetDrawingText()

    return svg, highlight_atoms


def find_boron_atoms(smiles: str):
    """Return the atom indices of boron atoms in a molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    boron_indices = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'B':
            boron_indices.append(atom.GetIdx())

    return boron_indices


def main():
    print("="*80)
    print(" Morgan Bit 690 Substructure Analysis")
    print("="*80)

    # 1. Load configuration and data
    print("\n[1/6] Loading configuration and data...")
    config = load_config()
    data_loader = DataLoader(config['paths'])

    # Load test set
    X_test, y_test, smiles_df = data_loader.load_all_data('test')
    smiles_test = smiles_df['Smiles']  # Extract SMILES column
    print(f"✓ Test set loaded: {len(X_test)} samples")

    # 2. Load SHAP values
    print("\n[2/6] Loading SHAP values...")
    shap_dir = Path(__file__).parent / 'outputs/shap_global'
    feature_importance_df = pd.read_csv(shap_dir / 'feature_importance_ranking.csv')

    # Find global importance of morgan_bit_690
    bit_690_row = feature_importance_df[
        feature_importance_df['feature_name'] == 'compound_morgan_bit_690'
    ]
    if len(bit_690_row) > 0:
        bit_690_importance = bit_690_row.iloc[0]['mean_abs_shap']
        bit_690_rank = bit_690_row.index[0] + 1
        print(f"✓ Morgan bit 690 global importance: {bit_690_importance:.4f} (rank: #{bit_690_rank})")

    # 3. Extract Morgan bit 690 feature column
    print("\n[3/6] Extracting Morgan bit 690 feature...")
    feature_mapper = FeatureMapper(config['features'])
    compound_morgan_start = feature_mapper.blocks['compound_morgan']['start']
    bit_690_col_idx = compound_morgan_start + 690

    bit_690_values = X_test[:, bit_690_col_idx]
    has_bit_690 = bit_690_values > 0

    print(f"✓ Molecules with bit 690: {has_bit_690.sum()} / {len(X_test)}")
    print(f"  Fraction: {has_bit_690.sum() / len(X_test) * 100:.1f}%")

    # 4. Load per-sample SHAP values and select molecules
    print("\n[4/6] Loading SHAP values and selecting molecules...")

    shap_values_path = shap_dir / 'shap_values.npy'
    if shap_values_path.exists():
        shap_values = np.load(shap_values_path)
        print(f"✓ SHAP value matrix loaded: {shap_values.shape}")

        # Extract bit 690 SHAP column
        bit_690_shap = shap_values[:, bit_690_col_idx]

        # Molecules with bit 690 — sorted by |SHAP| descending
        has_690_indices = np.where(has_bit_690)[0]
        has_690_shap = np.abs(bit_690_shap[has_690_indices])
        sorted_has_690 = has_690_indices[np.argsort(has_690_shap)[::-1]]

        # Molecules without bit 690 — sorted by total |SHAP| ascending
        no_690_indices = np.where(~has_bit_690)[0]
        total_shap_abs = np.abs(shap_values[no_690_indices]).sum(axis=1)
        sorted_no_690 = no_690_indices[np.argsort(total_shap_abs)]
    else:
        print("⚠ SHAP value matrix not found; recomputing from model...")

        from modules.utils import ModelLoader
        import shap

        model_path = Path(__file__).parent / config['paths']['model_path']
        model = ModelLoader.load_xgboost_model(str(model_path.resolve()))

        # Use training set as background
        X_train, _, smiles_train_df = data_loader.load_all_data('train')
        background_size = config['shap_analysis']['background_size']
        X_background = X_train[:background_size]

        explainer = shap.TreeExplainer(
            model,
            data=X_background,
            feature_perturbation='interventional'
        )

        print(f"  Computing SHAP values for {len(X_test)} test samples using {background_size} background samples...")
        shap_values = explainer.shap_values(X_test)

        np.save(shap_values_path, shap_values)
        print(f"✓ SHAP values saved to: {shap_values_path.name}")

        bit_690_shap = shap_values[:, bit_690_col_idx]

        has_690_indices = np.where(has_bit_690)[0]
        has_690_shap = np.abs(bit_690_shap[has_690_indices])
        sorted_has_690 = has_690_indices[np.argsort(has_690_shap)[::-1]]

        no_690_indices = np.where(~has_bit_690)[0]
        total_shap_abs = np.abs(shap_values[no_690_indices]).sum(axis=1)
        sorted_no_690 = no_690_indices[np.argsort(total_shap_abs)]

    # Select top 10 with bit 690 + top 5 without (lowest total SHAP)
    selected_has_690 = sorted_has_690[:10]
    selected_no_690 = sorted_no_690[:5]

    print(f"✓ Selected:")
    print(f"  - Molecules with bit 690: {len(selected_has_690)}")
    print(f"  - Molecules without bit 690: {len(selected_no_690)}")

    # 5. Generate SVG highlighted structure images
    print("\n[5/6] Generating SVG highlighted structure images...")

    output_dir = Path(__file__).parent / 'outputs/morgan_bit_690'
    output_dir.mkdir(parents=True, exist_ok=True)

    svg_dir = output_dir / 'svg_structures'
    svg_dir.mkdir(exist_ok=True)

    summary_data = []

    # Process molecules with bit 690
    print("  Processing molecules with bit 690...")
    for i, idx in enumerate(selected_has_690):
        smiles = str(smiles_test.iloc[idx])
        ppm = float(y_test[idx])
        bit_690_shap_val = float(bit_690_shap[idx])

        bit_info, fp = get_morgan_fingerprint_bits(smiles)

        if bit_info is None:
            print(f"    ⚠ Molecule {i+1}: invalid SMILES, skipping")
            continue

        svg_string, highlight_atoms = draw_molecule_with_highlights(
            smiles, bit_info, 690
        )

        if svg_string is None:
            continue

        svg_filename = f'has_bit690_{i+1:02d}.svg'
        svg_path = svg_dir / svg_filename
        with open(svg_path, 'w') as f:
            f.write(svg_string)

        boron_atoms = find_boron_atoms(smiles)

        summary_data.append({
            'molecule_id': f'has_690_{i+1:02d}',
            'has_bit_690': True,
            'smiles': smiles,
            'ppm': ppm,
            'bit_690_shap_value': bit_690_shap_val,
            'highlight_atoms': highlight_atoms,
            'boron_atoms': boron_atoms,
            'svg_file': svg_filename
        })

        print(f"    ✓ Molecule {i+1}/10: SHAP={bit_690_shap_val:.4f}, ppm={ppm:.2f}")

    # Process molecules without bit 690
    print("  Processing molecules without bit 690...")
    for i, idx in enumerate(selected_no_690):
        smiles = str(smiles_test.iloc[idx])
        ppm = float(y_test[idx])
        total_shap_val = float(total_shap_abs[np.where(no_690_indices == idx)[0][0]])

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        # Draw plain SVG (no highlight)
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.drawOptions().addAtomIndices = False
        drawer.drawOptions().bondLineWidth = 2
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_string = drawer.GetDrawingText()

        svg_filename = f'no_bit690_{i+1:02d}.svg'
        svg_path = svg_dir / svg_filename
        with open(svg_path, 'w') as f:
            f.write(svg_string)

        boron_atoms = find_boron_atoms(smiles)

        summary_data.append({
            'molecule_id': f'no_690_{i+1:02d}',
            'has_bit_690': False,
            'smiles': smiles,
            'ppm': ppm,
            'total_shap_value': total_shap_val,
            'highlight_atoms': [],
            'boron_atoms': boron_atoms,
            'svg_file': svg_filename
        })

        print(f"    ✓ Molecule {i+1}/5: total SHAP={total_shap_val:.4f}, ppm={ppm:.2f}")

    # 6. Save summary data
    print("\n[6/6] Saving summary data...")

    summary_json_path = output_dir / 'bit_690_analysis_summary.json'
    with open(summary_json_path, 'w', encoding='utf-8') as f:
        json.dump(summary_data, f, indent=2, ensure_ascii=False)

    summary_csv_path = output_dir / 'bit_690_analysis_summary.csv'
    df_summary = pd.DataFrame([
        {
            'molecule_id': d['molecule_id'],
            'has_bit_690': d['has_bit_690'],
            'smiles': d['smiles'],
            'ppm': d['ppm'],
            'shap_value': d.get('bit_690_shap_value', d.get('total_shap_value', 0)),
            'num_highlight_atoms': len(d['highlight_atoms']),
            'boron_atoms': ','.join(map(str, d['boron_atoms'])),
            'svg_file': d['svg_file']
        }
        for d in summary_data
    ])
    df_summary.to_csv(summary_csv_path, index=False)

    detail_path = output_dir / 'bit_690_highlight_atoms_detail.txt'
    with open(detail_path, 'w', encoding='utf-8') as f:
        f.write("Morgan Bit 690 Highlighted Atom Details\n")
        f.write("="*80 + "\n\n")

        for d in summary_data:
            f.write(f"Molecule ID: {d['molecule_id']}\n")
            f.write(f"SMILES: {d['smiles']}\n")
            f.write(f"Chemical shift: {d['ppm']:.2f} ppm\n")
            f.write(f"Has bit 690: {d['has_bit_690']}\n")

            if d['has_bit_690']:
                f.write(f"Bit 690 SHAP value: {d.get('bit_690_shap_value', 0):.4f}\n")
                f.write(f"Highlighted atom indices: {d['highlight_atoms']}\n")
            else:
                f.write(f"Total SHAP value: {d.get('total_shap_value', 0):.4f}\n")

            f.write(f"Boron atom indices: {d['boron_atoms']}\n")
            f.write(f"SVG file: {d['svg_file']}\n")
            f.write("-"*80 + "\n\n")

    print(f"✓ Summary data saved:")
    print(f"  - JSON: {summary_json_path.name}")
    print(f"  - CSV: {summary_csv_path.name}")
    print(f"  - Atom detail list: {detail_path.name}")

    print("\n" + "="*80)
    print(f"✓ Analysis complete!")
    print(f"✓ SVG structure images saved to: {svg_dir}")
    print(f"✓ Summary data saved to: {output_dir}")
    print("="*80)

    print("\nGenerated files:")
    print(f"  - SVG structure images: {len(list(svg_dir.glob('*.svg')))} files")
    print(f"  - Summary JSON: bit_690_analysis_summary.json")
    print(f"  - Summary CSV: bit_690_analysis_summary.csv")
    print(f"  - Atom detail list: bit_690_highlight_atoms_detail.txt")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
