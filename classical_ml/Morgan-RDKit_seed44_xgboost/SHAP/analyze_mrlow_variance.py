#!/usr/bin/env python3
"""
Analyze the causes of high MRLOW standard deviation
=====================================================
Investigate why MRLOW varies so widely within the same hybridization class.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Fragments
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent))
from modules.utils import load_config, DataLoader


def detect_boron_hybridization(smiles: str) -> str:
    """Detect boron hybridization state (same logic as analyze_mrlow_hybridization.py)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 'invalid'

    boron_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'B']
    if len(boron_atoms) == 0:
        return 'no_boron'

    b_atom = boron_atoms[0]
    degree = b_atom.GetDegree()
    total_valence = b_atom.GetTotalValence()
    formal_charge = b_atom.GetFormalCharge()

    if formal_charge == -1 and total_valence == 4:
        return 'sp3'
    if degree == 4:
        return 'sp3'
    if degree == 3 and formal_charge == -1 and total_valence == 4:
        return 'sp3'
    if degree == 2 and formal_charge == -1 and total_valence == 4:
        return 'sp3'
    if degree == 1 and formal_charge == -1 and total_valence == 4:
        return 'sp3'
    if degree == 3 and formal_charge == 0:
        return 'sp2'

    for neighbor in b_atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(b_atom.GetIdx(), neighbor.GetIdx())
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            return 'sp2'

    if degree == 2 and formal_charge == 0:
        return 'sp2'

    hybridization = b_atom.GetHybridization()
    if hybridization == Chem.HybridizationType.SP2:
        return 'sp2'
    elif hybridization == Chem.HybridizationType.SP3:
        return 'sp3'

    return 'unknown'


def analyze_boron_environment(smiles: str) -> dict:
    """Analyze the chemical environment around the boron atom."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    boron_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'B']
    if len(boron_atoms) == 0:
        return None

    b_atom = boron_atoms[0]

    # Count neighbor atom types
    neighbors = []
    for neighbor in b_atom.GetNeighbors():
        neighbors.append(neighbor.GetSymbol())

    # Tally each neighbor type
    neighbor_counts = {
        'C': neighbors.count('C'),
        'N': neighbors.count('N'),
        'O': neighbors.count('O'),
        'F': neighbors.count('F'),
        'Cl': neighbors.count('Cl'),
        'Br': neighbors.count('Br'),
        'H': neighbors.count('H'),
        'other': len([x for x in neighbors if x not in ['C', 'N', 'O', 'F', 'Cl', 'Br', 'H']])
    }

    # Detect special groups
    has_aromatic_ring = any([neighbor.GetIsAromatic() for neighbor in b_atom.GetNeighbors()])

    # Molecular-level properties
    try:
        mrlow = Descriptors.BCUT2D_MRLOW(mol)
        mol_wt = Descriptors.MolWt(mol)
        num_aromatic_rings = Descriptors.NumAromaticRings(mol)
        num_heteroatoms = Descriptors.NumHeteroatoms(mol)
        tpsa = Descriptors.TPSA(mol)

        # Electronic effect indicators
        num_nitro = Fragments.fr_nitro(mol)  # nitro group (strong electron-withdrawing)
        num_ether = Fragments.fr_ether(mol)  # ether linkage (weak electron-donating)
        num_amine = Fragments.fr_NH2(mol) + Fragments.fr_NH1(mol) + Fragments.fr_NH0(mol)

        return {
            'smiles': smiles,
            'mrlow': mrlow,
            'hybridization': detect_boron_hybridization(smiles),
            'b_neighbors_C': neighbor_counts['C'],
            'b_neighbors_N': neighbor_counts['N'],
            'b_neighbors_O': neighbor_counts['O'],
            'b_neighbors_F': neighbor_counts['F'],
            'b_neighbors_halogen': neighbor_counts['F'] + neighbor_counts['Cl'] + neighbor_counts['Br'],
            'b_has_aromatic_neighbor': has_aromatic_ring,
            'mol_wt': mol_wt,
            'num_aromatic_rings': num_aromatic_rings,
            'num_heteroatoms': num_heteroatoms,
            'tpsa': tpsa,
            'num_nitro': num_nitro,
            'num_ether': num_ether,
            'num_amine': num_amine,
        }
    except:
        return None


def main():
    print("="*70)
    print(" Analyzing causes of high MRLOW standard deviation")
    print("="*70)

    # Load data
    print("\n[Step 1/4] Loading data...")
    config = load_config()
    data_loader = DataLoader(config['paths'])
    X_test, y_test, smiles_df_test = data_loader.load_all_data('test')

    if 'Smiles' in smiles_df_test.columns:
        compound_smiles = smiles_df_test['Smiles'].tolist()
    else:
        compound_smiles = smiles_df_test.iloc[:, 0].tolist()

    print(f"✓ Loaded {len(compound_smiles)} compounds")

    # Analyze each molecule
    print("\n[Step 2/4] Analyzing chemical environments...")
    results = []
    for smiles in compound_smiles:
        env = analyze_boron_environment(smiles)
        if env:
            results.append(env)

    df = pd.DataFrame(results)

    # Separate sp2 and sp3
    df_sp2 = df[df['hybridization'] == 'sp2'].copy()
    df_sp3 = df[df['hybridization'] == 'sp3'].copy()

    print(f"✓ Analysis complete: sp2={len(df_sp2)}, sp3={len(df_sp3)}")

    # Correlate MRLOW with chemical environment factors
    print("\n" + "="*70)
    print("[Step 3/4] Correlation of MRLOW with chemical environment")
    print("="*70)

    # Correlations within sp2 group
    print("\n[SP2 group internal correlations]")
    print("-"*70)

    sp2_correlations = {}
    for col in ['b_neighbors_C', 'b_neighbors_N', 'b_neighbors_O',
                'b_neighbors_halogen', 'mol_wt', 'num_aromatic_rings',
                'num_heteroatoms', 'tpsa']:
        if col in df_sp2.columns:
            corr = df_sp2['mrlow'].corr(df_sp2[col])
            sp2_correlations[col] = corr

    sorted_sp2_corr = sorted(sp2_correlations.items(), key=lambda x: abs(x[1]), reverse=True)
    for feature, corr in sorted_sp2_corr:
        bar = '█' * int(abs(corr) * 50)
        sign = '+' if corr > 0 else '-'
        print(f"{feature:25s} {sign} {abs(corr):.3f}  {bar}")

    # Correlations within sp3 group
    print("\n[SP3 group internal correlations]")
    print("-"*70)

    sp3_correlations = {}
    for col in ['b_neighbors_C', 'b_neighbors_N', 'b_neighbors_O',
                'b_neighbors_halogen', 'mol_wt', 'num_aromatic_rings',
                'num_heteroatoms', 'tpsa']:
        if col in df_sp3.columns:
            corr = df_sp3['mrlow'].corr(df_sp3[col])
            sp3_correlations[col] = corr

    sorted_sp3_corr = sorted(sp3_correlations.items(), key=lambda x: abs(x[1]), reverse=True)
    for feature, corr in sorted_sp3_corr:
        bar = '█' * int(abs(corr) * 50)
        sign = '+' if corr > 0 else '-'
        print(f"{feature:25s} {sign} {abs(corr):.3f}  {bar}")

    # Group by neighbor atom type
    print("\n" + "="*70)
    print("[Step 4/4] Grouping by boron neighbor atom type")
    print("="*70)

    print("\n[SP2: grouped by number of C neighbors]")
    print("-"*70)
    sp2_by_c = df_sp2.groupby('b_neighbors_C')['mrlow'].agg(['count', 'mean', 'std'])
    print(sp2_by_c)

    print("\n[SP2: grouped by number of N neighbors]")
    print("-"*70)
    sp2_by_n = df_sp2.groupby('b_neighbors_N')['mrlow'].agg(['count', 'mean', 'std'])
    print(sp2_by_n)

    print("\n[SP2: grouped by number of O neighbors]")
    print("-"*70)
    sp2_by_o = df_sp2.groupby('b_neighbors_O')['mrlow'].agg(['count', 'mean', 'std'])
    print(sp2_by_o)

    print("\n[SP2: grouped by number of halogen neighbors]")
    print("-"*70)
    sp2_by_hal = df_sp2.groupby('b_neighbors_halogen')['mrlow'].agg(['count', 'mean', 'std'])
    print(sp2_by_hal)

    print("\n[SP2: aromatic vs non-aromatic neighbor]")
    print("-"*70)
    sp2_by_arom = df_sp2.groupby('b_has_aromatic_neighbor')['mrlow'].agg(['count', 'mean', 'std'])
    print(sp2_by_arom)

    # Visualization
    output_dir = Path(__file__).parent / 'outputs' / 'mrlow_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create distribution plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. MRLOW distribution histogram
    ax1 = axes[0, 0]
    ax1.hist(df_sp2['mrlow'], bins=30, alpha=0.6, label='sp2', color='blue', edgecolor='black')
    ax1.hist(df_sp3['mrlow'], bins=30, alpha=0.6, label='sp3', color='red', edgecolor='black')
    ax1.axvline(df_sp2['mrlow'].mean(), color='blue', linestyle='--', linewidth=2, label=f'sp2 mean={df_sp2["mrlow"].mean():.3f}')
    ax1.axvline(df_sp3['mrlow'].mean(), color='red', linestyle='--', linewidth=2, label=f'sp3 mean={df_sp3["mrlow"].mean():.3f}')
    ax1.set_xlabel('MRLOW', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('MRLOW Distribution by Hybridization', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(alpha=0.3)

    # 2. sp2 group: MRLOW vs number of C neighbors
    ax2 = axes[0, 1]
    sp2_c_groups = df_sp2.groupby('b_neighbors_C')['mrlow']
    positions = []
    data_to_plot = []
    labels = []
    for c_count, group in sp2_c_groups:
        if len(group) >= 3:  # need at least 3 samples for a box plot
            positions.append(c_count)
            data_to_plot.append(group.values)
            labels.append(f'{c_count}C\n(n={len(group)})')

    if data_to_plot:
        bp = ax2.boxplot(data_to_plot, positions=positions, labels=labels, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
        ax2.set_xlabel('Number of C neighbors', fontsize=12)
        ax2.set_ylabel('MRLOW', fontsize=12)
        ax2.set_title('SP2: MRLOW vs C Neighbors', fontsize=14, fontweight='bold')
        ax2.grid(alpha=0.3)

    # 3. sp2 group: MRLOW vs number of N neighbors
    ax3 = axes[1, 0]
    sp2_n_groups = df_sp2.groupby('b_neighbors_N')['mrlow']
    positions = []
    data_to_plot = []
    labels = []
    for n_count, group in sp2_n_groups:
        if len(group) >= 3:
            positions.append(n_count)
            data_to_plot.append(group.values)
            labels.append(f'{n_count}N\n(n={len(group)})')

    if data_to_plot:
        bp = ax3.boxplot(data_to_plot, positions=positions, labels=labels, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightgreen')
        ax3.set_xlabel('Number of N neighbors', fontsize=12)
        ax3.set_ylabel('MRLOW', fontsize=12)
        ax3.set_title('SP2: MRLOW vs N Neighbors', fontsize=14, fontweight='bold')
        ax3.grid(alpha=0.3)

    # 4. sp2 group: MRLOW vs aromatic neighbor
    ax4 = axes[1, 1]
    sp2_arom_groups = df_sp2.groupby('b_has_aromatic_neighbor')['mrlow']
    data_to_plot = []
    labels = []
    for has_arom, group in sp2_arom_groups:
        data_to_plot.append(group.values)
        label = 'Aromatic' if has_arom else 'Non-aromatic'
        labels.append(f'{label}\n(n={len(group)})')

    if len(data_to_plot) == 2:
        bp = ax4.boxplot(data_to_plot, labels=labels, patch_artist=True)
        for i, patch in enumerate(bp['boxes']):
            patch.set_facecolor(['lightyellow', 'lightcoral'][i])
        ax4.set_ylabel('MRLOW', fontsize=12)
        ax4.set_title('SP2: MRLOW vs Aromatic Neighbor', fontsize=14, fontweight='bold')
        ax4.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'mrlow_variance_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Plot saved: {output_dir / 'mrlow_variance_analysis.png'}")

    # Save detailed data
    csv_path = output_dir / 'mrlow_variance_detailed.csv'
    df.to_csv(csv_path, index=False)
    print(f"✓ Detailed data saved: {csv_path}")

    # Summary
    print("\n" + "="*70)
    print("Key Findings:")
    print("="*70)

    # Identify the strongest correlating factor
    top_sp2_factor = sorted_sp2_corr[0]
    top_sp3_factor = sorted_sp3_corr[0]

    print(f"\n1. Main factor driving MRLOW variance within the SP2 group:")
    print(f"   - {top_sp2_factor[0]}: correlation = {top_sp2_factor[1]:.3f}")

    print(f"\n2. Main factor driving MRLOW variance within the SP3 group:")
    print(f"   - {top_sp3_factor[0]}: correlation = {top_sp3_factor[1]:.3f}")

    print(f"\n3. Standard deviation analysis:")
    print(f"   - SP2 std dev: {df_sp2['mrlow'].std():.3f}")
    print(f"   - SP3 std dev: {df_sp3['mrlow'].std():.3f}")
    print(f"   - Even within the same hybridization class, variation in the boron chemical environment significantly affects MRLOW")

    print("\n" + "="*70)
    print("✓ Analysis complete!")
    print("="*70)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
