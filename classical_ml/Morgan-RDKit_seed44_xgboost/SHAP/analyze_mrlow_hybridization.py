#!/usr/bin/env python3
"""
Analyze the relationship between BCUT2D_MRLOW and boron hybridization state
============================================================================
Investigate the correlation between MRLOW feature values and
sp2/sp3 hybridization of the boron atom.

Hypothesis:
- sp2 boron (tricoordinate): typically in boronic acids and esters; planar geometry
- sp3 boron (tetracoordinate): typically in borate anions; tetrahedral geometry
- MRLOW may correlate with hybridization as it encodes electron distribution
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent))

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

from modules.utils import load_config, DataLoader


def detect_boron_hybridization(smiles: str) -> str:
    """
    Detect the hybridization state of the boron atom (corrected version).

    Args:
        smiles: SMILES string

    Returns:
        'sp2', 'sp3', 'unknown', or 'no_boron'

    Logic:
    - Detect BH3/BH2/BH complexes (implicit H causes inaccurate RDKit coordination count)
    - Tetracoordinate or negatively charged -> sp3
    - Tricoordinate and neutral -> sp2
    """
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
    num_implicit_hs = b_atom.GetNumImplicitHs()
    total_num_hs = b_atom.GetTotalNumHs()

    # Rule 1: BH3/BH2/BH complexes with negative charge and valence 4 -> sp3
    if formal_charge == -1 and total_valence == 4:
        return 'sp3'

    # Rule 2: tetracoordinate boron -> sp3
    if degree == 4:
        return 'sp3'

    # Rule 3: tricoordinate, negatively charged, total_valence=4 -> sp3 (has implicit H)
    if degree == 3 and formal_charge == -1 and total_valence == 4:
        return 'sp3'

    # Rule 4: dicoordinate, negatively charged, total_valence=4 -> sp3 (BH2 complex)
    if degree == 2 and formal_charge == -1 and total_valence == 4:
        return 'sp3'

    # Rule 5: monocoordinate, negatively charged, total_valence=4 -> sp3 (BH3 complex)
    if degree == 1 and formal_charge == -1 and total_valence == 4:
        return 'sp3'

    # Rule 6: tricoordinate and neutral -> sp2
    if degree == 3 and formal_charge == 0:
        return 'sp2'

    # Rule 7: B=N double bond present -> sp2
    for neighbor in b_atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(b_atom.GetIdx(), neighbor.GetIdx())
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            return 'sp2'

    # Rule 8: dicoordinate and neutral -> likely BN heterocycle (sp2)
    if degree == 2 and formal_charge == 0:
        return 'sp2'

    # Fallback: use RDKit hybridization
    hybridization = b_atom.GetHybridization()
    if hybridization == Chem.HybridizationType.SP2:
        return 'sp2'
    elif hybridization == Chem.HybridizationType.SP3:
        return 'sp3'

    return 'unknown'


def get_boron_coordination_info(smiles: str) -> dict:
    """
    Get detailed coordination information for the boron atom.

    Returns:
        {
            'hybridization': str,
            'degree': int,
            'formal_charge': int,
            'total_valence': int,
            'neighbors': list of str
        }
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    boron_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'B']
    if len(boron_atoms) == 0:
        return None

    b_atom = boron_atoms[0]

    neighbors = [neighbor.GetSymbol() for neighbor in b_atom.GetNeighbors()]

    return {
        'hybridization': str(b_atom.GetHybridization()),
        'degree': b_atom.GetDegree(),
        'formal_charge': b_atom.GetFormalCharge(),
        'total_valence': b_atom.GetTotalValence(),
        'neighbors': neighbors,
        'is_aromatic': b_atom.GetIsAromatic()
    }


def calculate_mrlow(smiles: str) -> float:
    """Calculate the BCUT2D_MRLOW descriptor."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.nan

    try:
        mrlow = Descriptors.BCUT2D_MRLOW(mol)
        return mrlow
    except:
        return np.nan


def main():
    print("="*70)
    print(" BCUT2D_MRLOW vs Boron Hybridization State Analysis")
    print("="*70)

    # Step 1: Load configuration and data
    print("\n[Step 1/4] Loading data...")
    config = load_config()
    data_loader = DataLoader(config['paths'])

    X_test, y_test, smiles_df_test = data_loader.load_all_data('test')

    if 'Smiles' in smiles_df_test.columns:
        compound_smiles = smiles_df_test['Smiles'].tolist()
    elif 'smiles' in smiles_df_test.columns:
        compound_smiles = smiles_df_test['smiles'].tolist()
    else:
        compound_smiles = smiles_df_test.iloc[:, 0].tolist()

    print(f"✓ Loaded {len(compound_smiles)} compounds")

    # Step 2: Analyze each compound
    print("\n[Step 2/4] Analyzing boron hybridization and MRLOW...")

    results = []

    for idx, smiles in enumerate(compound_smiles):
        if (idx + 1) % 50 == 0:
            print(f"  Progress: {idx+1}/{len(compound_smiles)}")

        hybridization = detect_boron_hybridization(smiles)
        coord_info = get_boron_coordination_info(smiles)
        mrlow = calculate_mrlow(smiles)
        ppm = y_test[idx]

        results.append({
            'smiles': smiles,
            'hybridization': hybridization,
            'mrlow': mrlow,
            'ppm': ppm,
            'degree': coord_info['degree'] if coord_info else None,
            'formal_charge': coord_info['formal_charge'] if coord_info else None,
            'neighbors': str(coord_info['neighbors']) if coord_info else None
        })

    df = pd.DataFrame(results)

    print(f"\n✓ Analysis complete")

    # Step 3: Statistical analysis
    print("\n[Step 3/4] Statistical analysis...")
    print("\n" + "="*70)
    print("Hybridization distribution:")
    print("="*70)

    hyb_counts = df['hybridization'].value_counts()
    print(hyb_counts)

    print("\n" + "="*70)
    print("MRLOW statistics by hybridization state:")
    print("="*70)

    valid_df = df[df['mrlow'].notna() & (df['hybridization'] != 'no_boron')].copy()

    valid_df['hybridization_simple'] = valid_df['hybridization'].apply(
        lambda x: 'sp2' if 'sp2' in x else ('sp3' if 'sp3' in x else 'other')
    )

    grouped = valid_df.groupby('hybridization_simple')['mrlow'].agg([
        'count', 'mean', 'std', 'min', 'max',
        ('Q25', lambda x: x.quantile(0.25)),
        ('median', 'median'),
        ('Q75', lambda x: x.quantile(0.75))
    ])

    print(grouped)

    print("\n" + "="*70)
    print("Statistical tests (sp2 vs sp3):")
    print("="*70)

    sp2_mrlow = valid_df[valid_df['hybridization_simple'] == 'sp2']['mrlow']
    sp3_mrlow = valid_df[valid_df['hybridization_simple'] == 'sp3']['mrlow']

    if len(sp2_mrlow) > 0 and len(sp3_mrlow) > 0:
        from scipy import stats

        t_stat, t_pval = stats.ttest_ind(sp2_mrlow, sp3_mrlow)
        print(f"\nt-test:")
        print(f"  t-statistic: {t_stat:.4f}")
        print(f"  p-value: {t_pval:.4e}")

        u_stat, u_pval = stats.mannwhitneyu(sp2_mrlow, sp3_mrlow)
        print(f"\nMann-Whitney U test:")
        print(f"  U-statistic: {u_stat:.4f}")
        print(f"  p-value: {u_pval:.4e}")

        cohen_d = (sp2_mrlow.mean() - sp3_mrlow.mean()) / np.sqrt(
            (sp2_mrlow.std()**2 + sp3_mrlow.std()**2) / 2
        )
        print(f"\nCohen's d (effect size): {cohen_d:.4f}")

        if abs(cohen_d) < 0.2:
            effect_interpretation = "negligible"
        elif abs(cohen_d) < 0.5:
            effect_interpretation = "small"
        elif abs(cohen_d) < 0.8:
            effect_interpretation = "medium"
        else:
            effect_interpretation = "large"

        print(f"  Interpretation: {effect_interpretation} effect")

    print("\n" + "="*70)
    print("Chemical shift statistics by hybridization state:")
    print("="*70)

    ppm_grouped = valid_df.groupby('hybridization_simple')['ppm'].agg([
        'count', 'mean', 'std', 'min', 'max', 'median'
    ])
    print(ppm_grouped)

    # Step 4: Visualization
    print("\n[Step 4/4] Generating visualizations...")

    output_dir = Path(__file__).parent / 'outputs' / 'mrlow_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1.1 Box plot
    ax = axes[0, 0]
    valid_df.boxplot(column='mrlow', by='hybridization_simple', ax=ax)
    ax.set_title('MRLOW Distribution by Hybridization State')
    ax.set_xlabel('Boron hybridization')
    ax.set_ylabel('BCUT2D_MRLOW')
    plt.sca(ax)
    plt.xticks(rotation=0)

    # 1.2 Violin plot
    ax = axes[0, 1]
    if len(valid_df['hybridization_simple'].unique()) > 1:
        sns.violinplot(data=valid_df, x='hybridization_simple', y='mrlow', ax=ax)
        ax.set_title('MRLOW Density Distribution by Hybridization State')
        ax.set_xlabel('Boron hybridization')
        ax.set_ylabel('BCUT2D_MRLOW')

    # 1.3 MRLOW vs chemical shift (sp2)
    ax = axes[1, 0]
    sp2_data = valid_df[valid_df['hybridization_simple'] == 'sp2']
    if len(sp2_data) > 0:
        ax.scatter(sp2_data['mrlow'], sp2_data['ppm'], alpha=0.5, s=30)
        ax.set_title('sp2 Boron: MRLOW vs Chemical Shift')
        ax.set_xlabel('BCUT2D_MRLOW')
        ax.set_ylabel('δ(B11) / ppm')
        ax.grid(True, alpha=0.3)

        if len(sp2_data) > 2:
            z = np.polyfit(sp2_data['mrlow'], sp2_data['ppm'], 1)
            p = np.poly1d(z)
            ax.plot(sp2_data['mrlow'].sort_values(),
                   p(sp2_data['mrlow'].sort_values()),
                   "r--", alpha=0.8, linewidth=2,
                   label=f'y={z[0]:.2f}x+{z[1]:.2f}')
            ax.legend()

    # 1.4 MRLOW vs chemical shift (sp3)
    ax = axes[1, 1]
    sp3_data = valid_df[valid_df['hybridization_simple'] == 'sp3']
    if len(sp3_data) > 0:
        ax.scatter(sp3_data['mrlow'], sp3_data['ppm'], alpha=0.5, s=30, color='orange')
        ax.set_title('sp3 Boron: MRLOW vs Chemical Shift')
        ax.set_xlabel('BCUT2D_MRLOW')
        ax.set_ylabel('δ(B11) / ppm')
        ax.grid(True, alpha=0.3)

        if len(sp3_data) > 2:
            z = np.polyfit(sp3_data['mrlow'], sp3_data['ppm'], 1)
            p = np.poly1d(z)
            ax.plot(sp3_data['mrlow'].sort_values(),
                   p(sp3_data['mrlow'].sort_values()),
                   "r--", alpha=0.8, linewidth=2,
                   label=f'y={z[0]:.2f}x+{z[1]:.2f}')
            ax.legend()

    plt.tight_layout()
    fig_path = output_dir / 'mrlow_hybridization_analysis.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Plot saved: {fig_path}")

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    colors = {'sp2': 'steelblue', 'sp3': 'coral', 'other': 'gray'}

    for hyb_type in valid_df['hybridization_simple'].unique():
        data = valid_df[valid_df['hybridization_simple'] == hyb_type]
        ax.scatter(data['mrlow'], data['ppm'],
                  alpha=0.6, s=50,
                  color=colors.get(hyb_type, 'gray'),
                  label=f'{hyb_type} (n={len(data)})')

    ax.set_xlabel('BCUT2D_MRLOW', fontsize=12)
    ax.set_ylabel('δ(B11) / ppm', fontsize=12)
    ax.set_title('MRLOW vs B11 Chemical Shift by Boron Hybridization', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig_path2 = output_dir / 'mrlow_vs_ppm_by_hybridization.png'
    plt.savefig(fig_path2, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Plot saved: {fig_path2}")

    csv_path = output_dir / 'mrlow_hybridization_data.csv'
    valid_df.to_csv(csv_path, index=False)
    print(f"✓ Data saved: {csv_path}")

    summary_path = output_dir / 'statistical_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("BCUT2D_MRLOW vs Boron Hybridization - Statistical Summary\n")
        f.write("="*70 + "\n\n")

        f.write("Hybridization distribution:\n")
        f.write(str(hyb_counts) + "\n\n")

        f.write("MRLOW statistics by hybridization state:\n")
        f.write(str(grouped) + "\n\n")

        if len(sp2_mrlow) > 0 and len(sp3_mrlow) > 0:
            f.write("Statistical tests (sp2 vs sp3):\n")
            f.write(f"  t-test p-value: {t_pval:.4e}\n")
            f.write(f"  Mann-Whitney U p-value: {u_pval:.4e}\n")
            f.write(f"  Cohen's d: {cohen_d:.4f} ({effect_interpretation} effect)\n\n")

        f.write("Chemical shift statistics by hybridization state:\n")
        f.write(str(ppm_grouped) + "\n")

    print(f"✓ Statistical summary saved: {summary_path}")

    print("\n" + "="*70)
    print("✓ Analysis complete!")
    print(f"✓ All results saved to: {output_dir}")
    print("="*70)

    print("\n" + "="*70)
    print("Key Findings:")
    print("="*70)

    if len(sp2_mrlow) > 0 and len(sp3_mrlow) > 0:
        print(f"\n1. sp2 boron: mean MRLOW = {sp2_mrlow.mean():.4f} ± {sp2_mrlow.std():.4f}")
        print(f"2. sp3 boron: mean MRLOW = {sp3_mrlow.mean():.4f} ± {sp3_mrlow.std():.4f}")
        print(f"3. Significance: p = {t_pval:.4e} {'(significant)' if t_pval < 0.05 else '(not significant)'}")
        print(f"4. Effect size: {effect_interpretation}")

        if t_pval < 0.05:
            if sp2_mrlow.mean() > sp3_mrlow.mean():
                print(f"\n✓ sp2 boron has significantly higher MRLOW than sp3 boron")
                print("  → Tricoordinate boron shows higher molar refractivity eigenvalue")
                print("  → May reflect the planar geometry of sp2 hybridization")
            else:
                print(f"\n✓ sp3 boron has significantly higher MRLOW than sp2 boron")
                print("  → Tetracoordinate boron shows higher molar refractivity eigenvalue")
                print("  → May reflect the tetrahedral geometry of sp3 hybridization")
        else:
            print(f"\n⚠️ No significant difference in MRLOW between sp2 and sp3 boron")
            print("  → MRLOW may not be the primary feature distinguishing hybridization states")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
