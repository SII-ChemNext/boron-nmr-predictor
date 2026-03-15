#!/usr/bin/env python3
"""
Run Morgan fingerprint substructure analysis
=============================================
Identify and visualize molecular substructures most influential for B11 NMR chemical shift prediction.

Usage:
    python run_morgan_analysis.py

Outputs:
    - outputs/morgan_analysis/morgan_bit_statistics.csv        (bit statistics)
    - outputs/morgan_analysis/bit_XXX_substructures.png        (substructure visualizations)
    - outputs/morgan_analysis/morgan_substructure_summary.md   (analysis report)
"""

import sys
from pathlib import Path
import time
import pandas as pd

# Add module path
sys.path.insert(0, str(Path(__file__).parent))

from modules.utils import load_config, DataLoader
from modules.morgan_substructure import MorganSubstructureAnalyzer


def main():
    """Main execution function."""
    print("="*70)
    print(" B11 NMR Chemical Shift Prediction Model - Morgan Fingerprint Substructure Analysis")
    print("="*70)

    start_time = time.time()

    # Step 1: Load configuration
    print("\n[Step 1/4] Loading configuration...")
    config = load_config()
    print("✓ Configuration loaded")

    # Step 2: Load data (retrieve SMILES)
    print("\n[Step 2/4] Loading compound data...")
    data_loader = DataLoader(config['paths'])

    # Load test set data for analysis
    X_test, y_test, smiles_df_test = data_loader.load_all_data('test')

    print(f"✓ Data loaded")
    print(f"  - Test set size: {len(smiles_df_test)}")
    print(f"  - Chemical shift range: [{y_test.min():.2f}, {y_test.max():.2f}] ppm")

    # Extract compound SMILES
    # smiles_df should have 'Smiles' (compound) and 'solvent' columns
    if 'Smiles' in smiles_df_test.columns:
        compound_smiles = smiles_df_test['Smiles'].tolist()
    elif 'smiles' in smiles_df_test.columns:
        compound_smiles = smiles_df_test['smiles'].tolist()
    else:
        # Fallback: assume first column is SMILES
        compound_smiles = smiles_df_test.iloc[:, 0].tolist()

    print(f"  - Extracted compound SMILES: {len(compound_smiles)}")

    # Step 3: Load SHAP feature importance ranking
    print("\n[Step 3/4] Loading SHAP feature importance ranking...")
    importance_path = Path(__file__).parent / 'outputs/shap_global/feature_importance_ranking.csv'

    if not importance_path.exists():
        print("❌ Error: SHAP feature importance file not found!")
        print("   Please run first: python run_shap_analysis.py")
        sys.exit(1)

    importance_df = pd.read_csv(importance_path)

    # Filter Morgan fingerprint features
    morgan_features = importance_df[
        importance_df['feature_block'] == 'compound_morgan'
    ].head(20)  # Take Top-20 Morgan bits

    print(f"✓ Feature importance loaded")
    print(f"  - Top-20 Morgan bits:")

    for idx, row in morgan_features.head(10).iterrows():
        # Extract bit index from feature name "compound_morgan_bit_XXX"
        bit_idx = int(row['feature_name'].split('_')[-1])
        print(f"    {idx+1:2d}. Bit {bit_idx:4d} | SHAP: {row['mean_abs_shap']:.4f}")

    # Step 4: Morgan substructure analysis
    print("\n[Step 4/4] Starting Morgan fingerprint substructure analysis...")
    output_dir = Path(__file__).parent / config['paths']['output_dir'] / 'morgan_analysis'

    analyzer = MorganSubstructureAnalyzer(
        output_dir=str(output_dir),
        radius=2,  # Morgan-2
        n_bits=1024
    )

    # Extract bit indices and SHAP values
    important_bits = []
    bit_names = []
    shap_values = []

    for _, row in morgan_features.iterrows():
        bit_idx = int(row['feature_name'].split('_')[-1])
        important_bits.append(bit_idx)
        bit_names.append(row['feature_name'])
        shap_values.append(row['mean_abs_shap'])

    # Analyze these bits
    stats_df = analyzer.analyze_important_bits(
        smiles_list=compound_smiles,
        important_bits=important_bits,
        bit_names=bit_names,
        shap_values=shap_values
    )

    # Visualize Top-10 bits
    print("\n" + "="*70)
    print("Visualizing substructures for Top-10 Morgan bits")
    print("="*70)

    for idx in range(min(10, len(important_bits))):
        bit_idx = important_bits[idx]
        print(f"\n[{idx+1}/10] ", end="")
        analyzer.visualize_bit_substructures(
            bit_idx=bit_idx,
            n_examples=6,
            mols_per_row=3
        )

    # Generate summary report
    print("\n" + "="*70)
    print("Generating analysis report")
    print("="*70)
    analyzer.generate_summary_report(stats_df)

    # Summary
    elapsed_time = time.time() - start_time
    print("\n" + "="*70)
    print(f"✓ Analysis complete! Total time: {elapsed_time:.1f} s")
    print(f"✓ All results saved to: {output_dir}")
    print("="*70)

    # List generated files
    print("\nGenerated files:")
    for file_path in sorted(output_dir.glob('*')):
        print(f"  - {file_path.name}")

    print("\nSuggested next steps:")
    print("  1. See morgan_bit_statistics.csv for bit statistics")
    print("  2. See bit_XXX_substructures.png for substructure visualizations")
    print("  3. See morgan_substructure_summary.md for the full report")
    print("  4. Interpret these substructures in terms of their effect on B11 NMR shift")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
