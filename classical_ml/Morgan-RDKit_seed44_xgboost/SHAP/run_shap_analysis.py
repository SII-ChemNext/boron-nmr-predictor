#!/usr/bin/env python3
"""
Run global SHAP analysis
=========================
Entry-point script for SHAP interpretability analysis of the B11 NMR model.

Usage:
    python run_shap_analysis.py

Outputs:
    - outputs/shap_global/feature_importance_ranking.csv  (full feature ranking)
    - outputs/shap_global/feature_block_summary.json      (feature block summary)
    - outputs/shap_global/shap_summary_plot.png           (SHAP beeswarm plot)
    - outputs/shap_global/shap_bar_plot.png               (feature importance bar chart)
    - outputs/shap_global/feature_block_pie_charts.png    (feature block pie charts)
"""

import sys
from pathlib import Path
import time

# Add module path
sys.path.insert(0, str(Path(__file__).parent))

from modules.utils import load_config, FeatureMapper, DataLoader, ModelLoader
from modules.shap_analysis import SHAPGlobalAnalyzer


def main():
    """Main execution function."""
    print("="*70)
    print(" B11 NMR Chemical Shift Prediction Model - SHAP Interpretability Analysis")
    print("="*70)

    start_time = time.time()

    # Step 1: Load configuration
    print("\n[Step 1/5] Loading configuration...")
    config = load_config()
    print("✓ Configuration loaded")

    # Step 2: Initialize tools
    print("\n[Step 2/5] Initializing feature mapper and data loader...")
    feature_mapper = FeatureMapper(config['features'])
    data_loader = DataLoader(config['paths'])
    print(f"✓ Feature mapper initialized (total features: {feature_mapper.total_dims})")

    # Step 3: Load data
    print("\n[Step 3/5] Loading data...")
    print("  Strategy: training set as background + test set as analysis target")
    X_train, y_train, smiles_train = data_loader.load_all_data('train')
    X_test, y_test, smiles_test = data_loader.load_all_data('test')
    print(f"✓ Data loaded")
    print(f"  - Train set (background): {len(X_train)} samples")
    print(f"  - Test set (analysis): {len(X_test)} samples")
    print(f"  - Feature dimensions: {X_train.shape[1]}")
    print(f"  - Train shift range: [{y_train.min():.2f}, {y_train.max():.2f}] ppm")
    print(f"  - Test shift range: [{y_test.min():.2f}, {y_test.max():.2f}] ppm")

    # Step 4: Load model
    print("\n[Step 4/5] Loading XGBoost model...")
    model_path = Path(__file__).parent / config['paths']['model_path']
    model = ModelLoader.load_xgboost_model(str(model_path.resolve()))
    print(f"✓ Model loaded")
    print(f"  - Model file: {model_path.name}")

    # Step 5: SHAP analysis
    print("\n[Step 5/5] Running global SHAP analysis...")
    output_dir = Path(__file__).parent / config['paths']['output_dir'] / 'shap_global'

    analyzer = SHAPGlobalAnalyzer(
        model=model,
        feature_mapper=feature_mapper,
        config=config['shap_analysis'],
        output_dir=str(output_dir)
    )

    # Run full analysis pipeline
    # Strategy: training set as background (baseline), test set as analysis target (generalization)
    # First fit background with training set, then analyze test set
    print("  Building SHAP background baseline from training set...")
    analyzer.fit(X_train, y_train)  # Use training set to establish background

    print("  Analyzing test set samples...")
    # Manually set analysis data to test set
    analyzer.analysis_data = X_test
    analyzer.analysis_labels = y_test
    # Recompute SHAP values for test set
    analyzer.shap_values = analyzer.explainer.shap_values(X_test)

    # Then run analysis
    top_features = analyzer.analyze_feature_importance(
        top_k=config['shap_analysis']['global']['top_k_features']
    )
    analyzer.analyze_feature_blocks()
    analyzer.plot_summary(max_display=20)
    analyzer.plot_bar(max_display=20)
    analyzer.plot_block_contribution_pie()

    # Print Top-10 features
    print("\n" + "="*70)
    print(" Top-10 Most Important Features")
    print("="*70)
    print(f"\n{'Rank':<6}{'Feature Block':<25}{'Feature Name':<45}{'SHAP':<10}")
    print("-"*90)

    for i, row in top_features.head(10).iterrows():
        print(f"{i+1:<6}{row['feature_block']:<25}{row['feature_name']:<45}{row['mean_abs_shap']:<10.5f}")

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
    print("  1. See feature_importance_ranking.csv for the full feature ranking")
    print("  2. See shap_summary_plot.png for feature value vs SHAP value relationships")
    print("  3. See feature_block_pie_charts.png for compound vs solvent contribution")
    print("  4. Increase analysis_sample_size to 500 in config for higher precision")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
