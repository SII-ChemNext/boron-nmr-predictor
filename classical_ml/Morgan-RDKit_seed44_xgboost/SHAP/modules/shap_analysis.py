"""
B11 NMR Interpretability Analysis - SHAP Global Analysis Module
===============================================================
Identify the most influential features for chemical shift prediction using SHAP values.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shap
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from .utils import FeatureMapper, DataLoader, ModelLoader, ensure_output_dir


class SHAPGlobalAnalyzer:
    """Global SHAP analyzer — identifies the most important features."""

    def __init__(
        self,
        model,
        feature_mapper: FeatureMapper,
        config: Dict,
        output_dir: str
    ):
        """
        Initialize the SHAP analyzer.

        Args:
            model: XGBoost model object
            feature_mapper: FeatureMapper instance
            config: SHAP analysis configuration dict
            output_dir: Output directory
        """
        self.model = model
        self.feature_mapper = feature_mapper
        self.config = config
        self.output_dir = Path(output_dir)
        ensure_output_dir(self.output_dir)

        # SHAP explainer (initialized during fit)
        self.explainer = None
        self.shap_values = None
        self.analysis_data = None

    def fit(self, X: np.ndarray, y: np.ndarray):
        """
        Compute SHAP values.

        Args:
            X: Feature matrix
            y: Labels (used for stratified sampling)
        """
        print("="*60)
        print("Global SHAP Analysis - Starting")
        print("="*60)

        # Step 1: Smart sampling
        X_background, X_analysis, y_analysis = self._smart_sampling(X, y)

        print(f"\n[1/3] Sampling complete:")
        print(f"  - Background dataset: {X_background.shape[0]} samples")
        print(f"  - Analysis dataset: {X_analysis.shape[0]} samples")
        print(f"  - Chemical shift range: [{y_analysis.min():.2f}, {y_analysis.max():.2f}] ppm")

        # Step 2: Initialize TreeExplainer
        print(f"\n[2/3] Initializing TreeExplainer...")
        # Use interventional method for stability without a large background set
        self.explainer = shap.TreeExplainer(
            self.model,
            data=X_background,
            feature_perturbation='interventional'
        )

        # Step 3: Compute SHAP values
        print(f"\n[3/3] Computing SHAP values...")
        self.shap_values = self.explainer.shap_values(X_analysis)
        self.analysis_data = X_analysis
        self.analysis_labels = y_analysis

        print(f"\n✓ SHAP values computed!")
        print(f"  - SHAP matrix shape: {self.shap_values.shape}")
        print(f"  - Expected value (baseline): {self.explainer.expected_value:.2f} ppm")

    def _smart_sampling(
        self, X: np.ndarray, y: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Stratified sampling strategy.

        Returns:
            (X_background, X_analysis, y_analysis)
        """
        n_samples = len(X)
        background_size = self.config['background_size']
        analysis_size = self.config['analysis_sample_size']

        # Use all data if dataset is small
        if n_samples <= analysis_size:
            print(f"  ⚠️  Dataset size ({n_samples}) <= analysis_sample_size ({analysis_size}), using all data")
            bg_indices = np.random.choice(n_samples, size=background_size, replace=False)
            return X[bg_indices], X, y

        # Stratified sampling by chemical shift tertiles
        bins = np.percentile(y, [0, 33, 67, 100])
        bin_labels = np.digitize(y, bins[1:-1])

        analysis_indices = []
        bg_indices = []

        for bin_id in range(3):
            bin_mask = (bin_labels == bin_id)
            bin_indices = np.where(bin_mask)[0]

            if len(bin_indices) == 0:
                continue

            # Analysis set sampling (proportional)
            n_from_bin = int(analysis_size * len(bin_indices) / n_samples)
            n_from_bin = min(n_from_bin, len(bin_indices))

            sampled = np.random.choice(
                bin_indices, size=n_from_bin, replace=False
            )
            analysis_indices.extend(sampled)

            # Background set sampling
            n_bg_from_bin = int(background_size * len(bin_indices) / n_samples)
            n_bg_from_bin = min(n_bg_from_bin, len(bin_indices))

            bg_sampled = np.random.choice(
                bin_indices, size=n_bg_from_bin, replace=False
            )
            bg_indices.extend(bg_sampled)

        analysis_indices = np.array(analysis_indices)
        bg_indices = np.array(bg_indices)

        return X[bg_indices], X[analysis_indices], y[analysis_indices]

    def analyze_feature_importance(self, top_k: int = 50) -> pd.DataFrame:
        """
        Analyze feature importance.

        Args:
            top_k: Number of top features to return

        Returns:
            Feature importance DataFrame
        """
        if self.shap_values is None:
            raise ValueError("Call fit() first to compute SHAP values.")

        print("\n" + "="*60)
        print("Feature Importance Analysis")
        print("="*60)

        # Compute mean absolute SHAP value per feature
        mean_abs_shap = np.abs(self.shap_values).mean(axis=0)

        importance_df = pd.DataFrame({
            'feature_idx': range(len(mean_abs_shap)),
            'mean_abs_shap': mean_abs_shap,
            'feature_name': self.feature_mapper.feature_names,
            'feature_block': [
                self.feature_mapper.get_block_name(i)
                for i in range(len(mean_abs_shap))
            ]
        })

        importance_df = importance_df.sort_values(
            'mean_abs_shap', ascending=False
        ).reset_index(drop=True)

        top_features = importance_df.head(top_k)

        print(f"\nTop-{top_k} most important features:")
        print("-" * 80)
        for i, row in top_features.head(10).iterrows():
            print(f"{i+1:2d}. [{row['feature_block']:20s}] "
                  f"{row['feature_name']:40s} | "
                  f"SHAP: {row['mean_abs_shap']:.4f}")

        save_path = self.output_dir / 'feature_importance_ranking.csv'
        importance_df.to_csv(save_path, index=False)
        print(f"\n✓ Full feature importance ranking saved to: {save_path}")

        return top_features

    def analyze_feature_blocks(self) -> Dict[str, float]:
        """
        Aggregate SHAP values by feature block.

        Returns:
            Mean absolute SHAP value per block
        """
        if self.shap_values is None:
            raise ValueError("Call fit() first to compute SHAP values.")

        print("\n" + "="*60)
        print("Feature Block Contribution Analysis")
        print("="*60)

        block_summary = self.feature_mapper.summarize_shap_by_block(self.shap_values)

        total_shap = sum(block_summary.values())
        block_percentages = {
            k: (v / total_shap * 100) for k, v in block_summary.items()
        }

        print("\nMean contribution per feature block:")
        print("-" * 60)
        for block_name, percentage in sorted(
            block_percentages.items(), key=lambda x: x[1], reverse=True
        ):
            print(f"  {block_name:25s}: {percentage:5.2f}%  "
                  f"(mean |SHAP|: {block_summary[block_name]:.4f})")

        compound_contribution = (
            block_summary['compound_morgan'] + block_summary['compound_rdkit']
        )
        solvent_contribution = (
            block_summary['solvent_morgan'] + block_summary['solvent_rdkit']
        )

        print("\nHigh-level summary:")
        print("-" * 60)
        print(f"  Compound features total: {compound_contribution/total_shap*100:.2f}%")
        print(f"  Solvent features total:  {solvent_contribution/total_shap*100:.2f}%")

        morgan_contribution = (
            block_summary['compound_morgan'] + block_summary['solvent_morgan']
        )
        rdkit_contribution = (
            block_summary['compound_rdkit'] + block_summary['solvent_rdkit']
        )

        print(f"\n  Morgan fingerprints total: {morgan_contribution/total_shap*100:.2f}%")
        print(f"  RDKit descriptors total:   {rdkit_contribution/total_shap*100:.2f}%")

        summary_data = {
            'block_summary': block_summary,
            'block_percentages': block_percentages,
            'high_level_summary': {
                'compound_vs_solvent': {
                    'compound': compound_contribution / total_shap * 100,
                    'solvent': solvent_contribution / total_shap * 100
                },
                'morgan_vs_rdkit': {
                    'morgan': morgan_contribution / total_shap * 100,
                    'rdkit': rdkit_contribution / total_shap * 100
                }
            }
        }

        import json
        with open(self.output_dir / 'feature_block_summary.json', 'w') as f:
            json.dump(summary_data, f, indent=2)

        return block_summary

    def plot_summary(self, max_display: int = 20):
        """
        Generate SHAP summary (beeswarm) plot.

        Args:
            max_display: Maximum number of features to display
        """
        if self.shap_values is None:
            raise ValueError("Call fit() first to compute SHAP values.")

        print(f"\nGenerating SHAP Summary Plot (top-{max_display} features)...")

        plt.figure(figsize=(12, 8))

        shap.summary_plot(
            self.shap_values,
            self.analysis_data,
            feature_names=self.feature_mapper.feature_names,
            max_display=max_display,
            show=False
        )

        plt.title('SHAP Summary Plot - B11 NMR Chemical Shift Prediction',
                  fontsize=14, pad=20)
        plt.xlabel('SHAP value (impact on prediction)', fontsize=12)
        plt.tight_layout()

        save_path = self.output_dir / 'shap_summary_plot.png'
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"✓ Summary plot saved to: {save_path}")

    def plot_bar(self, max_display: int = 20):
        """
        Generate SHAP bar plot.

        Args:
            max_display: Maximum number of features to display
        """
        if self.shap_values is None:
            raise ValueError("Call fit() first to compute SHAP values.")

        print(f"\nGenerating SHAP Bar Plot (top-{max_display} features)...")

        plt.figure(figsize=(10, 8))

        shap.summary_plot(
            self.shap_values,
            self.analysis_data,
            feature_names=self.feature_mapper.feature_names,
            max_display=max_display,
            plot_type='bar',
            show=False
        )

        plt.title('Feature Importance (mean |SHAP value|)', fontsize=14, pad=20)
        plt.xlabel('mean(|SHAP value|)', fontsize=12)
        plt.tight_layout()

        save_path = self.output_dir / 'shap_bar_plot.png'
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"✓ Bar plot saved to: {save_path}")

    def plot_block_contribution_pie(self):
        """Generate pie charts showing feature block contributions."""
        if self.shap_values is None:
            raise ValueError("Call fit() first to compute SHAP values.")

        print("\nGenerating feature block contribution pie charts...")

        block_summary = self.feature_mapper.summarize_shap_by_block(self.shap_values)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Pie chart 1: four feature blocks
        colors1 = ['#FF6B6B', '#FFA07A', '#4ECDC4', '#95E1D3']
        labels1 = [
            f"Compound\nMorgan\n({block_summary['compound_morgan']:.3f})",
            f"Compound\nRDKit\n({block_summary['compound_rdkit']:.3f})",
            f"Solvent\nMorgan\n({block_summary['solvent_morgan']:.3f})",
            f"Solvent\nRDKit\n({block_summary['solvent_rdkit']:.3f})"
        ]
        sizes1 = [
            block_summary['compound_morgan'],
            block_summary['compound_rdkit'],
            block_summary['solvent_morgan'],
            block_summary['solvent_rdkit']
        ]

        ax1.pie(sizes1, labels=labels1, colors=colors1, autopct='%1.1f%%',
                startangle=90, textprops={'fontsize': 10})
        ax1.set_title('Contribution by Feature Block', fontsize=14, pad=20)

        # Pie chart 2: compound vs solvent
        compound_total = block_summary['compound_morgan'] + block_summary['compound_rdkit']
        solvent_total = block_summary['solvent_morgan'] + block_summary['solvent_rdkit']

        colors2 = ['#FF6B6B', '#4ECDC4']
        labels2 = [
            f"Compound\nFeatures\n({compound_total:.3f})",
            f"Solvent\nFeatures\n({solvent_total:.3f})"
        ]
        sizes2 = [compound_total, solvent_total]

        ax2.pie(sizes2, labels=labels2, colors=colors2, autopct='%1.1f%%',
                startangle=90, textprops={'fontsize': 11})
        ax2.set_title('Compound vs Solvent Contribution', fontsize=14, pad=20)

        plt.tight_layout()

        save_path = self.output_dir / 'feature_block_pie_charts.png'
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"✓ Pie charts saved to: {save_path}")

    def run_full_analysis(self, X: np.ndarray, y: np.ndarray, top_k: int = 50):
        """
        Run the complete global SHAP analysis pipeline.

        Args:
            X: Feature matrix
            y: Labels
            top_k: Number of top features

        Returns:
            Feature importance DataFrame
        """
        # 1. Compute SHAP values
        self.fit(X, y)

        # 2. Feature importance analysis
        top_features = self.analyze_feature_importance(top_k=top_k)

        # 3. Feature block analysis
        block_summary = self.analyze_feature_blocks()

        # 4. Visualization
        self.plot_summary(max_display=20)
        self.plot_bar(max_display=20)
        self.plot_block_contribution_pie()

        print("\n" + "="*60)
        print("✓ Global SHAP analysis complete!")
        print(f"✓ All results saved to: {self.output_dir}")
        print("="*60)

        return top_features


def main():
    """Test and demonstration code."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent.parent))

    from modules.utils import load_config, FeatureMapper, DataLoader, ModelLoader

    config = load_config()

    feature_mapper = FeatureMapper(config['features'])
    data_loader = DataLoader(config['paths'])

    print("Loading data...")
    X_train, y_train, _ = data_loader.load_all_data('train')

    print("Loading model...")
    model_path = Path(__file__).parent.parent / config['paths']['model_path']
    model = ModelLoader.load_xgboost_model(str(model_path.resolve()))

    output_dir = Path(__file__).parent.parent / config['paths']['output_dir'] / 'shap_global'
    analyzer = SHAPGlobalAnalyzer(
        model=model,
        feature_mapper=feature_mapper,
        config=config['shap_analysis'],
        output_dir=str(output_dir)
    )

    top_features = analyzer.run_full_analysis(X_train, y_train, top_k=50)

    print("\nTop-10 most important features:")
    print(top_features.head(10))


if __name__ == "__main__":
    main()
