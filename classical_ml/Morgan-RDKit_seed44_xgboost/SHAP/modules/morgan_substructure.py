"""
B11 NMR Interpretability Analysis - Morgan Fingerprint Substructure Module
=============================================================================
Identify and visualize molecular substructures most influential for chemical shift prediction.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, DataStructs
from rdkit.Chem.Draw import IPythonConsole
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from .utils import ensure_output_dir


class MorganSubstructureAnalyzer:
    """Morgan fingerprint substructure analyzer."""

    def __init__(
        self,
        output_dir: str,
        radius: int = 2,
        n_bits: int = 1024
    ):
        """
        Initialize the Morgan substructure analyzer.

        Args:
            output_dir: Output directory
            radius: Morgan fingerprint radius (default 2, i.e. Morgan-2)
            n_bits: Fingerprint length (default 1024)
        """
        self.output_dir = Path(output_dir)
        ensure_output_dir(self.output_dir)

        self.radius = radius
        self.n_bits = n_bits

        # Store analysis results
        self.bit_info_dict = {}  # {bit_idx: {mol_idx: [atom_ids]}}

    def analyze_important_bits(
        self,
        smiles_list: List[str],
        important_bits: List[int],
        bit_names: List[str] = None,
        shap_values: List[float] = None
    ) -> pd.DataFrame:
        """
        Analyze substructures corresponding to important Morgan bits.

        Args:
            smiles_list: List of SMILES strings
            important_bits: List of bit indices to analyze
            bit_names: Display names for each bit (optional)
            shap_values: SHAP importance value for each bit (optional)

        Returns:
            Analysis results DataFrame
        """
        print("="*70)
        print("Morgan Fingerprint Substructure Analysis")
        print("="*70)
        print(f"\nAnalysis configuration:")
        print(f"  - Morgan radius: {self.radius}")
        print(f"  - Fingerprint length: {self.n_bits}")
        print(f"  - Bits to analyze: {len(important_bits)}")
        print(f"  - Total compounds: {len(smiles_list)}")

        # Initialize statistics
        results = []

        # Iterate over each important bit
        for idx, bit_idx in enumerate(important_bits):
            print(f"\n[{idx+1}/{len(important_bits)}] Analyzing bit_{bit_idx}...")

            bit_name = bit_names[idx] if bit_names else f"bit_{bit_idx}"
            shap_val = shap_values[idx] if shap_values else None

            # Find molecules containing this bit
            molecules_with_bit = []

            for mol_idx, smiles in enumerate(smiles_list):
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue

                # Compute Morgan fingerprint and get bit info
                bit_info = {}
                fp = AllChem.GetMorganFingerprintAsBitVect(
                    mol,
                    radius=self.radius,
                    nBits=self.n_bits,
                    bitInfo=bit_info
                )

                # Check if this bit is activated
                if fp[bit_idx] == 1:
                    molecules_with_bit.append({
                        'mol_idx': mol_idx,
                        'smiles': smiles,
                        'mol': mol,
                        'atom_ids': bit_info.get(bit_idx, [])
                    })

            # Statistics
            n_mols_with_bit = len(molecules_with_bit)
            frequency = n_mols_with_bit / len(smiles_list) * 100

            print(f"  - Frequency: {n_mols_with_bit}/{len(smiles_list)} ({frequency:.2f}%)")

            # Append to results
            results.append({
                'bit_idx': bit_idx,
                'bit_name': bit_name,
                'shap_value': shap_val,
                'n_molecules': n_mols_with_bit,
                'frequency_pct': frequency,
                'example_molecules': molecules_with_bit[:10]  # Save up to 10 examples
            })

            # Store bit_info for later visualization
            self.bit_info_dict[bit_idx] = molecules_with_bit

        # Convert to DataFrame
        results_df = pd.DataFrame([
            {
                'bit_idx': r['bit_idx'],
                'bit_name': r['bit_name'],
                'shap_value': r['shap_value'],
                'n_molecules': r['n_molecules'],
                'frequency_pct': r['frequency_pct']
            }
            for r in results
        ])

        # Save statistics
        save_path = self.output_dir / 'morgan_bit_statistics.csv'
        results_df.to_csv(save_path, index=False)
        print(f"\n✓ Statistics saved to: {save_path}")

        return results_df

    def visualize_bit_substructures(
        self,
        bit_idx: int,
        n_examples: int = 6,
        mols_per_row: int = 3,
        mol_size: Tuple[int, int] = (300, 250)
    ):
        """
        Visualize the substructures corresponding to a specific bit.

        Args:
            bit_idx: Morgan bit index
            n_examples: Number of example molecules to display
            mols_per_row: Molecules per row
            mol_size: Size of each molecule image
        """
        if bit_idx not in self.bit_info_dict:
            raise ValueError(f"Bit {bit_idx} not analyzed yet. Call analyze_important_bits() first.")

        molecules_with_bit = self.bit_info_dict[bit_idx]

        if len(molecules_with_bit) == 0:
            print(f"⚠️ Warning: bit_{bit_idx} not found in any molecule")
            return

        # Select molecules to display
        examples = molecules_with_bit[:n_examples]

        print(f"\nVisualizing bit_{bit_idx} ({len(molecules_with_bit)} molecules contain this bit)")
        print(f"  - Showing first {len(examples)} examples")

        # Prepare plot
        n_rows = (len(examples) + mols_per_row - 1) // mols_per_row
        fig = plt.figure(figsize=(mols_per_row * 4, n_rows * 3.5))

        for idx, example in enumerate(examples):
            mol = example['mol']
            atom_ids = example['atom_ids']

            # Get the atomic environment for this bit
            # atom_ids is a list of (center_atom, radius) tuples
            if len(atom_ids) > 0:
                # Take the first matching center atom
                center_atom, radius_info = atom_ids[0]

                # Get all atoms within the specified radius of the center atom
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, self.radius, center_atom)

                # Convert to atom index list
                atom_indices = set()
                for bond_idx in env:
                    bond = mol.GetBondWithIdx(bond_idx)
                    atom_indices.add(bond.GetBeginAtomIdx())
                    atom_indices.add(bond.GetEndAtomIdx())
                atom_indices.add(center_atom)

                highlight_atoms = list(atom_indices)
                highlight_bonds = list(env)
            else:
                highlight_atoms = []
                highlight_bonds = []

            # Draw molecule with highlighted substructure
            ax = fig.add_subplot(n_rows, mols_per_row, idx + 1)

            drawer = Draw.MolDraw2DCairo(mol_size[0], mol_size[1])
            drawer.DrawMolecule(
                mol,
                highlightAtoms=highlight_atoms,
                highlightBonds=highlight_bonds
            )
            drawer.FinishDrawing()

            # Convert image to matplotlib-compatible format
            img = drawer.GetDrawingText()
            from PIL import Image
            import io
            pil_img = Image.open(io.BytesIO(img))

            ax.imshow(pil_img)
            ax.axis('off')
            ax.set_title(f'Example {idx+1}\nSMILES: {example["smiles"][:30]}...',
                        fontsize=9)

        plt.suptitle(f'Morgan Bit {bit_idx} - Substructure Examples',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()

        # Save figure
        save_path = self.output_dir / f'bit_{bit_idx}_substructures.png'
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"✓ Substructure visualization saved to: {save_path}")

    def visualize_all_important_bits(
        self,
        top_k: int = 10,
        n_examples_per_bit: int = 4
    ):
        """
        Batch-visualize substructures for all important bits.

        Args:
            top_k: Number of top bits to visualize
            n_examples_per_bit: Example molecules per bit
        """
        print(f"\nBatch-visualizing top-{top_k} important Morgan bits...")

        # Sort bit_info_dict by frequency
        sorted_bits = sorted(
            self.bit_info_dict.keys(),
            key=lambda b: len(self.bit_info_dict[b]),
            reverse=True
        )[:top_k]

        for idx, bit_idx in enumerate(sorted_bits):
            print(f"\n[{idx+1}/{len(sorted_bits)}] ", end="")
            self.visualize_bit_substructures(
                bit_idx=bit_idx,
                n_examples=n_examples_per_bit,
                mols_per_row=2
            )

        print(f"\n✓ All visualizations complete!")

    def generate_summary_report(
        self,
        statistics_df: pd.DataFrame,
        output_name: str = 'morgan_substructure_summary.md'
    ):
        """
        Generate a summary report for the Morgan substructure analysis.

        Args:
            statistics_df: Statistics results DataFrame
            output_name: Output report filename
        """
        print(f"\nGenerating analysis report...")

        report = []
        report.append("# Morgan Fingerprint Substructure Analysis Report")
        report.append("\n## Overview\n")
        report.append(f"- **Morgan radius**: {self.radius}")
        report.append(f"- **Fingerprint length**: {self.n_bits}")
        report.append(f"- **Bits analyzed**: {len(statistics_df)}")
        report.append("")

        report.append("\n## Important Morgan Bits — Statistics\n")
        report.append("| Rank | Bit index | Bit name | SHAP value | # Molecules | Frequency (%) |")
        report.append("|------|-----------|----------|------------|-------------|---------------|")

        for idx, row in statistics_df.iterrows():
            rank = idx + 1
            shap_str = f"{row['shap_value']:.4f}" if pd.notna(row['shap_value']) else "N/A"
            report.append(
                f"| {rank} | {row['bit_idx']} | {row['bit_name']} | "
                f"{shap_str} | {row['n_molecules']} | {row['frequency_pct']:.2f} |"
            )

        report.append("\n## Substructure Visualizations\n")
        report.append("Example substructures for each important bit are saved as individual PNG files:\n")

        for idx, row in statistics_df.head(10).iterrows():
            bit_idx = row['bit_idx']
            img_path = f"bit_{bit_idx}_substructures.png"
            if (self.output_dir / img_path).exists():
                report.append(f"- **Bit {bit_idx}**: [{img_path}]({img_path})")

        report.append("\n## Chemical Interpretation\n")
        report.append("### High-frequency Substructure Patterns\n")

        high_freq = statistics_df[statistics_df['frequency_pct'] > 50]
        if len(high_freq) > 0:
            report.append("\n**High-frequency substructures** (frequency > 50%):\n")
            for _, row in high_freq.iterrows():
                report.append(f"- **Bit {row['bit_idx']}**: present in {row['frequency_pct']:.1f}% of molecules")

        low_freq_high_shap = statistics_df[
            (statistics_df['frequency_pct'] < 20) &
            (statistics_df['shap_value'] > 0.5)
        ]
        if len(low_freq_high_shap) > 0:
            report.append("\n**Rare but highly important substructures** (frequency < 20% and SHAP > 0.5):\n")
            for _, row in low_freq_high_shap.iterrows():
                report.append(
                    f"- **Bit {row['bit_idx']}**: frequency {row['frequency_pct']:.1f}%, "
                    f"SHAP={row['shap_value']:.4f}"
                )

        report.append("\n## Suggested Follow-up Analyses\n")
        report.append("1. Provide chemical interpretation for rare substructures with high SHAP values")
        report.append("2. Analyze the spatial relationship between these substructures and the boron atom")
        report.append("3. Investigate the directional effect (positive/negative) of each substructure on chemical shift")

        # Save report
        report_path = self.output_dir / output_name
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report))

        print(f"✓ Analysis report saved to: {report_path}")


def main():
    """Test and demonstration code."""
    # Example: analyze a few boron-containing compounds
    example_smiles = [
        "B(O)(O)c1ccccc1",            # phenylboronic acid
        "c1ccc(B2OC(C)(C)CO2)cc1",    # phenylboronic acid pinacol ester
        "Cc1ccc(B(O)O)cc1",           # 4-methylphenylboronic acid
    ]

    analyzer = MorganSubstructureAnalyzer(
        output_dir="outputs/morgan_analysis",
        radius=2,
        n_bits=1024
    )

    # Analyze bits 690, 271, 390
    important_bits = [690, 271, 390]

    stats_df = analyzer.analyze_important_bits(
        smiles_list=example_smiles,
        important_bits=important_bits,
        bit_names=[f"morgan_bit_{b}" for b in important_bits]
    )

    print("\nStatistics:")
    print(stats_df)

    # Visualize
    for bit_idx in important_bits:
        analyzer.visualize_bit_substructures(bit_idx, n_examples=3)


if __name__ == "__main__":
    main()
