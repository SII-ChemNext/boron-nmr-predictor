"""
B11 NMR Interpretability Analysis - Core Utility Module
========================================================
Data loading, feature mapping, and model loading utilities.
"""

import os
import json
import yaml
import numpy as np
import pandas as pd
import xgboost as xgb
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


class FeatureMapper:
    """Feature mapping utility — manages the 2482-dimensional feature blocks and naming."""

    def __init__(self, config: Dict):
        """
        Initialize the feature mapper.

        Args:
            config: Feature configuration dict containing total_dims and blocks definition.
        """
        self.config = config
        self.total_dims = config['total_dims']
        self.blocks = config['blocks']

        # Generate feature name mapping
        self.feature_names = self._generate_feature_names()

    def _generate_feature_names(self) -> List[str]:
        """Generate the name list for all features."""
        names = []

        # Compound Morgan fingerprint
        compound_morgan = self.blocks['compound_morgan']
        for i in range(compound_morgan['start'], compound_morgan['end']):
            names.append(f"compound_morgan_bit_{i}")

        # Compound RDKit descriptors
        compound_rdkit = self.blocks['compound_rdkit']
        rdkit_desc_names = self._get_rdkit_descriptor_names()
        for desc_name in rdkit_desc_names:
            names.append(f"compound_{desc_name}")

        # Solvent Morgan fingerprint
        solvent_morgan = self.blocks['solvent_morgan']
        for i in range(solvent_morgan['end'] - solvent_morgan['start']):
            names.append(f"solvent_morgan_bit_{i}")

        # Solvent RDKit descriptors
        for desc_name in rdkit_desc_names:
            names.append(f"solvent_{desc_name}")

        return names

    def _get_rdkit_descriptor_names(self) -> List[str]:
        """Return the list of RDKit descriptor names (217 descriptors)."""
        descriptor_list = []
        for name, func in Descriptors.descList:
            descriptor_list.append(name)
        return descriptor_list

    def get_block_indices(self, block_name: str) -> Tuple[int, int]:
        """
        Get the index range for a feature block.

        Args:
            block_name: Block name ('compound_morgan', 'compound_rdkit', etc.)

        Returns:
            (start_idx, end_idx) tuple
        """
        if block_name not in self.blocks:
            raise ValueError(f"Unknown block name: {block_name}")

        block = self.blocks[block_name]
        return block['start'], block['end']

    def get_block_features(self, X: np.ndarray, block_name: str) -> np.ndarray:
        """
        Extract features for a specific block.

        Args:
            X: Full feature matrix (n_samples, 2482)
            block_name: Block name

        Returns:
            Feature sub-matrix for the specified block
        """
        start, end = self.get_block_indices(block_name)
        return X[:, start:end]

    def get_feature_name(self, idx: int) -> str:
        """Return the feature name for a given index."""
        if idx < 0 or idx >= self.total_dims:
            raise ValueError(f"Feature index {idx} out of range [0, {self.total_dims})")
        return self.feature_names[idx]

    def get_block_name(self, idx: int) -> str:
        """Return the block name for a given feature index."""
        for block_name, block_info in self.blocks.items():
            if block_info['start'] <= idx < block_info['end']:
                return block_name
        raise ValueError(f"Feature index {idx} not in any block")

    def summarize_shap_by_block(self, shap_values: np.ndarray) -> Dict[str, float]:
        """
        Aggregate SHAP values by feature block.

        Args:
            shap_values: SHAP value matrix (n_samples, n_features)

        Returns:
            Mean absolute SHAP value per block
        """
        summary = {}

        for block_name, block_info in self.blocks.items():
            start, end = block_info['start'], block_info['end']
            block_shap = shap_values[:, start:end]
            summary[block_name] = np.abs(block_shap).mean()

        return summary


class DataLoader:
    """Data loader — handles loading of features, labels, and SMILES."""

    def __init__(self, paths_config: Dict):
        """
        Initialize the data loader.

        Args:
            paths_config: Path configuration dict
        """
        self.paths = paths_config
        self.base_dir = Path(__file__).parent.parent

    def _resolve_path(self, relative_path: str) -> Path:
        """Resolve a relative path to an absolute path."""
        return (self.base_dir / relative_path).resolve()

    def load_features_and_labels(
        self, split: str = 'train'
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Load features and labels.

        Args:
            split: Dataset split ('train', 'val', 'test')

        Returns:
            (X, y) tuple
        """
        features_path = self._resolve_path(self.paths[f'{split}_features'])
        labels_path = self._resolve_path(self.paths[f'{split}_labels'])

        X = pd.read_csv(features_path).values
        y = pd.read_csv(labels_path).values.ravel()

        return X, y

    def load_smiles(self, split: str = 'train') -> pd.DataFrame:
        """
        Load raw SMILES data.

        Args:
            split: Dataset split

        Returns:
            DataFrame with Smiles and solvent columns
        """
        smiles_path = self._resolve_path(self.paths[f'{split}_smiles'])
        return pd.read_csv(smiles_path)

    def load_all_data(
        self, split: str = 'train'
    ) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
        """
        Load features, labels, and SMILES in one call.

        Args:
            split: Dataset split

        Returns:
            (X, y, smiles_df) tuple
        """
        X, y = self.load_features_and_labels(split)
        smiles_df = self.load_smiles(split)

        return X, y, smiles_df


class ModelLoader:
    """Model loader."""

    @staticmethod
    def load_xgboost_model(model_path: str) -> xgb.Booster:
        """
        Load an XGBoost model.

        Args:
            model_path: Path to the model file (.json format)

        Returns:
            XGBoost Booster object
        """
        booster = xgb.Booster()
        booster.load_model(model_path)
        return booster


class MorganBitAnalyzer:
    """Morgan fingerprint bit analyzer — extracts and visualizes important substructures."""

    def __init__(self, radius: int = 2, n_bits: int = 1024):
        """
        Initialize.

        Args:
            radius: Morgan fingerprint radius
            n_bits: Fingerprint length
        """
        self.radius = radius
        self.n_bits = n_bits

    def get_bit_info(self, smiles: str) -> Dict[int, Tuple]:
        """
        Get Morgan fingerprint bit information for a SMILES string.

        Args:
            smiles: Molecular SMILES string

        Returns:
            Dict {bit_idx: [(atom_idx, radius), ...]}
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}

        bit_info = {}
        _ = AllChem.GetMorganFingerprintAsBitVect(
            mol, self.radius, nBits=self.n_bits, bitInfo=bit_info
        )

        return bit_info

    def get_substructure_for_bit(
        self, smiles: str, bit_idx: int
    ) -> Optional[Tuple[Chem.Mol, List[int]]]:
        """
        Get the substructure that activates a specific bit.

        Args:
            smiles: Molecular SMILES
            bit_idx: Morgan bit index

        Returns:
            (mol, atom_list) or None;
            atom_list contains the atom indices that activate this bit
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        bit_info = self.get_bit_info(smiles)

        if bit_idx not in bit_info:
            return None

        # bit_info[bit_idx] is a list of [(atom_idx, radius), ...]
        # Take the first substructure environment
        center_atom, _ = bit_info[bit_idx][0]

        # Get all atoms within the environment
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, self.radius, center_atom)
        atom_list = list(set([center_atom] + [bond for bond in env]))

        return mol, atom_list

    def batch_analyze_bits(
        self, smiles_list: List[str], important_bits: List[int]
    ) -> Dict[int, List[str]]:
        """
        Batch analysis of important bits across multiple molecules.

        Args:
            smiles_list: List of molecular SMILES
            important_bits: List of important bit indices

        Returns:
            Dict {bit_idx: [smiles_containing_this_bit]}
        """
        bit_to_smiles = {bit: [] for bit in important_bits}

        for smiles in smiles_list:
            bit_info = self.get_bit_info(smiles)
            for bit in important_bits:
                if bit in bit_info:
                    bit_to_smiles[bit].append(smiles)

        return bit_to_smiles


def load_config(config_path: str = None) -> Dict:
    """
    Load the configuration file.

    Args:
        config_path: Path to config file; uses default path if None.

    Returns:
        Configuration dict
    """
    if config_path is None:
        config_path = Path(__file__).parent.parent / 'configs' / 'analysis_config.yaml'

    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)

    return config


def ensure_output_dir(output_dir: str):
    """Ensure the output directory exists."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)


def save_results(results: Dict, output_path: str):
    """Save analysis results as JSON."""
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    # Quick test
    print("Testing FeatureMapper...")

    config = load_config()
    mapper = FeatureMapper(config['features'])

    print(f"Total features: {mapper.total_dims}")
    print(f"Feature blocks: {list(mapper.blocks.keys())}")
    print(f"First 5 feature names: {mapper.feature_names[:5]}")
    print(f"Feature 1024 name: {mapper.get_feature_name(1024)}")
    print(f"Feature 1024 block: {mapper.get_block_name(1024)}")

    print("\nTesting MorganBitAnalyzer...")
    analyzer = MorganBitAnalyzer()
    test_smiles = "CC(C)C(B1Nc2cccc3cccc(c23)N1)B1OC(C)(C)C(C)(C)O1"
    bit_info = analyzer.get_bit_info(test_smiles)
    print(f"Morgan bits activated: {len(bit_info)}")
