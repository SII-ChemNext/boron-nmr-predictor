"""11B NMR Predictor - encapsulates V3 model prediction logic"""

import torch
import os
import uuid
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from torch_geometric.data import Data, Batch

from .model import BoronNMRNet_V3
from .features import get_atom_features, get_bond_features, SOLVENT_FP_SIZE
from .ml_features import compute_and_normalize, load_scaler, N_ML_FEATURES
from utils.exceptions import PredictionError

# Solvent name -> solvent ID mapping (consistent with training)
SOLVENT_NAME_TO_ID = {
    'CDCl3': 0,
    'C6D6': 1,
    'd6-DMSO': 2,
    'CD3COCD3': 3,
    'CD3CN': 4,
    'CD3OD': 5,
    'CD2Cl2': 6,
    'd8-THF': 7,
    'd8-Toluene': 8,
    'D2O': 9,
}
UNKNOWN_SOLVENT_ID = 10


class BoronNMRPredictor:
    """
    11B NMR Chemical Shift Predictor

    Encapsulates the loading and ensemble prediction logic for 5-fold V3 models.
    """

    def __init__(self, model_dir, device='cpu', hidden_dim=256,
                 dropout=0.012558398103042557, solvent_dim=32,
                 ml_feature_dim=20, ml_hidden_dim=64):
        """
        Initialize the predictor.

        Args:
            model_dir (str): directory containing model files
            device (str): compute device ('cpu' or 'cuda')
            hidden_dim (int): GNN hidden layer dimension
            dropout (float): dropout probability
            solvent_dim (int): solvent embedding dimension
            ml_feature_dim (int): ML global feature dimension
            ml_hidden_dim (int): encoded dimension of ML features
        """
        self.device = torch.device(device)
        self.hidden_dim = hidden_dim
        self.dropout = dropout
        self.solvent_dim = solvent_dim
        self.ml_feature_dim = ml_feature_dim
        self.ml_hidden_dim = ml_hidden_dim
        self.models = []

        # Automatically infer feature dimensions
        self.node_dim, self.edge_dim = self._get_feature_dims()

        # Load ML feature scaler
        scaler_path = os.path.join(model_dir, 'ml_feature_scaler.pkl')
        if os.path.exists(scaler_path):
            self.ml_scaler, self.ml_medians = load_scaler(scaler_path)
            print(f"  ✓ ML feature scaler loaded: {scaler_path}")
        else:
            self.ml_scaler = None
            self.ml_medians = None
            print(f"  ⚠ ML feature scaler not found, ML global features will not be used")

        # Load 5-fold models
        self.models = self._load_models(model_dir)
        print(f"✓ Successfully loaded {len(self.models)} models on device: {device}")

    def _get_feature_dims(self):
        """Automatically infer node and edge feature dimensions"""
        try:
            m_tmp = Chem.MolFromSmiles('CB')
            AllChem.ComputeGasteigerCharges(m_tmp)
            node_dim = get_atom_features(m_tmp.GetAtoms()[0]).shape[0]
            edge_dim = get_bond_features(m_tmp.GetBonds()[0]).shape[0]
            return node_dim, edge_dim
        except Exception as e:
            print(f"⚠ Warning: cannot infer feature dimensions, using defaults - {e}")
            return 58, 10  # V3 default values

    def _load_models(self, model_dir):
        """Load 5-fold V3 models"""
        models = []

        # Determine whether to use ML features
        effective_ml_dim = self.ml_feature_dim if self.ml_scaler is not None else 0

        for fold_idx in range(1, 6):
            model_path = os.path.join(model_dir, f'model_v3_fold_{fold_idx}.pth')

            if not os.path.exists(model_path):
                raise FileNotFoundError(f"Model file not found: {model_path}")

            try:
                model = BoronNMRNet_V3(
                    node_in_dim=self.node_dim,
                    edge_in_dim=self.edge_dim,
                    num_solvents=11,
                    solvent_dim=self.solvent_dim,
                    hidden_dim=self.hidden_dim,
                    dropout=self.dropout,
                    num_heads=4,
                    ml_feature_dim=effective_ml_dim,
                    ml_hidden_dim=self.ml_hidden_dim
                ).to(self.device)

                # Load weights
                try:
                    model.load_state_dict(
                        torch.load(model_path, map_location=self.device, weights_only=True)
                    )
                except TypeError:
                    model.load_state_dict(
                        torch.load(model_path, map_location=self.device)
                    )

                model.eval()
                models.append(model)
                print(f"  ✓ Loaded model_v3_fold_{fold_idx}.pth")

            except Exception as e:
                raise Exception(f"Failed to load model_v3_fold_{fold_idx}.pth: {e}")

        return models

    def _get_solvent_id(self, solvent_name):
        """Map solvent name to integer ID"""
        return SOLVENT_NAME_TO_ID.get(solvent_name, UNKNOWN_SOLVENT_ID)

    def _smiles_to_data(self, mol_smiles, solvent_name):
        """
        Convert SMILES to a PyG Data object (V3 format).

        Args:
            mol_smiles (str): molecule SMILES
            solvent_name (str): solvent name

        Returns:
            tuple: (data, mol, canonical_smiles)
        """
        try:
            # ========== STEP 1: Molecule canonicalization ==========
            temp_mol = Chem.MolFromSmiles(mol_smiles)
            if temp_mol is None:
                raise PredictionError(f"Invalid molecule SMILES: {mol_smiles}")

            canonical_mol_smiles = Chem.MolToSmiles(
                temp_mol, canonical=True, isomericSmiles=True
            )
            mol = Chem.MolFromSmiles(canonical_mol_smiles)

            # ========== STEP 2: Compute Gasteiger charges ==========
            AllChem.ComputeGasteigerCharges(mol)

            # ========== STEP 3: Extract node features ==========
            atom_feats = []
            boron_mask = []

            for atom in mol.GetAtoms():
                atom_feats.append(get_atom_features(atom))
                boron_mask.append(atom.GetSymbol() == 'B')

            x = torch.stack(atom_feats)
            mask_b = torch.tensor(boron_mask, dtype=torch.bool)

            if mask_b.sum() == 0:
                raise PredictionError("No boron atom found in molecule!")

            # ========== STEP 4: Extract edge features ==========
            edge_indices = []
            edge_attrs = []

            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()

                edge_indices.append([i, j])
                edge_indices.append([j, i])

                b_feat = get_bond_features(bond)
                edge_attrs.append(b_feat)
                edge_attrs.append(b_feat)

            if edge_indices:
                edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
                edge_attr = torch.stack(edge_attrs)
            else:
                edge_index = torch.empty((2, 0), dtype=torch.long)
                edge_attr = torch.empty((0, self.edge_dim), dtype=torch.float)

            # ========== STEP 5: Solvent ID ==========
            solvent_id = torch.tensor([self._get_solvent_id(solvent_name)], dtype=torch.long)

            # ========== STEP 6: ML global features ==========
            if self.ml_scaler is not None:
                ml_feats = compute_and_normalize(canonical_mol_smiles, self.ml_scaler, self.ml_medians)
                ml_global_features = torch.tensor(ml_feats, dtype=torch.float).unsqueeze(0)  # [1, 20]
            else:
                ml_global_features = None

            # ========== STEP 7: Pack PyG Data object ==========
            data = Data(
                x=x,
                edge_index=edge_index,
                edge_attr=edge_attr,
                mask_b=mask_b,
                solvent_id=solvent_id,
            )

            if ml_global_features is not None:
                data.ml_global_features = ml_global_features

            return data, mol, canonical_mol_smiles

        except PredictionError:
            raise
        except Exception as e:
            raise PredictionError(f"Data processing failed: {e}")

    def predict(self, mol_smiles, solvent_name):
        """
        Run ensemble prediction.

        Args:
            mol_smiles (str): molecule SMILES
            solvent_name (str): solvent name (e.g. 'CDCl3')

        Returns:
            dict: prediction results
        """
        try:
            # ========== STEP 1: Data preprocessing ==========
            data, mol, canonical_smiles = self._smiles_to_data(mol_smiles, solvent_name)
            data = data.to(self.device)

            # ========== STEP 2: Simulate batch ==========
            batch = Batch.from_data_list([data])

            # ========== STEP 3: Ensemble inference (5 models) ==========
            fold_preds = []

            with torch.no_grad():
                for i, model in enumerate(self.models):
                    try:
                        # Prepare ML global features
                        ml_feats = getattr(batch, 'ml_global_features', None)

                        pred = model(
                            batch.x,
                            batch.edge_index,
                            batch.edge_attr,
                            batch.solvent_id,
                            batch.mask_b,
                            batch.batch,
                            ml_global_features=ml_feats
                        )

                        if pred.ndim == 0:
                            pred = pred.view(1)

                        fold_preds.append(pred)

                    except Exception as e:
                        raise PredictionError(f"Fold {i+1} prediction failed: {e}")

            # ========== STEP 4: Average predictions across 5 models ==========
            if not fold_preds:
                raise PredictionError("No successful model predictions")

            avg_pred = torch.stack(fold_preds).mean(dim=0)

            if avg_pred.ndim == 0:
                avg_pred = avg_pred.view(1)

            # ========== STEP 5: Map predictions to atom indices ==========
            boron_indices = torch.where(data.mask_b)[0].cpu().numpy()
            predictions_array = avg_pred.cpu().numpy()

            predictions = []
            for idx, ppm in zip(boron_indices, predictions_array):
                atom = mol.GetAtomWithIdx(int(idx))
                predictions.append({
                    'atom_index': int(idx),
                    'element': atom.GetSymbol(),
                    'ppm': float(ppm)
                })

            return {
                'canonical_smiles': canonical_smiles,
                'predictions': predictions,
                'num_borons': len(predictions),
                'mol_object': mol
            }

        except PredictionError:
            raise
        except Exception as e:
            raise PredictionError(f"Prediction process error: {e}")

    def generate_molecule_image(self, mol, boron_indices, predictions, output_dir):
        """Generate molecule image with highlighted boron atoms"""
        try:
            os.makedirs(output_dir, exist_ok=True)

            filename = f"prediction_{uuid.uuid4().hex[:8]}.png"
            filepath = os.path.join(output_dir, filename)

            highlight_atoms = list(boron_indices)
            atom_colors = {idx: (0.8, 0.9, 1.0) for idx in highlight_atoms}

            from rdkit.Chem.Draw import rdMolDraw2D

            drawer = rdMolDraw2D.MolDraw2DCairo(800, 600)

            opts = drawer.drawOptions()
            for pred in predictions:
                idx = pred['atom_index']
                opts.atomLabels[idx] = str(idx)

            drawer.DrawMolecule(
                mol,
                highlightAtoms=highlight_atoms,
                highlightAtomColors=atom_colors
            )
            drawer.FinishDrawing()

            with open(filepath, 'wb') as f:
                f.write(drawer.GetDrawingText())

            return filepath

        except Exception as e:
            raise Exception(f"Molecule image generation failed: {e}")
