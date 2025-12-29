"""硼核磁预测器 - 封装预测逻辑"""

import torch
import os
import uuid
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from torch_geometric.data import Data, Batch
import numpy as np

from .model import BoronNMRNet
from .features import get_atom_features, get_bond_features, get_solvent_features, SOLVENT_FP_SIZE
from utils.exceptions import PredictionError


class BoronNMRPredictor:
    """
    硼核磁预测器 - 单例模式

    封装了 5 个 fold 模型的加载和集成预测逻辑
    """

    def __init__(self, model_dir, device='cpu', hidden_dim=256, dropout=0.0472):
        """
        初始化预测器

        Args:
            model_dir (str): 模型文件目录
            device (str): 计算设备（'cpu' 或 'cuda'）
            hidden_dim (int): GNN 隐藏层维度
            dropout (float): Dropout 概率
        """
        self.device = torch.device(device)
        self.hidden_dim = hidden_dim
        self.dropout = dropout
        self.models = []

        # 自动推断特征维度
        self.node_dim, self.edge_dim = self._get_feature_dims()

        # 加载 5 个 fold 模型
        self.models = self._load_models(model_dir)
        print(f"✓ 成功加载 {len(self.models)} 个模型，使用设备: {device}")

    def _get_feature_dims(self):
        """
        自动推断节点和边特征维度

        Returns:
            tuple: (node_dim, edge_dim)
        """
        try:
            m_tmp = Chem.MolFromSmiles('CB')
            AllChem.ComputeGasteigerCharges(m_tmp)
            node_dim = get_atom_features(m_tmp.GetAtoms()[0]).shape[0]
            edge_dim = get_bond_features(m_tmp.GetBonds()[0]).shape[0]
            return node_dim, edge_dim
        except Exception as e:
            print(f"⚠ 警告: 无法自动推断特征维度，使用默认值 - {e}")
            return 134, 14  # 默认值

    def _load_models(self, model_dir):
        """
        加载 5 个 fold 模型

        Args:
            model_dir (str): 模型目录

        Returns:
            list: 5 个模型对象列表

        Raises:
            FileNotFoundError: 模型文件不存在
            Exception: 模型加载失败
        """
        models = []

        for fold_idx in range(1, 6):
            model_path = os.path.join(model_dir, f'model_fold_{fold_idx}.pth')

            if not os.path.exists(model_path):
                raise FileNotFoundError(f"模型文件不存在: {model_path}")

            try:
                model = BoronNMRNet(
                    node_in_dim=self.node_dim,
                    edge_in_dim=self.edge_dim,
                    solvent_dim=SOLVENT_FP_SIZE,
                    hidden_dim=self.hidden_dim,
                    dropout=self.dropout
                ).to(self.device)

                # 加载权重，支持新旧 PyTorch 版本
                try:
                    model.load_state_dict(
                        torch.load(model_path, map_location=self.device, weights_only=True)
                    )
                except TypeError:
                    # 旧版本 PyTorch 不支持 weights_only 参数
                    model.load_state_dict(
                        torch.load(model_path, map_location=self.device)
                    )

                model.eval()
                models.append(model)
                print(f"  ✓ 加载 model_fold_{fold_idx}.pth")

            except Exception as e:
                raise Exception(f"加载 model_fold_{fold_idx}.pth 失败: {e}")

        return models

    def _smiles_to_data(self, mol_smiles, solvent_smiles):
        """
        将 SMILES 转换为 PyG Data 对象

        Args:
            mol_smiles (str): 分子 SMILES
            solvent_smiles (str): 溶剂 SMILES

        Returns:
            tuple: (data, mol, canonical_smiles)

        Raises:
            PredictionError: 处理失败
        """
        try:
            # ========== STEP 1: 分子标准化 ==========
            temp_mol = Chem.MolFromSmiles(mol_smiles)
            if temp_mol is None:
                raise PredictionError(f"无效的分子 SMILES: {mol_smiles}")

            # 转换为标准 SMILES（保留立体化学）
            canonical_mol_smiles = Chem.MolToSmiles(
                temp_mol, canonical=True, isomericSmiles=True
            )

            # 基于标准 SMILES 重建分子
            mol = Chem.MolFromSmiles(canonical_mol_smiles)

            # ========== STEP 2: 溶剂标准化 ==========
            temp_solv = Chem.MolFromSmiles(solvent_smiles)
            if temp_solv:
                solvent_smiles = Chem.MolToSmiles(temp_solv, canonical=True)

            # ========== STEP 3: 计算 Gasteiger 电荷 ==========
            # 【关键步骤】必须在特征提取前计算！
            AllChem.ComputeGasteigerCharges(mol)

            # ========== STEP 4: 提取节点特征 ==========
            atom_feats = []
            boron_mask = []

            for atom in mol.GetAtoms():
                atom_feats.append(get_atom_features(atom))
                boron_mask.append(atom.GetSymbol() == 'B')

            x = torch.stack(atom_feats)  # [N_atoms, node_in_dim]
            mask_b = torch.tensor(boron_mask, dtype=torch.bool)

            if mask_b.sum() == 0:
                raise PredictionError("分子中不含硼原子！")

            # ========== STEP 5: 提取边特征 ==========
            edge_indices = []
            edge_attrs = []

            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()

                # 无向图：双向添加
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

            # ========== STEP 6: 提取溶剂特征 ==========
            solvent_feat = get_solvent_features(solvent_smiles).unsqueeze(0)

            # ========== STEP 7: 打包 PyG Data 对象 ==========
            data = Data(
                x=x,
                edge_index=edge_index,
                edge_attr=edge_attr,
                mask_b=mask_b,
                solvent_x=solvent_feat
            )

            return data, mol, canonical_mol_smiles

        except PredictionError:
            raise
        except Exception as e:
            raise PredictionError(f"数据处理失败: {e}")

    def predict(self, mol_smiles, solvent_smiles):
        """
        执行集成预测

        Args:
            mol_smiles (str): 分子 SMILES
            solvent_smiles (str): 溶剂 SMILES

        Returns:
            dict: {
                'canonical_smiles': str,
                'predictions': [
                    {'atom_index': int, 'element': str, 'ppm': float},
                    ...
                ],
                'num_borons': int,
                'mol_object': Mol
            }

        Raises:
            PredictionError: 预测失败
        """
        try:
            # ========== STEP 1: 数据预处理 ==========
            data, mol, canonical_smiles = self._smiles_to_data(mol_smiles, solvent_smiles)
            data = data.to(self.device)

            # ========== STEP 2: 模拟 Batch（PyG 要求） ==========
            batch = Batch.from_data_list([data])

            # ========== STEP 3: 集成推理（5 个模型） ==========
            fold_preds = []

            with torch.no_grad():
                for i, model in enumerate(self.models):
                    try:
                        pred = model(
                            batch.x,
                            batch.edge_index,
                            batch.edge_attr,
                            batch.solvent_x,
                            batch.mask_b,
                            batch.batch
                        )

                        # 确保输出为 1D 向量
                        if pred.ndim == 0:
                            pred = pred.view(1)

                        fold_preds.append(pred)

                    except Exception as e:
                        raise PredictionError(f"Fold {i+1} 预测失败: {e}")

            # ========== STEP 4: 取 5 个模型的平均值 ==========
            if not fold_preds:
                raise PredictionError("没有成功的模型预测")

            avg_pred = torch.stack(fold_preds).mean(dim=0)

            if avg_pred.ndim == 0:
                avg_pred = avg_pred.view(1)

            # ========== STEP 5: 映射预测值到原子索引 ==========
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
            raise PredictionError(f"预测过程异常: {e}")

    def generate_molecule_image(self, mol, boron_indices, predictions, output_dir):
        """
        生成高亮硼原子的分子图

        Args:
            mol (Mol): RDKit 分子对象
            boron_indices (list): 硼原子索引列表
            predictions (list): 预测结果列表
            output_dir (str): 输出目录

        Returns:
            str: 图片文件路径

        Raises:
            Exception: 图片生成失败
        """
        try:
            os.makedirs(output_dir, exist_ok=True)

            # 生成唯一的文件名
            filename = f"prediction_{uuid.uuid4().hex[:8]}.png"
            filepath = os.path.join(output_dir, filename)

            # 高亮原子列表和颜色
            highlight_atoms = list(boron_indices)
            atom_colors = {idx: (0.8, 0.9, 1.0) for idx in highlight_atoms}  # 浅蓝色

            # 使用 MolDraw2DCairo 绘制
            from rdkit.Chem.Draw import rdMolDraw2D

            drawer = rdMolDraw2D.MolDraw2DCairo(800, 600)

            # 设置原子标签
            opts = drawer.drawOptions()
            for pred in predictions:
                idx = pred['atom_index']
                # 仅显示原子序号（简洁清晰）
                opts.atomLabels[idx] = str(idx)

            drawer.DrawMolecule(
                mol,
                highlightAtoms=highlight_atoms,
                highlightAtomColors=atom_colors
            )
            drawer.FinishDrawing()

            # 保存图片
            with open(filepath, 'wb') as f:
                f.write(drawer.GetDrawingText())

            return filepath

        except Exception as e:
            raise Exception(f"分子图生成失败: {e}")
