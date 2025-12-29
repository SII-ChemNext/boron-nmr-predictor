import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import TransformerConv, BatchNorm

class BoronNMRNet(nn.Module):
    def __init__(self, node_in_dim, edge_in_dim, solvent_dim, hidden_dim=128, dropout=0.2):
        """
        参数说明:
        node_in_dim:  节点特征输入维度 (来自 features.py)
        edge_in_dim:  边特征输入维度 (来自 features.py)
        solvent_dim:  溶剂指纹维度 (1024)
        hidden_dim:   GNN 内部隐藏层宽度 (支持 Optuna 调参)
        dropout:      Dropout 概率 (支持 Optuna 调参)
        """
        super(BoronNMRNet, self).__init__()
        
        self.dropout_rate = dropout
        
        # 1. 特征编码层 (Embedding)
        # 将不同维度的原始特征统一映射到 hidden_dim
        self.node_encoder = nn.Linear(node_in_dim, hidden_dim)
        self.edge_encoder = nn.Linear(edge_in_dim, hidden_dim)
        
        # 2. GNN 主干 (Backbone)
        # 使用 TransformerConv，因为它擅长捕捉长距离依赖 (Long-range dependency)
        # heads=4, concat=False 意味着多头注意力融合后维度保持不变
        self.conv1 = TransformerConv(hidden_dim, hidden_dim, heads=4, edge_dim=hidden_dim, concat=False)
        self.bn1 = BatchNorm(hidden_dim)
        
        self.conv2 = TransformerConv(hidden_dim, hidden_dim, heads=4, edge_dim=hidden_dim, concat=False)
        self.bn2 = BatchNorm(hidden_dim)
        
        self.conv3 = TransformerConv(hidden_dim, hidden_dim, heads=4, edge_dim=hidden_dim, concat=False)
        self.bn3 = BatchNorm(hidden_dim)
        
        # 3. 融合与预测头 (Prediction Head)
        # 最终特征 = GNN提取的原子特征 + 溶剂的全局指纹
        fusion_dim = hidden_dim + solvent_dim
        
        self.mlp = nn.Sequential(
            nn.Linear(fusion_dim, 256),
            nn.ReLU(),
            nn.Dropout(self.dropout_rate), # 这里使用传入的 dropout
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, 1) # 输出最终的 ppm 值
        )

    def forward(self, x, edge_index, edge_attr, solvent_x, mask_b, batch_index):
        """
        前向传播函数
        """
        # --- A. 初始编码 ---
        x = self.node_encoder(x)
        edge_attr = self.edge_encoder(edge_attr)
        
        # --- B. GNN 消息传递 (带残差连接) ---
        
        # 第一层
        x_in = x # 保存输入用于残差
        x = self.conv1(x, edge_index, edge_attr)
        x = self.bn1(x)
        x = F.relu(x)
        x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x = x + x_in # 残差连接 1
        
        # 第二层
        x_in = x
        x = self.conv2(x, edge_index, edge_attr)
        x = self.bn2(x)
        x = F.relu(x)
        x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x = x + x_in # 残差连接 2
        
        # 第三层
        x_in = x
        x = self.conv3(x, edge_index, edge_attr)
        x = self.bn3(x)
        x = F.relu(x)
        x = x + x_in # 残差连接 3
        
        # --- C. 提取目标原子 (Boron) ---
        # 只取 mask_b 为 True 的节点
        b_features = x[mask_b] 
        
        # --- D. 融合溶剂信息 ---
        # 1. 找出这些 B 原子分别属于 Batch 中的哪个分子
        # batch_index 的形状是 [Num_Nodes]，只取 B 原子的部分
        b_batch_idx = batch_index[mask_b]
        
        # 2. 根据索引，将对应的溶剂特征“广播”给 B 原子
        # solvent_x: [Batch_Size, 1024] -> b_solvent_features: [Num_B_Atoms, 1024]
        b_solvent_features = solvent_x[b_batch_idx]
        
        # 3. 拼接 (Concatenate)
        combined = torch.cat([b_features, b_solvent_features], dim=1)
        
        # --- E. MLP 预测 ---
        out = self.mlp(combined)
        
        return out.squeeze(-1) # 压缩维度，输出 [Num_B_Atoms]