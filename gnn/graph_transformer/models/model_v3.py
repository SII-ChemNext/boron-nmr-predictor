import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import TransformerConv, BatchNorm, global_add_pool

class BoronNMRNet_V3(nn.Module):
    """
    Improved 11B NMR prediction model (multi-head concatenation version)

    Improvements:
    1. Virtual Node mechanism
    2. Learnable solvent embedding (64-dim)
    3. Support for attention weight extraction (interpretability)
    4. [Modified] TransformerConv uses concat=True (multi-head concatenation)
    """

    def __init__(self, node_in_dim, edge_in_dim,
                 num_solvents=11, solvent_dim=64,
                 hidden_dim=256, dropout=0.05, num_heads=4,
                 ml_feature_dim=0, ml_hidden_dim=64):
        super().__init__()

        self.dropout_rate = dropout
        self.hidden_dim = hidden_dim
        self.num_heads = num_heads
        self.head_dim = hidden_dim // num_heads  # 256 / 4 = 64

        # ========================================
        # 1. Encoding layers
        # ========================================
        self.node_encoder = nn.Linear(node_in_dim, hidden_dim)
        self.edge_encoder = nn.Linear(edge_in_dim, hidden_dim)

        # Solvent embedding layer
        self.solvent_embedding = nn.Embedding(num_solvents, solvent_dim)

        # Virtual node embedding
        self.virtual_node_encoder = nn.Embedding(1, hidden_dim)

        # ========================================
        # 2. GNN backbone + virtual node update layers
        # concat=True: each head outputs head_dim, concatenated to hidden_dim
        # ========================================
        self.conv1 = TransformerConv(hidden_dim, self.head_dim, heads=num_heads,
                                     edge_dim=hidden_dim, concat=True)
        self.bn1 = BatchNorm(hidden_dim)

        self.conv2 = TransformerConv(hidden_dim, self.head_dim, heads=num_heads,
                                     edge_dim=hidden_dim, concat=True)
        self.bn2 = BatchNorm(hidden_dim)

        self.conv3 = TransformerConv(hidden_dim, self.head_dim, heads=num_heads,
                                     edge_dim=hidden_dim, concat=True)
        self.bn3 = BatchNorm(hidden_dim)

        # Virtual node update MLP (one per layer)
        self.vn_mlp_list = nn.ModuleList([
            nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ) for _ in range(3)
        ])

        # ========================================
        # 3. ML global feature encoding layer (SHAP Top-20)
        # ========================================
        self.ml_feature_dim = ml_feature_dim
        if ml_feature_dim > 0:
            self.ml_encoder = nn.Sequential(
                nn.Linear(ml_feature_dim, ml_hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            )
            fusion_dim = hidden_dim + solvent_dim + ml_hidden_dim
        else:
            self.ml_encoder = None
            fusion_dim = hidden_dim + solvent_dim

        # ========================================
        # 4. Prediction head
        # ========================================
        self.mlp = nn.Sequential(
            nn.Linear(fusion_dim, 256),
            nn.ReLU(),
            nn.Dropout(self.dropout_rate),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, 1)
        )

    def forward(self, x, edge_index, edge_attr, solvent_ids, mask_b, batch_index,
                ml_global_features=None, return_attention=False):
        """
        Forward pass.

        Args:
            x: node features [num_nodes, node_in_dim]
            edge_index: edge indices [2, num_edges]
            edge_attr: edge features [num_edges, edge_in_dim]
            solvent_ids: solvent IDs [batch_size]
            mask_b: boron atom mask [num_nodes]
            batch_index: batch assignment for each node [num_nodes]
            ml_global_features: ML global features [batch_size, ml_feature_dim] (optional)
            return_attention: whether to return attention weights (for interpretability)

        Returns:
            out: predicted chemical shifts [num_b_atoms]
            (optional) attention_weights: attention weights
        """
        # ========================================
        # A. Initial encoding
        # ========================================
        x = self.node_encoder(x)
        edge_attr = self.edge_encoder(edge_attr)

        # Initialize virtual nodes
        num_graphs = batch_index.max().item() + 1
        virtual_node_feat = self.virtual_node_encoder(
            torch.zeros(num_graphs, dtype=torch.long, device=x.device)
        )  # [num_graphs, hidden_dim]

        # ========================================
        # B. GNN message passing + virtual node interaction
        # ========================================
        attention_weights = [] if return_attention else None

        # === Layer 1 ===
        x_in = x
        x = x + virtual_node_feat[batch_index]  # virtual node -> real nodes

        if return_attention:
            x, (edge_idx, attn) = self.conv1(x, edge_index, edge_attr,
                                            return_attention_weights=True)
            attention_weights.append((edge_idx, attn))
        else:
            x = self.conv1(x, edge_index, edge_attr)

        x = self.bn1(x)
        x = F.relu(x)
        x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x = x + x_in  # residual connection

        virtual_node_feat = self.vn_mlp_list[0](
            global_add_pool(x, batch_index) + virtual_node_feat
        )  # real nodes -> virtual node

        # === Layer 2 ===
        x_in = x
        x = x + virtual_node_feat[batch_index]

        if return_attention:
            x, (edge_idx, attn) = self.conv2(x, edge_index, edge_attr,
                                            return_attention_weights=True)
            attention_weights.append((edge_idx, attn))
        else:
            x = self.conv2(x, edge_index, edge_attr)

        x = self.bn2(x)
        x = F.relu(x)
        x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x = x + x_in

        virtual_node_feat = self.vn_mlp_list[1](
            global_add_pool(x, batch_index) + virtual_node_feat
        )

        # === Layer 3 ===
        x_in = x
        x = x + virtual_node_feat[batch_index]

        if return_attention:
            x, (edge_idx, attn) = self.conv3(x, edge_index, edge_attr,
                                            return_attention_weights=True)
            attention_weights.append((edge_idx, attn))
        else:
            x = self.conv3(x, edge_index, edge_attr)

        x = self.bn3(x)
        x = F.relu(x)
        x = x + x_in

        virtual_node_feat = self.vn_mlp_list[2](
            global_add_pool(x, batch_index) + virtual_node_feat
        )

        # ========================================
        # C. Extract boron atom features + fuse solvent + ML global features
        # ========================================
        b_features = x[mask_b]
        b_batch_idx = batch_index[mask_b]

        # Get solvent features
        solvent_features = self.solvent_embedding(solvent_ids)  # [batch_size, solvent_dim]

        # Ensure solvent_features is 2D (handles batch_size=1 case)
        if solvent_features.dim() == 1:
            solvent_features = solvent_features.unsqueeze(0)  # [1, solvent_dim]

        b_solvent_features = solvent_features[b_batch_idx]      # [num_b_atoms, solvent_dim]

        # Fuse ML global features (Strategy A: concatenation at prediction head)
        if self.ml_encoder is not None and ml_global_features is not None:
            ml_encoded = self.ml_encoder(ml_global_features)    # [batch_size, ml_hidden_dim]
            if ml_encoded.dim() == 1:
                ml_encoded = ml_encoded.unsqueeze(0)
            b_ml_features = ml_encoded[b_batch_idx]             # [num_b_atoms, ml_hidden_dim]
            combined = torch.cat([b_features, b_solvent_features, b_ml_features], dim=1)
        else:
            combined = torch.cat([b_features, b_solvent_features], dim=1)

        # ========================================
        # D. MLP prediction
        # ========================================
        out = self.mlp(combined).squeeze(-1)

        if return_attention:
            return out, attention_weights
        else:
            return out
