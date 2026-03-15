import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATv2Conv, BatchNorm, global_add_pool

class BoronNMRNet_GATv2(nn.Module):
    """
    GATv2-based 11B NMR prediction model (with ML Feature Top-20 concatenation)

    Features:
    1. GATv2Conv (improved graph attention network)
    2. Virtual Node mechanism
    3. Learnable solvent embedding (32-dim)
    4. 4 attention heads
    5. [Added] ML global feature fusion (SHAP Top-20)
    """

    def __init__(self, node_in_dim, edge_in_dim,
                 num_solvents=11, solvent_dim=32,
                 hidden_dim=256, dropout=0.0355, heads=4,
                 ml_feature_dim=0, ml_hidden_dim=64):
        super().__init__()

        self.dropout_rate = dropout
        self.hidden_dim = hidden_dim
        self.heads = heads

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
        # 2. GNN backbone (GATv2Conv + virtual node)
        # ========================================
        self.conv1 = GATv2Conv(hidden_dim, hidden_dim, heads=heads,
                              edge_dim=hidden_dim, concat=False)
        self.bn1 = BatchNorm(hidden_dim)

        self.conv2 = GATv2Conv(hidden_dim, hidden_dim, heads=heads,
                              edge_dim=hidden_dim, concat=False)
        self.bn2 = BatchNorm(hidden_dim)

        self.conv3 = GATv2Conv(hidden_dim, hidden_dim, heads=heads,
                              edge_dim=hidden_dim, concat=False)
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
                ml_global_features=None):
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

        Returns:
            out: predicted chemical shifts [num_b_atoms]
        """
        # ========================================
        # A. Initial encoding
        # ========================================
        x = self.node_encoder(x)
        edge_attr = self.edge_encoder(edge_attr)

        # Initialize virtual nodes
        num_graphs = batch_index.max().item() + 1
        vn = self.virtual_node_encoder(torch.zeros(num_graphs, dtype=torch.long, device=x.device))

        # ========================================
        # B. GNN message passing + virtual node update (3 layers)
        # ========================================
        conv_layers = [(self.conv1, self.bn1), (self.conv2, self.bn2), (self.conv3, self.bn3)]

        for i, (conv, bn) in enumerate(conv_layers):
            # 1) Broadcast virtual node to all nodes
            x = x + vn[batch_index]

            # 2) GATv2 convolution
            x_new = conv(x, edge_index, edge_attr=edge_attr)
            x_new = bn(x_new)
            x_new = F.relu(x_new)
            x_new = F.dropout(x_new, p=self.dropout_rate, training=self.training)

            # 3) Residual connection
            x = x + x_new

            # 4) Update virtual node (global pooling + MLP)
            vn_update = global_add_pool(x, batch_index)
            vn = vn + self.vn_mlp_list[i](vn_update)

        # ========================================
        # C. Boron atom prediction
        # ========================================
        # Extract boron atom features
        b_features = x[mask_b]
        b_batch_idx = batch_index[mask_b]

        # Get solvent features
        solvent_features = self.solvent_embedding(solvent_ids)

        # Ensure solvent_features is 2D (handles batch_size=1 case)
        if solvent_features.dim() == 1:
            solvent_features = solvent_features.unsqueeze(0)

        b_solvent_features = solvent_features[b_batch_idx]

        # Fuse ML global features
        if self.ml_encoder is not None and ml_global_features is not None:
            ml_encoded = self.ml_encoder(ml_global_features)    # [batch_size, ml_hidden_dim]
            if ml_encoded.dim() == 1:
                ml_encoded = ml_encoded.unsqueeze(0)
            b_ml_features = ml_encoded[b_batch_idx]             # [num_b_atoms, ml_hidden_dim]
            b_combined = torch.cat([b_features, b_solvent_features, b_ml_features], dim=1)
        else:
            b_combined = torch.cat([b_features, b_solvent_features], dim=1)

        # Predict
        out = self.mlp(b_combined).squeeze(-1)

        return out
