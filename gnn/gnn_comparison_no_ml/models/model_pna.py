import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import PNAConv, BatchNorm, global_add_pool

class BoronNMRNet_PNA(nn.Module):
    """
    PNA-based 11B NMR prediction model

    Features:
    1. PNAConv (Principal Neighbourhood Aggregation)
    2. Virtual Node mechanism
    3. Learnable solvent embedding (32-dim)
    4. Multiple aggregator combination
    """

    def __init__(self, node_in_dim, edge_in_dim,
                 num_solvents=11, solvent_dim=32,
                 hidden_dim=256, dropout=0.0355,
                 deg=None):  # deg is the degree statistics
        super().__init__()

        self.dropout_rate = dropout
        self.hidden_dim = hidden_dim

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
        # 2. GNN backbone (PNAConv + virtual node)
        # ========================================
        # PNA aggregator configuration
        aggregators = ['mean', 'min', 'max', 'std']
        scalers = ['identity', 'amplification', 'attenuation']

        self.conv1 = PNAConv(hidden_dim, hidden_dim,
                            aggregators=aggregators,
                            scalers=scalers,
                            deg=deg,
                            edge_dim=hidden_dim,
                            towers=1,
                            pre_layers=1,
                            post_layers=1,
                            divide_input=False)
        self.bn1 = BatchNorm(hidden_dim)

        self.conv2 = PNAConv(hidden_dim, hidden_dim,
                            aggregators=aggregators,
                            scalers=scalers,
                            deg=deg,
                            edge_dim=hidden_dim,
                            towers=1,
                            pre_layers=1,
                            post_layers=1,
                            divide_input=False)
        self.bn2 = BatchNorm(hidden_dim)

        self.conv3 = PNAConv(hidden_dim, hidden_dim,
                            aggregators=aggregators,
                            scalers=scalers,
                            deg=deg,
                            edge_dim=hidden_dim,
                            towers=1,
                            pre_layers=1,
                            post_layers=1,
                            divide_input=False)
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
        # 3. Prediction head
        # ========================================
        fusion_dim = hidden_dim + solvent_dim

        self.mlp = nn.Sequential(
            nn.Linear(fusion_dim, 256),
            nn.ReLU(),
            nn.Dropout(self.dropout_rate),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, 1)
        )

    def forward(self, x, edge_index, edge_attr, solvent_ids, mask_b, batch_index):
        """
        Forward pass.

        Args:
            x: node features [num_nodes, node_in_dim]
            edge_index: edge indices [2, num_edges]
            edge_attr: edge features [num_edges, edge_in_dim]
            solvent_ids: solvent IDs [batch_size]
            mask_b: boron atom mask [num_nodes]
            batch_index: batch assignment for each node [num_nodes]

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

            # 2) PNA convolution
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

        # Concatenate features
        b_combined = torch.cat([b_features, b_solvent_features], dim=1)

        # Predict
        out = self.mlp(b_combined).squeeze(-1)

        return out
