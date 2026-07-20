"""Fixed-coordinate EGNN for boron node-level NMR shift prediction."""

from __future__ import annotations

import torch
from torch import nn
from torch_geometric.nn import global_mean_pool
from torch_geometric.utils import scatter


class GaussianRBF(nn.Module):
    def __init__(self, num_rbf: int, cutoff: float):
        super().__init__()
        centers = torch.linspace(0.0, cutoff, num_rbf)
        self.register_buffer("centers", centers)
        spacing = float(centers[1] - centers[0]) if num_rbf > 1 else cutoff
        self.gamma = 1.0 / max(spacing * spacing, 1.0e-6)
        self.cutoff = cutoff

    def forward(self, distances: torch.Tensor) -> torch.Tensor:
        rbf = torch.exp(-self.gamma * (distances[:, None] - self.centers) ** 2)
        envelope = 0.5 * (torch.cos(torch.pi * distances / self.cutoff) + 1.0)
        envelope = envelope * (distances < self.cutoff)
        return rbf * envelope[:, None]


class FixedCoordinateEGNNLayer(nn.Module):
    """EGNN message layer whose coordinate map is the identity.

    Coordinates are read to construct invariant distances but are never updated.
    Consequently the node states and final scalar predictions are invariant to
    translation, rotation and reflection of the input conformer.
    """

    def __init__(
        self,
        hidden_dim: int,
        edge_dim: int,
        num_rbf: int,
        cutoff: float,
        dropout: float,
    ):
        super().__init__()
        self.rbf = GaussianRBF(num_rbf=num_rbf, cutoff=cutoff)
        message_input_dim = 2 * hidden_dim + edge_dim + num_rbf
        self.message_mlp = nn.Sequential(
            nn.Linear(message_input_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
        )
        self.node_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.norm = nn.LayerNorm(hidden_dim)

    def forward(
        self,
        h: torch.Tensor,
        pos: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
    ) -> torch.Tensor:
        src, dst = edge_index
        distances = torch.linalg.vector_norm(pos[src] - pos[dst], dim=-1)
        messages = self.message_mlp(
            torch.cat([h[src], h[dst], self.rbf(distances), edge_attr], dim=-1)
        )
        aggregated = scatter(
            messages, dst, dim=0, dim_size=h.shape[0], reduce="mean"
        )
        update = self.node_mlp(torch.cat([h, aggregated], dim=-1))
        return self.norm(h + update)


class BoronEGNN(nn.Module):
    def __init__(
        self,
        node_in_dim: int = 58,
        edge_in_dim: int = 11,
        hidden_dim: int = 128,
        num_layers: int = 4,
        num_rbf: int = 32,
        cutoff: float = 5.0,
        dropout: float = 0.05,
        num_solvents: int = 11,
        solvent_dim: int = 32,
        ml_feature_dim: int = 20,
        ml_hidden_dim: int = 64,
    ):
        super().__init__()
        self.node_encoder = nn.Sequential(
            nn.Linear(node_in_dim, hidden_dim),
            nn.SiLU(),
            nn.LayerNorm(hidden_dim),
        )
        self.layers = nn.ModuleList(
            [
                FixedCoordinateEGNNLayer(
                    hidden_dim=hidden_dim,
                    edge_dim=edge_in_dim,
                    num_rbf=num_rbf,
                    cutoff=cutoff,
                    dropout=dropout,
                )
                for _ in range(num_layers)
            ]
        )
        self.solvent_embedding = nn.Embedding(num_solvents, solvent_dim)
        self.ml_encoder = nn.Sequential(
            nn.Linear(ml_feature_dim, ml_hidden_dim),
            nn.SiLU(),
            nn.LayerNorm(ml_hidden_dim),
        )
        fusion_dim = 2 * hidden_dim + solvent_dim + ml_hidden_dim
        self.prediction_head = nn.Sequential(
            nn.Linear(fusion_dim, 256),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(256, 128),
            nn.SiLU(),
            nn.Linear(128, 1),
        )

    def forward(self, data, return_node_features: bool = False):
        h = self.node_encoder(data.x)
        # data.pos is deliberately passed unchanged to every layer.
        for layer in self.layers:
            h = layer(h, data.pos, data.edge_index, data.edge_attr)

        if return_node_features:
            return h

        graph_features = global_mean_pool(h, data.batch)
        b_features = h[data.mask_b]
        b_batch = data.batch[data.mask_b]
        solvent_features = self.solvent_embedding(data.solvent_id.view(-1))
        ml_features = self.ml_encoder(data.ml_global_features)
        combined = torch.cat(
            [
                b_features,
                graph_features[b_batch],
                solvent_features[b_batch],
                ml_features[b_batch],
            ],
            dim=-1,
        )
        return self.prediction_head(combined).squeeze(-1)
