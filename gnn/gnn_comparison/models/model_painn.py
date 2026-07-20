"""PaiNN-style equivariant network for boron node-level NMR prediction."""

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


def vector_linear(layer: nn.Linear, vectors: torch.Tensor) -> torch.Tensor:
    """Mix vector channels without mixing Cartesian xyz components."""
    return layer(vectors.transpose(1, 2)).transpose(1, 2)


class PaiNNInteraction(nn.Module):
    def __init__(
        self,
        hidden_dim: int,
        edge_dim: int,
        num_rbf: int,
        cutoff: float,
    ):
        super().__init__()
        self.rbf = GaussianRBF(num_rbf=num_rbf, cutoff=cutoff)
        self.scalar_context = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, 3 * hidden_dim),
        )
        self.filter_net = nn.Sequential(
            nn.Linear(num_rbf + edge_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, 3 * hidden_dim),
        )
        self.scalar_norm = nn.LayerNorm(hidden_dim)

    def forward(
        self,
        scalars: torch.Tensor,
        vectors: torch.Tensor,
        pos: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        src, dst = edge_index
        relative = pos[src] - pos[dst]
        distances = torch.linalg.vector_norm(relative, dim=-1).clamp_min(1.0e-8)
        directions = relative / distances[:, None]

        context = self.scalar_context(scalars[src])
        filters = self.filter_net(
            torch.cat([self.rbf(distances), edge_attr], dim=-1)
        )
        message_scalar, gate_vector, gate_direction = (context * filters).chunk(3, dim=-1)
        message_vector = (
            gate_vector[:, :, None] * vectors[src]
            + gate_direction[:, :, None] * directions[:, None, :]
        )

        aggregated_scalar = scatter(
            message_scalar,
            dst,
            dim=0,
            dim_size=scalars.shape[0],
            reduce="mean",
        )
        aggregated_vector = scatter(
            message_vector,
            dst,
            dim=0,
            dim_size=scalars.shape[0],
            reduce="mean",
        )
        scalars = self.scalar_norm(scalars + aggregated_scalar)
        vectors = vectors + aggregated_vector
        return scalars, vectors


class PaiNNUpdate(nn.Module):
    def __init__(self, hidden_dim: int, dropout: float):
        super().__init__()
        self.vector_left = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.vector_right = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.update_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 3 * hidden_dim),
        )
        self.scalar_norm = nn.LayerNorm(hidden_dim)

    def forward(
        self, scalars: torch.Tensor, vectors: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        vector_left = vector_linear(self.vector_left, vectors)
        vector_right = vector_linear(self.vector_right, vectors)
        vector_norm = torch.sqrt(
            torch.sum(vector_right * vector_right, dim=-1) + 1.0e-8
        )
        delta_scalar, scalar_vector_gate, vector_gate = self.update_mlp(
            torch.cat([scalars, vector_norm], dim=-1)
        ).chunk(3, dim=-1)
        invariant_dot = torch.sum(vector_left * vector_right, dim=-1)
        scalars = self.scalar_norm(
            scalars + delta_scalar + scalar_vector_gate * invariant_dot
        )
        vectors = vectors + vector_gate[:, :, None] * vector_left
        return scalars, vectors


class PaiNNBlock(nn.Module):
    def __init__(
        self,
        hidden_dim: int,
        edge_dim: int,
        num_rbf: int,
        cutoff: float,
        dropout: float,
    ):
        super().__init__()
        self.interaction = PaiNNInteraction(
            hidden_dim=hidden_dim,
            edge_dim=edge_dim,
            num_rbf=num_rbf,
            cutoff=cutoff,
        )
        self.update = PaiNNUpdate(hidden_dim=hidden_dim, dropout=dropout)

    def forward(self, scalars, vectors, pos, edge_index, edge_attr):
        scalars, vectors = self.interaction(
            scalars, vectors, pos, edge_index, edge_attr
        )
        return self.update(scalars, vectors)


class BoronPaiNN(nn.Module):
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
        self.blocks = nn.ModuleList(
            [
                PaiNNBlock(
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
        # B scalar + B vector norm + graph scalar + solvent + ML features.
        fusion_dim = 3 * hidden_dim + solvent_dim + ml_hidden_dim
        self.prediction_head = nn.Sequential(
            nn.Linear(fusion_dim, 256),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(256, 128),
            nn.SiLU(),
            nn.Linear(128, 1),
        )

    def forward(self, data, return_node_features: bool = False):
        scalars = self.node_encoder(data.x)
        vectors = torch.zeros(
            (data.x.shape[0], scalars.shape[1], 3),
            dtype=scalars.dtype,
            device=scalars.device,
        )
        for block in self.blocks:
            scalars, vectors = block(
                scalars, vectors, data.pos, data.edge_index, data.edge_attr
            )

        if return_node_features:
            return scalars, vectors

        vector_norms = torch.sqrt(torch.sum(vectors * vectors, dim=-1) + 1.0e-8)
        graph_features = global_mean_pool(scalars, data.batch)
        b_batch = data.batch[data.mask_b]
        solvent_features = self.solvent_embedding(data.solvent_id.view(-1))
        ml_features = self.ml_encoder(data.ml_global_features)
        combined = torch.cat(
            [
                scalars[data.mask_b],
                vector_norms[data.mask_b],
                graph_features[b_batch],
                solvent_features[b_batch],
                ml_features[b_batch],
            ],
            dim=-1,
        )
        return self.prediction_head(combined).squeeze(-1)
