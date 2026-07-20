"""Task-adapted MACE ablation without the 20 prior ML features."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import torch
from e3nn import o3
from torch import nn
from torch_geometric.nn import global_mean_pool


SCRIPT_DIR = Path(__file__).resolve().parent
BASELINE_MODEL_DIR = SCRIPT_DIR.parent.parent / "gnn_comparison" / "models"
_SPEC = importlib.util.spec_from_file_location(
    "_boron_with_ml_model_mace", BASELINE_MODEL_DIR / "model_mace.py"
)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError("Cannot load the shared MACE layers")
_CORE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_CORE)
GaussianRBF = _CORE.GaussianRBF
MACEInteractionBlock = _CORE.MACEInteractionBlock


class BoronMACEWithoutPriorML(nn.Module):
    """MACE node model with the prior-ML fusion branch removed entirely."""

    def __init__(
        self,
        node_in_dim: int = 58,
        edge_in_dim: int = 11,
        hidden_channels: int = 32,
        num_interactions: int = 2,
        correlation: int = 3,
        lmax: int = 2,
        num_rbf: int = 32,
        radial_hidden_dim: int = 128,
        cutoff: float = 5.0,
        dropout: float = 0.05,
        node_scalar_dim: int = 128,
        num_solvents: int = 11,
        solvent_dim: int = 32,
    ):
        super().__init__()
        if lmax < 1:
            raise ValueError("lmax must be at least 1")
        self.cutoff = float(cutoff)
        self.hidden_channels = int(hidden_channels)
        self.irreps_hidden = o3.Irreps(
            " + ".join(
                f"{hidden_channels}x{ell}{'e' if ell % 2 == 0 else 'o'}"
                for ell in range(lmax + 1)
            )
        )
        self.irreps_sh = o3.Irreps.spherical_harmonics(lmax)
        self.node_encoder = nn.Sequential(
            nn.Linear(node_in_dim, hidden_channels),
            nn.SiLU(),
            nn.LayerNorm(hidden_channels),
        )
        self.rbf = GaussianRBF(num_rbf=num_rbf, cutoff=cutoff)
        self.interactions = nn.ModuleList(
            [
                MACEInteractionBlock(
                    irreps_hidden=self.irreps_hidden,
                    irreps_sh=self.irreps_sh,
                    edge_dim=edge_in_dim,
                    num_rbf=num_rbf,
                    radial_hidden_dim=radial_hidden_dim,
                    correlation=correlation,
                    dropout=dropout,
                )
                for _ in range(num_interactions)
            ]
        )
        self.atomic_readouts = nn.ModuleList(
            [
                o3.Linear(self.irreps_hidden, f"{node_scalar_dim}x0e")
                for _ in range(num_interactions)
            ]
        )
        self.atomic_norm = nn.LayerNorm(node_scalar_dim)
        self.solvent_embedding = nn.Embedding(num_solvents, solvent_dim)
        fusion_dim = 2 * node_scalar_dim + solvent_dim
        self.prediction_head = nn.Sequential(
            nn.Linear(fusion_dim, 256),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(256, 128),
            nn.SiLU(),
            nn.Linear(128, 1),
        )

    def _lift_scalar_nodes(self, scalar_features: torch.Tensor) -> torch.Tensor:
        non_scalar_dim = self.irreps_hidden.dim - self.hidden_channels
        zeros = scalar_features.new_zeros((scalar_features.shape[0], non_scalar_dim))
        return torch.cat([scalar_features, zeros], dim=-1)

    def forward(self, data, return_node_features: bool = False):
        relative = data.pos[data.edge_index[0]] - data.pos[data.edge_index[1]]
        distances = torch.linalg.vector_norm(relative, dim=-1).clamp_min(1.0e-8)
        spherical_harmonics = o3.spherical_harmonics(
            self.irreps_sh,
            relative,
            normalize=True,
            normalization="component",
        )
        radial_edge_features = torch.cat(
            [self.rbf(distances), data.edge_attr], dim=-1
        )

        node_features = self._lift_scalar_nodes(self.node_encoder(data.x))
        scalar_readouts = []
        for interaction, readout in zip(self.interactions, self.atomic_readouts):
            node_features = interaction(
                node_features,
                data.edge_index,
                spherical_harmonics,
                radial_edge_features,
            )
            scalar_readouts.append(readout(node_features))

        scalar_nodes = self.atomic_norm(torch.stack(scalar_readouts).sum(dim=0))
        if return_node_features:
            return node_features, scalar_nodes

        graph_features = global_mean_pool(scalar_nodes, data.batch)
        boron_features = scalar_nodes[data.mask_b]
        boron_batch = data.batch[data.mask_b]
        solvent_features = self.solvent_embedding(data.solvent_id.view(-1))
        combined = torch.cat(
            [
                boron_features,
                graph_features[boron_batch],
                solvent_features[boron_batch],
            ],
            dim=-1,
        )
        return self.prediction_head(combined).squeeze(-1)
