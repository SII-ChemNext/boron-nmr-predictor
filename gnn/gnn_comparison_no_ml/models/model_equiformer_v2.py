"""Compact EquiformerV2 ablation without the 20 prior ML features."""

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
    "_boron_with_ml_model_equiformer_v2",
    BASELINE_MODEL_DIR / "model_equiformer_v2.py",
)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError("Cannot load the shared EquiformerV2 layers")
_CORE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_CORE)
EquiformerV2Block = _CORE.EquiformerV2Block
EquivariantRMSNorm = _CORE.EquivariantRMSNorm
GaussianRBF = _CORE.GaussianRBF


class BoronEquiformerV2WithoutPriorML(nn.Module):
    """EquiformerV2 with the prior-ML fusion branch removed entirely."""

    def __init__(
        self,
        node_in_dim: int = 58,
        edge_in_dim: int = 11,
        sphere_channels: int = 32,
        num_layers: int = 4,
        num_heads: int = 4,
        lmax: int = 2,
        num_rbf: int = 32,
        edge_hidden_dim: int = 128,
        ffn_channels: int = 64,
        grid_resolution: int = 32,
        cutoff: float = 5.0,
        attention_dropout: float = 0.05,
        projection_dropout: float = 0.05,
        use_attention_renorm: bool = True,
        node_scalar_dim: int = 128,
        num_solvents: int = 11,
        solvent_dim: int = 32,
    ):
        super().__init__()
        if lmax < 1:
            raise ValueError("lmax must be at least 1")
        self.cutoff = float(cutoff)
        self.sphere_channels = int(sphere_channels)
        self.irreps_hidden = o3.Irreps(
            " + ".join(
                f"{sphere_channels}x{ell}{'e' if ell % 2 == 0 else 'o'}"
                for ell in range(lmax + 1)
            )
        )
        self.irreps_sh = o3.Irreps.spherical_harmonics(lmax)
        self.node_encoder = nn.Sequential(
            nn.Linear(node_in_dim, sphere_channels),
            nn.SiLU(),
            nn.LayerNorm(sphere_channels),
        )
        self.rbf = GaussianRBF(num_rbf=num_rbf, cutoff=cutoff)
        self.blocks = nn.ModuleList(
            [
                EquiformerV2Block(
                    irreps_hidden=self.irreps_hidden,
                    irreps_sh=self.irreps_sh,
                    lmax=lmax,
                    edge_dim=edge_in_dim,
                    num_rbf=num_rbf,
                    edge_hidden_dim=edge_hidden_dim,
                    num_heads=num_heads,
                    ffn_channels=ffn_channels,
                    sphere_channels=sphere_channels,
                    attention_dropout=attention_dropout,
                    projection_dropout=projection_dropout,
                    use_attention_renorm=use_attention_renorm,
                    grid_resolution=grid_resolution,
                )
                for _ in range(num_layers)
            ]
        )
        self.final_norm = EquivariantRMSNorm(self.irreps_hidden)
        self.scalar_readout = o3.Linear(
            self.irreps_hidden, f"{node_scalar_dim}x0e"
        )
        self.scalar_norm = nn.LayerNorm(node_scalar_dim)
        self.solvent_embedding = nn.Embedding(num_solvents, solvent_dim)
        fusion_dim = 2 * node_scalar_dim + solvent_dim
        self.prediction_head = nn.Sequential(
            nn.Linear(fusion_dim, 256),
            nn.SiLU(),
            nn.Dropout(projection_dropout),
            nn.Linear(256, 128),
            nn.SiLU(),
            nn.Linear(128, 1),
        )

    def _lift_scalar_nodes(self, scalar_features: torch.Tensor) -> torch.Tensor:
        non_scalar_dim = self.irreps_hidden.dim - self.sphere_channels
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
        for block in self.blocks:
            node_features = block(
                node_features,
                data.edge_index,
                spherical_harmonics,
                radial_edge_features,
            )
        node_features = self.final_norm(node_features)
        scalar_nodes = self.scalar_norm(self.scalar_readout(node_features))
        if return_node_features:
            return node_features, scalar_nodes

        graph_features = global_mean_pool(scalar_nodes, data.batch)
        boron_batch = data.batch[data.mask_b]
        solvent_features = self.solvent_embedding(data.solvent_id.view(-1))
        combined = torch.cat(
            [
                scalar_nodes[data.mask_b],
                graph_features[boron_batch],
                solvent_features[boron_batch],
            ],
            dim=-1,
        )
        return self.prediction_head(combined).squeeze(-1)
