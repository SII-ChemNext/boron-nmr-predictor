"""PaiNN ablation without the 20 prior ML features."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import torch
from torch import nn
from torch_geometric.nn import global_mean_pool


SCRIPT_DIR = Path(__file__).resolve().parent
BASELINE_MODEL_DIR = SCRIPT_DIR.parent.parent / "gnn_comparison" / "models"
_SPEC = importlib.util.spec_from_file_location(
    "_boron_with_ml_model_painn", BASELINE_MODEL_DIR / "model_painn.py"
)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError("Cannot load the shared PaiNN layers")
_CORE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_CORE)
PaiNNBlock = _CORE.PaiNNBlock


class BoronPaiNNWithoutPriorML(nn.Module):
    """PaiNN with the prior-ML encoder and fusion branch removed entirely."""

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
        fusion_dim = 3 * hidden_dim + solvent_dim
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
        boron_batch = data.batch[data.mask_b]
        solvent_features = self.solvent_embedding(data.solvent_id.view(-1))
        combined = torch.cat(
            [
                scalars[data.mask_b],
                vector_norms[data.mask_b],
                graph_features[boron_batch],
                solvent_features[boron_batch],
            ],
            dim=-1,
        )
        return self.prediction_head(combined).squeeze(-1)
