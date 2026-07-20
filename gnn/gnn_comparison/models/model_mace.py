"""Task-adapted MACE network for boron node-level NMR shift prediction.

The model keeps coordinates fixed.  Relative vectors are expanded in spherical
harmonics and combined with node irreducible representations using e3nn tensor
products.  Repeated on-site tensor products provide the higher-body-order
product basis characteristic of MACE.  Only even scalar (l=0) channels are
passed to the final boron-node readout, so predicted shifts are O(3)-invariant.
"""

from __future__ import annotations

import math

import torch
from torch import nn
import torch.nn.functional as F
from torch_geometric.nn import global_mean_pool
from torch_geometric.utils import scatter

try:
    from e3nn import o3
    from e3nn.nn import BatchNorm, Dropout, NormActivation
except ImportError as error:  # pragma: no cover - gives users an actionable error
    raise ImportError(
        "03-MACE requires e3nn==0.5.8. Install with: "
        "python -m pip install --index-url https://pypi.org/simple e3nn==0.5.8"
    ) from error


class GaussianRBF(nn.Module):
    """Gaussian radial basis with a smooth cosine cutoff envelope."""

    def __init__(self, num_rbf: int, cutoff: float):
        super().__init__()
        centers = torch.linspace(0.0, cutoff, num_rbf)
        self.register_buffer("centers", centers)
        spacing = float(centers[1] - centers[0]) if num_rbf > 1 else cutoff
        self.gamma = 1.0 / max(spacing * spacing, 1.0e-6)
        self.cutoff = float(cutoff)

    def forward(self, distances: torch.Tensor) -> torch.Tensor:
        rbf = torch.exp(-self.gamma * (distances[:, None] - self.centers) ** 2)
        envelope = 0.5 * (
            torch.cos(torch.pi * distances / self.cutoff) + 1.0
        )
        envelope = envelope * (distances < self.cutoff)
        return rbf * envelope[:, None]


class RadialTensorProduct(nn.Module):
    """Distance/edge-conditioned equivariant tensor product for edge messages."""

    def __init__(
        self,
        irreps_node: o3.Irreps,
        irreps_sh: o3.Irreps,
        edge_dim: int,
        num_rbf: int,
        radial_hidden_dim: int,
    ):
        super().__init__()
        self.tensor_product = o3.FullyConnectedTensorProduct(
            irreps_node,
            irreps_sh,
            irreps_node,
            shared_weights=False,
            internal_weights=False,
            irrep_normalization="component",
            path_normalization="element",
        )
        self.radial_mlp = nn.Sequential(
            nn.Linear(num_rbf + edge_dim, radial_hidden_dim),
            nn.SiLU(),
            nn.Linear(radial_hidden_dim, radial_hidden_dim),
            nn.SiLU(),
            nn.Linear(radial_hidden_dim, self.tensor_product.weight_numel),
        )

    def forward(
        self,
        source_features: torch.Tensor,
        spherical_harmonics: torch.Tensor,
        radial_edge_features: torch.Tensor,
    ) -> torch.Tensor:
        weights = self.radial_mlp(radial_edge_features)
        return self.tensor_product(source_features, spherical_harmonics, weights)


class MACEProductBasis(nn.Module):
    """Higher-order on-site product basis up to a requested correlation order.

    MACE's efficient symmetric contraction is represented here by successive
    Clebsch-Gordan tensor products.  With correlation=3 the block contains
    linear, quadratic and cubic equivariant terms of the aggregated neighbor
    density while retaining the same output irreps after every contraction.
    """

    def __init__(self, irreps: o3.Irreps, correlation: int):
        super().__init__()
        if correlation < 1:
            raise ValueError("correlation must be at least 1")
        self.products = nn.ModuleList(
            [
                o3.FullyConnectedTensorProduct(
                    irreps,
                    irreps,
                    irreps,
                    irrep_normalization="component",
                    path_normalization="element",
                )
                for _ in range(correlation - 1)
            ]
        )
        self.mix = o3.Linear(irreps, irreps)
        self.scale = 1.0 / math.sqrt(float(correlation))

    def forward(self, density: torch.Tensor) -> torch.Tensor:
        result = density
        order_term = density
        for tensor_product in self.products:
            order_term = tensor_product(order_term, density)
            result = result + order_term
        return self.mix(result * self.scale)


class MACEInteractionBlock(nn.Module):
    """One equivariant neighbor-density and higher-order product-basis block."""

    def __init__(
        self,
        irreps_hidden: o3.Irreps,
        irreps_sh: o3.Irreps,
        edge_dim: int,
        num_rbf: int,
        radial_hidden_dim: int,
        correlation: int,
        dropout: float,
    ):
        super().__init__()
        self.edge_messages = RadialTensorProduct(
            irreps_node=irreps_hidden,
            irreps_sh=irreps_sh,
            edge_dim=edge_dim,
            num_rbf=num_rbf,
            radial_hidden_dim=radial_hidden_dim,
        )
        self.message_linear = o3.Linear(irreps_hidden, irreps_hidden)
        self.product_basis = MACEProductBasis(
            irreps=irreps_hidden, correlation=correlation
        )
        self.skip = o3.Linear(irreps_hidden, irreps_hidden)
        self.norm = BatchNorm(
            irreps_hidden,
            normalization="component",
            instance=False,
            affine=True,
            reduce="mean",
        )
        self.activation = NormActivation(
            irreps_hidden,
            scalar_nonlinearity=F.silu,
            normalize=True,
            epsilon=1.0e-8,
            bias=False,
        )
        self.dropout = Dropout(irreps_hidden, p=dropout)

    def forward(
        self,
        node_features: torch.Tensor,
        edge_index: torch.Tensor,
        spherical_harmonics: torch.Tensor,
        radial_edge_features: torch.Tensor,
    ) -> torch.Tensor:
        src, dst = edge_index
        messages = self.edge_messages(
            node_features[src], spherical_harmonics, radial_edge_features
        )
        density = scatter(
            messages,
            dst,
            dim=0,
            dim_size=node_features.shape[0],
            reduce="mean",
        )
        density = self.message_linear(density)
        update = self.product_basis(density)
        output = self.norm(self.skip(node_features) + update)
        return self.dropout(self.activation(output))


class BoronMACE(nn.Module):
    """MACE node encoder with an invariant B-atom chemical-shift readout."""

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
        ml_feature_dim: int = 20,
        ml_hidden_dim: int = 64,
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
        # MACE uses an atomic readout after each interaction.  Their sum is an
        # invariant node representation used here instead of an energy sum.
        self.atomic_readouts = nn.ModuleList(
            [
                o3.Linear(self.irreps_hidden, f"{node_scalar_dim}x0e")
                for _ in range(num_interactions)
            ]
        )
        self.atomic_norm = nn.LayerNorm(node_scalar_dim)
        self.solvent_embedding = nn.Embedding(num_solvents, solvent_dim)
        self.ml_encoder = nn.Sequential(
            nn.Linear(ml_feature_dim, ml_hidden_dim),
            nn.SiLU(),
            nn.LayerNorm(ml_hidden_dim),
        )
        fusion_dim = 2 * node_scalar_dim + solvent_dim + ml_hidden_dim
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
        ml_features = self.ml_encoder(data.ml_global_features)
        combined = torch.cat(
            [
                boron_features,
                graph_features[boron_batch],
                solvent_features[boron_batch],
                ml_features[boron_batch],
            ],
            dim=-1,
        )
        return self.prediction_head(combined).squeeze(-1)
