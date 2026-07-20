"""Compact EquiformerV2 architecture for boron node-level NMR prediction.

This task-adapted implementation retains the defining symmetry mechanisms used
by EquiformerV2: higher-degree SO(3)/O(3) node representations, spherical
harmonic edge bases, distance-conditioned equivariant tensor products,
invariant multi-head attention with attention re-normalization, pre-norm
residual blocks, and a separable gated equivariant feed-forward network.

Coordinates are read-only.  The B-node output is built from even l=0 channels
and is therefore invariant to rigid translations, rotations and reflections.
"""

from __future__ import annotations

import math

import torch
from torch import nn
import torch.nn.functional as F
from torch_geometric.nn import global_mean_pool
from torch_geometric.utils import scatter, softmax

try:
    from e3nn import o3
    from e3nn.nn import Dropout, S2Activation
except ImportError as error:  # pragma: no cover
    raise ImportError(
        "04-EquiformerV2 requires e3nn==0.5.8. Install with: "
        "python -m pip install --index-url https://pypi.org/simple e3nn==0.5.8"
    ) from error


class GaussianRBF(nn.Module):
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


class EquivariantRMSNorm(nn.Module):
    """Per-node RMS normalization with one gain per irrep copy.

    The RMS is invariant because e3nn representations are orthogonal.  Gains
    are shared over all m components of each irrep copy, preserving exact
    equivariance unlike an ordinary feature-wise LayerNorm.
    """

    def __init__(self, irreps: o3.Irreps, eps: float = 1.0e-8):
        super().__init__()
        self.irreps = o3.Irreps(irreps)
        gain_index = []
        copy_index = 0
        for multiplicity, irrep in self.irreps:
            for _ in range(multiplicity):
                gain_index.extend([copy_index] * irrep.dim)
                copy_index += 1
        self.weight = nn.Parameter(torch.ones(copy_index))
        self.register_buffer("gain_index", torch.tensor(gain_index, dtype=torch.long))
        self.eps = float(eps)

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        rms = torch.sqrt(torch.mean(features * features, dim=-1, keepdim=True) + self.eps)
        gains = self.weight[self.gain_index].view(1, -1)
        return features * gains / rms


class RadialTensorProduct(nn.Module):
    def __init__(
        self,
        irreps_node: o3.Irreps,
        irreps_sh: o3.Irreps,
        edge_dim: int,
        num_rbf: int,
        edge_hidden_dim: int,
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
        self.radial = nn.Sequential(
            nn.Linear(num_rbf + edge_dim, edge_hidden_dim),
            nn.SiLU(),
            nn.Linear(edge_hidden_dim, edge_hidden_dim),
            nn.SiLU(),
            nn.Linear(edge_hidden_dim, self.tensor_product.weight_numel),
        )

    def forward(
        self,
        source_features: torch.Tensor,
        spherical_harmonics: torch.Tensor,
        radial_edge_features: torch.Tensor,
    ) -> torch.Tensor:
        return self.tensor_product(
            source_features,
            spherical_harmonics,
            self.radial(radial_edge_features),
        )


class SO3MultiheadAttention(nn.Module):
    """Invariant attention weights applied to equivariant value messages."""

    def __init__(
        self,
        irreps_hidden: o3.Irreps,
        irreps_sh: o3.Irreps,
        edge_dim: int,
        num_rbf: int,
        edge_hidden_dim: int,
        num_heads: int,
        attention_dropout: float,
        use_attention_renorm: bool,
    ):
        super().__init__()
        if num_heads < 1:
            raise ValueError("num_heads must be positive")
        self.num_heads = int(num_heads)
        self.edge_tensor_product = RadialTensorProduct(
            irreps_node=irreps_hidden,
            irreps_sh=irreps_sh,
            edge_dim=edge_dim,
            num_rbf=num_rbf,
            edge_hidden_dim=edge_hidden_dim,
        )
        self.query = o3.Linear(irreps_hidden, irreps_hidden)
        self.key = o3.Linear(irreps_hidden, irreps_hidden)
        # Contracting query and key to 0e produces rotation-invariant logits.
        self.logit_product = o3.FullyConnectedTensorProduct(
            irreps_hidden,
            irreps_hidden,
            f"{num_heads}x0e",
            irrep_normalization="component",
            path_normalization="element",
        )
        self.radial_bias = nn.Sequential(
            nn.Linear(num_rbf + edge_dim, edge_hidden_dim),
            nn.SiLU(),
            nn.Linear(edge_hidden_dim, num_heads),
        )
        self.attention_renorm = (
            nn.LayerNorm(num_heads) if use_attention_renorm and num_heads > 1 else nn.Identity()
        )
        self.attention_dropout = nn.Dropout(attention_dropout)
        self.value_heads = nn.ModuleList(
            [o3.Linear(irreps_hidden, irreps_hidden) for _ in range(num_heads)]
        )
        self.output = o3.Linear(irreps_hidden, irreps_hidden)
        self.logit_scale = 1.0 / math.sqrt(float(irreps_hidden.dim))

    def forward(
        self,
        node_features: torch.Tensor,
        edge_index: torch.Tensor,
        spherical_harmonics: torch.Tensor,
        radial_edge_features: torch.Tensor,
    ) -> torch.Tensor:
        src, dst = edge_index
        edge_features = self.edge_tensor_product(
            node_features[src], spherical_harmonics, radial_edge_features
        )
        query = self.query(node_features)[dst]
        key = self.key(edge_features)
        logits = self.logit_product(query, key) * self.logit_scale
        logits = self.attention_renorm(logits + self.radial_bias(radial_edge_features))
        weights = softmax(logits, dst, num_nodes=node_features.shape[0])
        weights = self.attention_dropout(weights)

        aggregated = node_features.new_zeros(node_features.shape)
        for head, value_projection in enumerate(self.value_heads):
            values = value_projection(edge_features)
            head_messages = weights[:, head : head + 1] * values
            aggregated = aggregated + scatter(
                head_messages,
                dst,
                dim=0,
                dim_size=node_features.shape[0],
                reduce="sum",
            )
        return self.output(aggregated / math.sqrt(float(self.num_heads)))


class ChannelwiseS2Activation(nn.Module):
    """Apply an S² grid activation independently to every spherical channel."""

    def __init__(self, channels: int, lmax: int, resolution: int):
        super().__init__()
        self.channels = int(channels)
        self.lmax = int(lmax)
        self.single_channel_irreps = o3.Irreps.spherical_harmonics(lmax)
        self.activation = S2Activation(
            self.single_channel_irreps,
            F.silu,
            res=resolution,
            normalization="component",
        )

    def _to_channel_major(self, features: torch.Tensor) -> torch.Tensor:
        pieces = []
        start = 0
        for ell in range(self.lmax + 1):
            irrep_dim = 2 * ell + 1
            length = self.channels * irrep_dim
            pieces.append(
                features[:, start : start + length].reshape(
                    features.shape[0], self.channels, irrep_dim
                )
            )
            start += length
        return torch.cat(pieces, dim=-1).reshape(
            features.shape[0] * self.channels, -1
        )

    def _from_channel_major(
        self, features: torch.Tensor, number_of_nodes: int
    ) -> torch.Tensor:
        features = features.reshape(number_of_nodes, self.channels, -1)
        pieces = []
        start = 0
        for ell in range(self.lmax + 1):
            irrep_dim = 2 * ell + 1
            pieces.append(
                features[:, :, start : start + irrep_dim].reshape(
                    number_of_nodes, self.channels * irrep_dim
                )
            )
            start += irrep_dim
        return torch.cat(pieces, dim=-1)

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        channel_major = self._to_channel_major(features)
        activated = self.activation(channel_major)
        return self._from_channel_major(activated, features.shape[0])


class SeparableS2FeedForward(nn.Module):
    """EquiformerV2 separable scalar/S² feed-forward nonlinearity."""

    def __init__(
        self,
        irreps_hidden: o3.Irreps,
        lmax: int,
        input_channels: int,
        ffn_channels: int,
        dropout: float,
        grid_resolution: int,
    ):
        super().__init__()
        self.ffn_channels = int(ffn_channels)
        self.ffn_irreps = o3.Irreps(
            " + ".join(
                f"{ffn_channels}x{ell}{'e' if ell % 2 == 0 else 'o'}"
                for ell in range(lmax + 1)
            )
        )
        self.input = o3.Linear(irreps_hidden, self.ffn_irreps)
        self.scalar_path = nn.Linear(input_channels, ffn_channels)
        self.s2_activation = ChannelwiseS2Activation(
            channels=ffn_channels,
            lmax=lmax,
            resolution=grid_resolution,
        )
        self.output = o3.Linear(self.ffn_irreps, irreps_hidden)
        self.dropout = Dropout(irreps_hidden, p=dropout)

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        projected = self.s2_activation(self.input(features))
        # EquiformerV2's separable path derives the new l=0 coefficients
        # directly from the input scalar channels rather than from the grid.
        scalar_features = F.silu(self.scalar_path(features[:, : self.scalar_path.in_features]))
        projected = torch.cat(
            [scalar_features, projected[:, self.ffn_channels :]], dim=-1
        )
        return self.dropout(self.output(projected))


class EquiformerV2Block(nn.Module):
    def __init__(
        self,
        irreps_hidden: o3.Irreps,
        irreps_sh: o3.Irreps,
        lmax: int,
        edge_dim: int,
        num_rbf: int,
        edge_hidden_dim: int,
        num_heads: int,
        ffn_channels: int,
        sphere_channels: int,
        attention_dropout: float,
        projection_dropout: float,
        use_attention_renorm: bool,
        grid_resolution: int,
    ):
        super().__init__()
        self.norm_attention = EquivariantRMSNorm(irreps_hidden)
        self.attention = SO3MultiheadAttention(
            irreps_hidden=irreps_hidden,
            irreps_sh=irreps_sh,
            edge_dim=edge_dim,
            num_rbf=num_rbf,
            edge_hidden_dim=edge_hidden_dim,
            num_heads=num_heads,
            attention_dropout=attention_dropout,
            use_attention_renorm=use_attention_renorm,
        )
        self.attention_dropout = Dropout(irreps_hidden, p=projection_dropout)
        self.norm_ffn = EquivariantRMSNorm(irreps_hidden)
        self.feed_forward = SeparableS2FeedForward(
            irreps_hidden=irreps_hidden,
            lmax=lmax,
            input_channels=sphere_channels,
            ffn_channels=ffn_channels,
            dropout=projection_dropout,
            grid_resolution=grid_resolution,
        )

    def forward(
        self,
        node_features: torch.Tensor,
        edge_index: torch.Tensor,
        spherical_harmonics: torch.Tensor,
        radial_edge_features: torch.Tensor,
    ) -> torch.Tensor:
        attention_update = self.attention(
            self.norm_attention(node_features),
            edge_index,
            spherical_harmonics,
            radial_edge_features,
        )
        node_features = node_features + self.attention_dropout(attention_update)
        return node_features + self.feed_forward(self.norm_ffn(node_features))


class BoronEquiformerV2(nn.Module):
    """Higher-degree equivariant Transformer with a B-node scalar readout."""

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
        ml_feature_dim: int = 20,
        ml_hidden_dim: int = 64,
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
        self.ml_encoder = nn.Sequential(
            nn.Linear(ml_feature_dim, ml_hidden_dim),
            nn.SiLU(),
            nn.LayerNorm(ml_hidden_dim),
        )
        fusion_dim = 2 * node_scalar_dim + solvent_dim + ml_hidden_dim
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
        ml_features = self.ml_encoder(data.ml_global_features)
        combined = torch.cat(
            [
                scalar_nodes[data.mask_b],
                graph_features[boron_batch],
                solvent_features[boron_batch],
                ml_features[boron_batch],
            ],
            dim=-1,
        )
        return self.prediction_head(combined).squeeze(-1)
