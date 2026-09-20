"""E(n)-Equivariant Graph Neural Network (Satorras et al., ICML 2021).

Implementation notes (from lit review 02, arXiv:2102.09844):
- EGCL update:  m_ij = phi_e(h_i, h_j, ||x_i - x_j||^2, a_ij)
                 x_i' = x_i + C * sum_j (x_i - x_j) * phi_x(m_ij)
                 h_i' = phi_h(h_i, sum_j m_ij)
- Coordinates are tensors here (not scalars); keep them centered (zero-mean)
  and use the normalized displacement (x_i - x_j) / (||x_i - x_j|| + 1) for
  numerical stability (lucidrains/egnn-pytorch practice).
- For our scoring task we need E(3)-INVARIANT outputs: pool node features
  (and optionally coordinate norms) into a graph-level scalar.

Usage:
    net = EGNN(node_dim=..., edge_dim=..., depth=4, hidden=64)
    logit = net(h, x, edge_index, edge_attr, batch)   # graph-level scalar
"""

from __future__ import annotations

import torch
import torch.nn as nn
from torch import Tensor


class EGCL(nn.Module):
    """One equivariant graph convolution layer."""

    def __init__(self, node_dim: int, edge_dim: int, hidden: int = 64):
        super().__init__()
        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * node_dim + 1 + edge_dim, hidden),
            nn.SiLU(),
            nn.Linear(hidden, hidden),
            nn.SiLU(),
            nn.Linear(hidden, node_dim + 1),  # last unit -> coordinate update scalar
        )
        self.node_mlp = nn.Sequential(
            nn.Linear(node_dim + hidden, hidden),
            nn.SiLU(),
            nn.Linear(hidden, node_dim),
        )

    def forward(
        self,
        h: Tensor,            # (N, node_dim)
        x: Tensor,            # (N, 3) coordinates
        edge_index: Tensor,   # (2, E)
        edge_attr: Tensor,    # (E, edge_dim)
    ) -> tuple[Tensor, Tensor]:
        row, col = edge_index
        rel_x = x[row] - x[col]                      # (E, 3)
        dist2 = (rel_x ** 2).sum(dim=-1, keepdim=True)  # (E, 1)
        msg = self.edge_mlp(torch.cat([h[row], h[col], dist2, edge_attr], dim=-1))
        msg_feat, coord_scale = msg[:, :-1], msg[:, -1:]  # (E,node_dim), (E,1)

        # aggregate to nodes
        agg = torch.zeros(h.shape[0], msg_feat.shape[1], device=h.device)
        agg.index_add_(0, row, msg_feat)

        # coordinate update (equivariant): normalized displacement
        eps = rel_x / (rel_x.norm(dim=-1, keepdim=True) + 1.0)
        coord_update = (eps * coord_scale).float()
        agg_x = torch.zeros(x.shape[0], 3, device=x.device, dtype=coord_update.dtype)
        agg_x.index_add_(0, row, coord_update)
        x_new = x + agg_x / max(1, h.shape[0] - 1)   # C = 1/(M-1)

        h_new = self.node_mlp(torch.cat([h, agg], dim=-1))
        return h_new, x_new


class EGNN(nn.Module):
    """Graph-level invariant scorer over an interface graph."""

    def __init__(
        self,
        node_dim: int,
        edge_dim: int,
        depth: int = 4,
        hidden: int = 64,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.embed = nn.Linear(node_dim, hidden)
        self.layers = nn.ModuleList(
            [EGCL(hidden, edge_dim, hidden) for _ in range(depth)]
        )
        self.dropout = nn.Dropout(dropout)
        self.readout = nn.Sequential(
            nn.Linear(2 * hidden + 1, hidden),  # mean_h | max_h | coord-norm
            nn.SiLU(),
            nn.Linear(hidden, 1),
        )

    def forward(
        self,
        h: Tensor,
        x: Tensor,
        edge_index: Tensor,
        edge_attr: Tensor,
        batch: Tensor,
    ) -> Tensor:
        """Return graph-level logits, shape (num_graphs, 1)."""
        # center each graph's coordinates at its centroid (translation invariance)
        centroid = torch.zeros(batch.max().item() + 1, 3, device=x.device, dtype=x.dtype)
        centroid.index_add_(0, batch, x)
        counts = torch.bincount(batch).clamp(min=1).float().unsqueeze(-1)
        centroid = centroid / counts
        x = x - centroid[batch]

        h = self.embed(h)
        for layer in self.layers:
            h, x = layer(h, x, edge_index, edge_attr)
            h = self.dropout(h)

        # invariant pooling: mean + max node features + mean coordinate norm
        n_graphs = batch.max().item() + 1
        node_feat = torch.zeros(n_graphs, h.shape[1], device=h.device)
        node_feat.index_add_(0, batch, h)
        coord_feat = torch.zeros(n_graphs, 1, device=x.device)
        coord_feat.index_add_(0, batch, x.norm(dim=-1, keepdim=True))
        # max pooling via scatter (fill with -inf so padding never wins)
        max_feat = torch.full((n_graphs, h.shape[1]), -1e9, device=h.device)
        max_feat.scatter_reduce_(0, batch.unsqueeze(1).expand_as(h), h,
                                 reduce="amax", include_self=True)
        pooled = torch.cat(
            [node_feat / counts, max_feat, coord_feat / counts], dim=-1
        )
        return self.readout(pooled)
