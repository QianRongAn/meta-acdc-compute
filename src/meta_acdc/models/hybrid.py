"""Hybrid sequence+graph model (graft-decoy diagnosis, v7).

Rationale: the pure-EGNN cannot learn chemical compatibility from 149
structure examples (each groove appears once; peptide node features are the
only signal). Sequence information is the missing chemistry channel:
- Branch A: peptide sequence encoder (1D CNN over one-hot, fixed len 15)
- Branch B: EGNN over the interface graph (geometry)
- Late fusion -> classifier.

The peptide sequence is extracted from peptide-node one-hot features ordered
by (resid, icode) — no dataset change needed.
"""

from __future__ import annotations

import torch
import torch.nn as nn
from torch import Tensor

from meta_acdc.models.egnn import EGNN

PEPTIDE_MAXLEN = 15
N_AA = 20

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {aa: i for i, aa in enumerate(AA_ORDER)}


class PeptideSeqEncoder(nn.Module):
    """1D CNN over one-hot peptide sequences (length-padded to 15)."""

    def __init__(self, hidden: int = 128):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv1d(N_AA, hidden, kernel_size=3, padding=1),
            nn.SiLU(),
            nn.Conv1d(hidden, hidden, kernel_size=3, padding=1),
            nn.SiLU(),
        )
        self.pool = nn.AdaptiveMaxPool1d(1)

    def forward(self, seq: Tensor) -> Tensor:  # (B, L, 20)
        return self.pool(self.conv(seq.transpose(1, 2))).squeeze(-1)


class HybridModel(nn.Module):
    def __init__(self, node_dim: int = 25, edge_dim: int = 12,
                 hidden: int = 128, egnn_depth: int = 6):
        super().__init__()
        self.egnn = EGNN(node_dim=node_dim, edge_dim=edge_dim,
                         depth=egnn_depth, hidden=hidden)
        self.seq_enc = PeptideSeqEncoder(hidden=hidden)
        self.head = nn.Sequential(
            nn.Linear(2 * hidden, hidden),
            nn.SiLU(),
            nn.Linear(hidden, 1),
        )

    def forward(self, h: Tensor, x: Tensor, edge_index: Tensor,
                edge_attr: Tensor, batch: Tensor, seq: Tensor) -> Tensor:
        g = self.egnn(h, x, edge_index, edge_attr, batch)  # (B, 1)
        s = self.seq_enc(seq)                               # (B, hidden)
        return self.head(torch.cat([s, g], dim=-1))
