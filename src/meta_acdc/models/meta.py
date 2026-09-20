"""Meta-learning framework: Prototypical Networks + MAML++ (proposal Module 1).

Lit-review-driven decisions (02-ml-tcr-prediction.md section 5):
- Task definition: "predicting the binding affinity of a specific TCR against
  diverse peptides" is a Task. At v0 we use epitope-cluster tasks (PanPep-style);
  structure-cluster tasks arrive with the AF3/TCRmodel2 structure pipeline.
- Avoid vanilla MAML (task conflict / over-memorization, RR-ADS 2025).
  Implement ProtoNet first (stable, cheap), then MAML++ (multi-step loss,
  per-layer learnable LRs) for the EGNN encoder.
- Support-set construction follows the proposal's 1:4 hard-negative ratio.

This module depends on the EGNN encoder interface:
    encode(batch) -> (num_graphs, hidden) embeddings.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor


class PrototypicalNet(nn.Module):
    """Prototypical network head over an encoder.

    Encode support/query graphs, compute per-class prototypes from support
    embeddings, classify queries by distance to prototypes.
    """

    def __init__(self, encoder: nn.Module, hidden: int):
        super().__init__()
        self.encoder = encoder
        self.proj = nn.Linear(hidden, hidden)  # learnable metric projection

    def encode(self, batch: dict) -> Tensor:
        return self.proj(self.encoder(**batch))

    def prototypes(self, support: Tensor, labels: Tensor, n_ways: int) -> Tensor:
        """support: (n_support, hidden); labels: (n_support,) class ids."""
        proto = torch.zeros(n_ways, support.shape[1], device=support.device)
        counts = torch.zeros(n_ways, device=support.device)
        proto.index_add_(0, labels, support)
        counts.index_add_(0, labels, torch.ones_like(labels, dtype=torch.float32))
        return proto / counts.clamp(min=1).unsqueeze(-1)

    def forward(
        self,
        support: dict,
        support_labels: Tensor,
        query: dict,
        n_ways: int,
    ) -> Tensor:
        """Return query logits (n_query, n_ways)."""
        s_emb = self.encode(support)
        q_emb = self.encode(query)
        proto = self.prototypes(s_emb, support_labels, n_ways)
        # squared Euclidean distance -> logits
        dist = (
            q_emb.unsqueeze(1) - proto.unsqueeze(0)
        ).pow(2).sum(dim=-1)  # (n_query, n_ways)
        return -dist


class MAMLPP(nn.Module):
    """MAML++-style wrapper (Antoniou et al. 2019) around an EGNN encoder.

    Implements the training loop helpers only; the full bi-level optimization
    (multi-step loss, per-layer learnable LRs) is assembled in the training
    script. This class holds the per-parameter meta-LRs.

    NOTE: v0 placeholder — implemented and tuned at benchmark time (KN-4).
    """

    def __init__(self, encoder: nn.Module):
        super().__init__()
        self.encoder = encoder
        self.meta_lrs = nn.ParameterDict({
            name: nn.Parameter(torch.full_like(p, 1e-3))
            for name, p in encoder.named_parameters()
        })

    def inner_step(self, loss: Tensor) -> None:
        """One inner-loop (adaptation) step with per-layer learned LRs."""
        grads = torch.autograd.grad(loss, self.encoder.parameters(),
                                    create_graph=True)
        with torch.no_grad():
            for (name, p), g in zip(self.encoder.named_parameters(), grads):
                p.sub_(self.meta_lrs[name] * g)


def hard_negative_ratio(support_labels: Tensor, positive_id: int = 0) -> float:
    """Report neg/pos ratio of a support set (proposal target: 1:4)."""
    n_pos = (support_labels == positive_id).sum().item()
    n_neg = (support_labels != positive_id).sum().item()
    return n_neg / max(n_pos, 1)


__all__ = ["PrototypicalNet", "MAMLPP", "hard_negative_ratio"]
