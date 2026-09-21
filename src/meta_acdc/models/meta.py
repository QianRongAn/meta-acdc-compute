"""Meta-learning framework: Prototypical Networks + MAML++ (proposal Module 1).

Lit-review-driven decisions (02-ml-tcr-prediction.md section 5):
- Task definition: "does this TCR bind epitope E?" is a Task (PanPep-style).
- Avoid vanilla MAML (task conflict / over-memorization, RR-ADS 2025);
  use MAML++ stabilizations: per-layer learnable inner LRs + multi-step loss.
- Support sets follow the proposal's 1:4 hard-negative ratio.

MAML++ implementation uses torch.func.functional_call so adapted parameters
remain graph-connected (outer-loop gradients flow through inner steps).
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
from torch.func import functional_call


class PrototypicalNet(nn.Module):
    """Prototypical network head over an encoder."""

    def __init__(self, encoder: nn.Module, hidden: int):
        super().__init__()
        self.encoder = encoder
        self.proj = nn.Linear(hidden, hidden)  # learnable metric projection

    def encode(self, batch: dict) -> Tensor:
        return self.proj(self.encoder(**batch))

    def prototypes(self, support: Tensor, labels: Tensor, n_ways: int) -> Tensor:
        proto = torch.zeros(n_ways, support.shape[1], device=support.device)
        counts = torch.zeros(n_ways, device=support.device)
        proto.index_add_(0, labels, support)
        counts.index_add_(0, labels, torch.ones_like(labels, dtype=torch.float32))
        return proto / counts.clamp(min=1).unsqueeze(-1)

    def forward(self, support, support_labels, query, n_ways) -> Tensor:
        s_emb = self.encode(support)
        q_emb = self.encode(query)
        proto = self.prototypes(s_emb, support_labels, n_ways)
        dist = (q_emb.unsqueeze(1) - proto.unsqueeze(0)).pow(2).sum(dim=-1)
        return -dist


class MAMLPP(nn.Module):
    """MAML++ wrapper (Antoniou et al. 2019) around an encoder.

    Holds per-parameter learnable inner-loop learning rates and exposes
    `adapt()` which returns graph-connected adapted parameters for
    functional_call. The training loop (train_maml.py) chains inner steps
    and applies multi-step loss weights.
    """

    def __init__(self, encoder: nn.Module, init_lr: float = 1e-3):
        super().__init__()
        self.encoder = encoder
        self.param_names = list(dict(encoder.named_parameters()).keys())
        self.params = list(encoder.parameters())
        self.meta_lrs = nn.ParameterList(
            [nn.Parameter(torch.full_like(p, init_lr)) for p in self.params]
        )

    def call_encoder(self, params: dict[str, Tensor], x: Tensor) -> Tensor:
        """Forward the encoder with the given (adapted) parameters."""
        return functional_call(self.encoder, params, (x,))

    def adapt(self, loss: Tensor,
              params: dict[str, Tensor] | None = None) -> dict[str, Tensor]:
        """One inner-loop step: returns adapted params as graph expressions.

        With params=None, adapt from the base (meta-)parameters; pass the
        previous step's dict to chain inner steps.
        """
        if params is None:
            cur = {n: p for n, p in zip(self.param_names, self.params)}
        else:
            cur = params
        ps = [cur[n] for n in self.param_names]
        grads = torch.autograd.grad(loss, ps, create_graph=True,
                                    allow_unused=True)
        adapted: dict[str, Tensor] = {}
        for n, g, lr in zip(self.param_names, grads, self.meta_lrs):
            adapted[n] = cur[n] - lr * (g if g is not None else 0.0)
        return adapted

    def lr_stats(self) -> tuple[float, float]:
        vals = torch.cat([lr.data.flatten() for lr in self.meta_lrs])
        return float(vals.min()), float(vals.max())


def hard_negative_ratio(support_labels: Tensor, positive_id: int = 0) -> float:
    """Report neg/pos ratio of a support set (proposal target: 1:4)."""
    n_pos = (support_labels == positive_id).sum().item()
    n_neg = (support_labels != positive_id).sum().item()
    return n_neg / max(n_pos, 1)


__all__ = ["PrototypicalNet", "MAMLPP", "hard_negative_ratio"]
