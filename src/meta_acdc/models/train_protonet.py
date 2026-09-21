"""Prototypical-network few-shot training on epitope tasks (KN-6).

Encoder: physicochemical profile of CDR3 (fixed-length, padded to 20).
Task: "does TCR bind epitope E?" — binary, 1:4 hard-negative support set.
Evaluation: few-shot AUROC/AUPRC per held-out epitope, mean over epitopes.
Baseline floor: the logistic sequence baseline reached AUROC 0.609 on the
epitope-split; a useful ProtoNet must beat it in the few-shot regime.

Usage:
    .venv/bin/python src/meta_acdc/models/train_protonet.py \
        --data data/processed/vdjdb.clean.tsv --episodes 3000
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import average_precision_score, roc_auc_score

from meta_acdc.models.meta_tasks import TaskSampler, load_vdjdb, split_epitopes

CDR3_MAXLEN = 20
FEAT_PER_AA = 3  # charge, hydrophobicity, volume

AA_PHYS = {
    "A": [0.0, 1.8, 88.6], "R": [1.0, -4.5, 173.4], "N": [0.0, -3.5, 114.1],
    "D": [-1.0, -3.5, 111.1], "C": [0.0, 2.5, 108.5], "Q": [0.0, -3.5, 143.8],
    "E": [-1.0, -3.5, 138.4], "G": [0.0, -0.4, 60.1], "H": [0.1, -3.2, 153.2],
    "I": [0.0, 4.5, 166.7], "L": [0.0, 3.8, 166.7], "K": [1.0, -3.9, 168.6],
    "M": [0.0, 1.9, 162.9], "F": [0.0, 2.8, 189.9], "P": [0.0, -1.6, 112.7],
    "S": [0.0, -0.8, 89.0], "T": [0.0, -0.7, 116.1], "W": [0.0, -0.9, 227.8],
    "Y": [0.0, -1.3, 193.6], "V": [0.0, 4.2, 140.0], "X": [0.0, 0.0, 150.0],
}


def encode_cdr3(cdr3: str) -> np.ndarray:
    vec = np.zeros((CDR3_MAXLEN, FEAT_PER_AA), dtype=np.float32)
    for i, aa in enumerate(cdr3[:CDR3_MAXLEN]):
        vec[i] = AA_PHYS.get(aa, AA_PHYS["X"])
    # z-score per position across the vector (standardize)
    v = vec.ravel()
    return (v - v.mean()) / (v.std() + 1e-6)


class CDR3Encoder(nn.Module):
    def __init__(self, in_dim: int = CDR3_MAXLEN * FEAT_PER_AA, hidden: int = 256):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(in_dim, hidden), nn.SiLU(),
            nn.Linear(hidden, hidden), nn.SiLU(),
            nn.Linear(hidden, hidden),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.mlp(x)


def task_tensors(task, device: str) -> tuple[torch.Tensor, torch.Tensor]:
    sup = torch.tensor(np.vstack([encode_cdr3(t) for t in task.support_pos + task.support_neg]),
                       dtype=torch.float32, device=device)
    qry = torch.tensor(np.vstack([encode_cdr3(t) for t in task.query]),
                       dtype=torch.float32, device=device)
    sup_labels = torch.tensor([1] * len(task.support_pos) + [0] * len(task.support_neg),
                              dtype=torch.long, device=device)
    return sup, qry, sup_labels


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("data/processed/vdjdb.clean.tsv"))
    ap.add_argument("--episodes", type=int, default=10000)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    epitope_tcrs, tcr_epitopes = load_vdjdb(args.data)
    sampler = TaskSampler(epitope_tcrs, tcr_epitopes, seed=42)
    print(f"epitopes with >=10 TCRs: {len(sampler.epitopes)}", flush=True)
    train_eps, val_eps = split_epitopes(sampler, seed=7)
    print(f"train epitopes: {len(train_eps)}, val epitopes: {len(val_eps)}", flush=True)

    encoder = CDR3Encoder().to(args.device)
    opt = torch.optim.Adam(encoder.parameters(), lr=3e-4, weight_decay=1e-4)

    def sample_task(eps: set):
        e = sampler.rng.choice(sorted(eps))
        return sampler.sample_task(epitope=e)

    def proto_loss(sup, qry, sup_labels, device):
        emb_s = encoder(sup)
        emb_q = encoder(qry)
        pos_mask = sup_labels == 1
        proto_pos = emb_s[pos_mask].mean(0)
        proto_neg = emb_s[~pos_mask].mean(0)
        d_pos = (emb_q - proto_pos).pow(2).sum(-1)
        d_neg = (emb_q - proto_neg).pow(2).sum(-1)
        # class 0 = non-binder (-d_neg), class 1 = binder (-d_pos)
        logits = torch.stack([-d_neg, -d_pos], dim=1)  # (Q, 2)
        # query labels: [1]*n_pos + [0]*n_neg
        n_q = qry.shape[0] // 2
        ql = torch.tensor([1] * n_q + [0] * n_q, dtype=torch.long, device=device)
        return F.cross_entropy(logits, ql), (d_neg - d_pos)  # score: >0 = binder

    print("training ...", flush=True)
    for ep in range(args.episodes):
        task = sample_task(train_eps)
        sup, qry, sup_labels = task_tensors(task, args.device)
        opt.zero_grad()
        loss, _ = proto_loss(sup, qry, sup_labels, args.device)
        loss.backward()
        opt.step()
        if ep % 500 == 0 or ep == args.episodes - 1:
            print(f"episode {ep}: loss={loss.item():.4f}", flush=True)

    # few-shot evaluation on held-out epitopes
    print("evaluating few-shot on held-out epitopes ...", flush=True)
    aurocs, auprcs = [], []
    for e in sorted(val_eps):
        task = sampler.sample_task(epitope=e)
        sup, qry, sup_labels = task_tensors(task, args.device)
        with torch.no_grad():
            _, scores = proto_loss(sup, qry, sup_labels, args.device)
        y = np.array(task.query_labels)
        s = scores.cpu().numpy()
        aurocs.append(roc_auc_score(y, s))
        auprcs.append(average_precision_score(y, s))

    print(f"few-shot over {len(aurocs)} held-out epitopes:")
    print(f"  AUROC mean={np.mean(aurocs):.3f} (+/- {np.std(aurocs):.3f})")
    print(f"  AUPRC mean={np.mean(auprcs):.3f} (+/- {np.std(auprcs):.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
