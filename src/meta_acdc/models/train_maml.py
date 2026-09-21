"""MAML++ training on epitope few-shot tasks (KN-6).

MAML++ stabilizations over vanilla MAML:
- per-parameter learnable inner-loop learning rates (meta.MAMLPP)
- multi-step loss: query loss after EACH inner step, weighted
- 3 chained inner adaptation steps (graph-connected via functional_call)

Task regime identical to train_protonet.py (epitope-split, held-out
epitopes, 1:4 hard-negative support sets) for direct comparison.

Usage:
    .venv/bin/python src/meta_acdc/models/train_maml.py \
        --data data/processed/vdjdb.clean.tsv --episodes 6000
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

from meta_acdc.models.meta import MAMLPP
from meta_acdc.models.meta_tasks import TaskSampler, load_vdjdb, split_epitopes
from meta_acdc.models.train_protonet import CDR3Encoder, encode_cdr3, task_tensors

INNER_STEPS = 3
MULTI_STEP_WEIGHTS = [0.3, 0.3, 0.4]


def proto_losses(model: MAMLPP, params: dict, sup, qry, sup_labels, q_labels):
    """Return (query loss, support loss) under the given parameters."""
    emb_s = model.call_encoder(params, sup)
    emb_q = model.call_encoder(params, qry)
    pos_mask = sup_labels == 1
    proto_pos = emb_s[pos_mask].mean(0)
    proto_neg = emb_s[~pos_mask].mean(0)

    d_pos = (emb_q - proto_pos).pow(2).sum(-1)
    d_neg = (emb_q - proto_neg).pow(2).sum(-1)
    q_loss = F.cross_entropy(torch.stack([-d_neg, -d_pos], dim=1), q_labels)

    d_pos_s = (emb_s - proto_pos).pow(2).sum(-1)
    d_neg_s = (emb_s - proto_neg).pow(2).sum(-1)
    s_loss = F.cross_entropy(torch.stack([-d_neg_s, -d_pos_s], dim=1), sup_labels)
    return q_loss, s_loss


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("data/processed/vdjdb.clean.tsv"))
    ap.add_argument("--episodes", type=int, default=6000)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    epitope_tcrs, tcr_epitopes = load_vdjdb(args.data)
    sampler = TaskSampler(epitope_tcrs, tcr_epitopes, seed=42)
    train_eps, val_eps = split_epitopes(sampler, seed=7)
    print(f"train epitopes {len(train_eps)}, val {len(val_eps)}", flush=True)

    model = MAMLPP(CDR3Encoder(hidden=256)).to(args.device)
    outer_opt = torch.optim.Adam(model.parameters(), lr=3e-4, weight_decay=1e-4)

    def sample_task(eps):
        e = sampler.rng.choice(sorted(eps))
        return sampler.sample_task(epitope=e)

    print("training MAML++ ...", flush=True)
    for ep in range(args.episodes):
        task = sample_task(train_eps)
        sup, qry, sup_labels = task_tensors(task, args.device)
        n_q = qry.shape[0] // 2
        q_labels = torch.tensor([1] * n_q + [0] * n_q, dtype=torch.long,
                                device=args.device)
        outer_opt.zero_grad()
        params = dict(model.encoder.named_parameters())
        step_losses = []
        for _ in range(INNER_STEPS):
            q_loss, s_loss = proto_losses(model, params, sup, qry, sup_labels,
                                          q_labels)
            step_losses.append(q_loss)
            params = model.adapt(s_loss, params)
        total = sum(w * l for w, l in zip(MULTI_STEP_WEIGHTS, step_losses))
        total.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        outer_opt.step()
        if ep % 1000 == 0 or ep == args.episodes - 1:
            lo, hi = model.lr_stats()
            print(f"episode {ep}: loss={total.item():.4f} meta-lr=[{lo:.2e},{hi:.2e}]",
                  flush=True)

    # few-shot eval on held-out epitopes
    print("evaluating few-shot ...", flush=True)
    aurocs, auprcs = [], []
    for e in sorted(val_eps):
        task = sampler.sample_task(epitope=e)
        sup, qry, sup_labels = task_tensors(task, args.device)
        n_q = qry.shape[0] // 2
        q_labels = torch.tensor([1] * n_q + [0] * n_q, dtype=torch.long,
                                device=args.device)
        params = dict(model.encoder.named_parameters())
        for _ in range(INNER_STEPS):
            q_loss, s_loss = proto_losses(model, params, sup, qry, sup_labels,
                                          q_labels)
            params = model.adapt(s_loss, params)
        with torch.no_grad():
            emb_q = model.call_encoder(params, qry)
            emb_s = model.call_encoder(params, sup)
            pos_mask = sup_labels == 1
            proto_pos = emb_s[pos_mask].mean(0)
            proto_neg = emb_s[~pos_mask].mean(0)
            d_pos = (emb_q - proto_pos).pow(2).sum(-1)
            d_neg = (emb_q - proto_neg).pow(2).sum(-1)
            scores = (d_neg - d_pos).cpu().numpy()
        y = np.array(task.query_labels)
        aurocs.append(roc_auc_score(y, scores))
        auprcs.append(average_precision_score(y, scores))

    print(f"few-shot over {len(aurocs)} held-out epitopes:")
    print(f"  AUROC mean={np.mean(aurocs):.3f} (+/- {np.std(aurocs):.3f})")
    print(f"  AUPRC mean={np.mean(auprcs):.3f} (+/- {np.std(auprcs):.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
