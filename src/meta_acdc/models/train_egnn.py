"""EGNN prototype training on interface graphs with graft-decoy negatives.

This is the KN-4 pipeline-validation run: can the EGNN distinguish native
TCR-pMHC interfaces from grafted decoys (foreign peptide in the groove)?

Usage:
    .venv/bin/python src/meta_acdc/models/train_egnn.py \
        --data data/processed/egnn_dataset.pt --epochs 30
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import average_precision_score, roc_auc_score

from meta_acdc.models.egnn import EGNN

NODE_DIM = 4   # charge, hydrophobicity, side-chain volume, pLDDT
EDGE_DIM = 12  # RBF distance encoding


def collate(graphs: list[tuple]) -> dict:
    """Collate a list of (node_feats, coords, edge_index, edge_feats) into one batch."""
    hs, xs, eis, eas, batch = [], [], [], [], []
    offset = 0
    for nf, coords, ei, ef in graphs:
        n = len(nf)
        hs.append(torch.tensor(nf, dtype=torch.float32))
        xs.append(torch.tensor(coords, dtype=torch.float32))
        eis.append(torch.tensor(ei, dtype=torch.long).t() + offset)
        eas.append(torch.tensor(ef, dtype=torch.float32))
        batch.append(torch.full((n,), len(batch), dtype=torch.long))
        offset += n
    return {
        "h": torch.cat(hs),
        "x": torch.cat(xs),
        "edge_index": torch.cat(eis, dim=1),
        "edge_attr": torch.cat(eas),
        "batch": torch.cat(batch),
    }


def train_val_split(n: int, val_frac: float = 0.2, seed: int = 42) -> tuple[list, list]:
    rng = random.Random(seed)
    idx = list(range(n))
    rng.shuffle(idx)
    n_val = int(n * val_frac)
    return idx[n_val:], idx[:n_val]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("data/processed/egnn_dataset.pt"))
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    data = torch.load(args.data, weights_only=False)
    graphs, labels = data["graphs"], torch.tensor(data["labels"], dtype=torch.float32)
    print(f"dataset: {len(graphs)} graphs ({int(labels.sum())} pos / {len(labels) - int(labels.sum())} neg)")

    tr_idx, va_idx = train_val_split(len(graphs))
    tr_batch = collate([graphs[i] for i in tr_idx])
    va_batch = collate([graphs[i] for i in va_idx])
    tr_batch = {k: v.to(args.device) for k, v in tr_batch.items()}
    va_batch = {k: v.to(args.device) for k, v in va_batch.items()}
    y_tr = labels[tr_idx].to(args.device)
    y_va = labels[va_idx].to(args.device)

    torch.manual_seed(0)
    net = EGNN(node_dim=NODE_DIM, edge_dim=EDGE_DIM, depth=4, hidden=64).to(args.device)
    opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-4)
    loss_fn = nn.BCEWithLogitsLoss()
    # class imbalance 1:0.94 — mild; use plain BCE for the prototype

    for epoch in range(args.epochs):
        net.train()
        opt.zero_grad()
        logits = net(**tr_batch).squeeze()
        loss = loss_fn(logits, y_tr)
        loss.backward()
        opt.step()
        if epoch % 5 == 0 or epoch == args.epochs - 1:
            net.eval()
            with torch.no_grad():
                prob = torch.sigmoid(net(**va_batch).squeeze()).cpu().numpy()
                y = y_va.cpu().numpy()
            print(f"epoch {epoch:3d}  loss={loss.item():.4f}  "
                  f"val AUROC={roc_auc_score(y, prob):.3f} AUPRC={average_precision_score(y, prob):.3f}")

    net.eval()
    with torch.no_grad():
        prob_tr = torch.sigmoid(net(**tr_batch).squeeze()).cpu().numpy()
    print(f"train AUROC={roc_auc_score(y_tr.cpu().numpy(), prob_tr):.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
