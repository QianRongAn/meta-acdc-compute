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

NODE_DIM = 5   # charge, hydrophobicity, side-chain volume, pLDDT, is-peptide
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
    h = torch.cat(hs)
    # batch-level z-score normalization (mixed physical scales)
    h = (h - h.mean(0, keepdim=True)) / (h.std(0, keepdim=True) + 1e-6)
    return {
        "h": h,
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
    pdbs = data["pdbs"]
    print(f"dataset: {len(graphs)} graphs ({int(labels.sum())} pos / {len(labels) - int(labels.sum())} neg)",
          flush=True)

    # PDB-grouped split: all variants (native + decoys) of one PDB share a fold
    def pdb_group(name: str) -> str:
        return name.split("-disp")[0].split("<-")[0]

    groups: dict[str, list[int]] = {}
    for i, p in enumerate(pdbs):
        groups.setdefault(pdb_group(p), []).append(i)
    rng = random.Random(42)
    gids = list(groups.keys())
    rng.shuffle(gids)
    n_val_groups = max(1, int(0.2 * len(gids)))
    val_groups = set(gids[:n_val_groups])
    tr_idx = [i for g in gids if g not in val_groups for i in groups[g]]
    va_idx = [i for g in gids if g in val_groups for i in groups[g]]
    print(f"split: {len(gids)} PDB groups -> train={len(tr_idx)} val={len(va_idx)}", flush=True)

    va_batch = collate([graphs[i] for i in va_idx])
    va_batch = {k: v.to(args.device) for k, v in va_batch.items()}
    y_va = labels[va_idx].to(args.device)

    torch.manual_seed(0)
    net = EGNN(node_dim=NODE_DIM, edge_dim=EDGE_DIM, depth=6, hidden=128).to(args.device)
    opt = torch.optim.Adam(net.parameters(), lr=3e-4, weight_decay=1e-4)
    loss_fn = nn.BCEWithLogitsLoss()

    batch_size = 32
    rng = random.Random(42)
    for epoch in range(args.epochs):
        net.train()
        order = tr_idx[:]
        rng.shuffle(order)
        total_loss = 0.0
        for start in range(0, len(order), batch_size):
            idx = order[start : start + batch_size]
            batch = collate([graphs[i] for i in idx])
            batch = {k: v.to(args.device) for k, v in batch.items()}
            y = labels[idx].to(args.device)
            opt.zero_grad()
            logits = net(**batch).squeeze(-1)
            loss = loss_fn(logits, y)
            loss.backward()
            opt.step()
            total_loss += loss.item()
        if epoch % 5 == 0 or epoch == args.epochs - 1:
            net.eval()
            with torch.no_grad():
                prob = torch.sigmoid(net(**va_batch).squeeze(-1)).cpu().numpy()
                y = y_va.cpu().numpy()
            print(f"epoch {epoch:3d}  loss={total_loss / max(1, len(order) // batch_size):.4f}  "
                  f"val AUROC={roc_auc_score(y, prob):.3f} AUPRC={average_precision_score(y, prob):.3f}",
                  flush=True)

    net.eval()
    with torch.no_grad():
        tr_batch = collate([graphs[i] for i in tr_idx[:100]])
        tr_batch = {k: v.to(args.device) for k, v in tr_batch.items()}
        prob_tr = torch.sigmoid(net(**tr_batch).squeeze(-1)).cpu().numpy()
    print(f"train(subset) AUROC={roc_auc_score(labels[tr_idx[:100]].numpy(), prob_tr):.3f}",
          flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
