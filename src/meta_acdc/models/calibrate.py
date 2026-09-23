"""Temperature calibration for EGNN compatibility scores.

The raw sigmoid scores are not calibrated binding probabilities (BCE training
only constrains within-group ordering). Temperature scaling fits a single
scalar T on a held-out labeled set (crystal decoy val split) so that
sigmoid(logit / T) is calibrated — a prerequisite for any absolute score
threshold in a safety screen. AUROC is invariant to T; calibration affects
the probability scale and expected calibration error (ECE).

Usage:
    .venv/bin/python src/meta_acdc/models/calibrate.py \
        --model data/processed/da_seed0.model.pt \
        --data data/processed/egnn_dataset.pt \
        --out data/processed/da_seed0.temperature.json
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score

from meta_acdc.models.egnn import EGNN
from meta_acdc.models.train_egnn import collate


def mask_plddt(graphs):
    out = []
    for nf, coords, ei, ef in graphs:
        nf = np.array(nf, dtype=np.float32).copy()
        nf[:, 3] = 0.0
        out.append((nf.tolist(), coords, ei, ef))
    return out


def ece(probs: np.ndarray, labels: np.ndarray, n_bins: int = 10) -> float:
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    total = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (probs > lo) & (probs <= hi)
        if m.sum() == 0:
            continue
        conf = probs[m].mean()
        acc = labels[m].mean()
        total += m.mean() * abs(conf - acc)
    return float(total)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--data", type=Path,
                    default=Path("data/processed/egnn_dataset.pt"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--device",
                    default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    data = torch.load(args.data, weights_only=False)
    graphs = mask_plddt(list(data["graphs"]))
    labels = np.asarray(data["labels"])
    pdbs = list(data["pdbs"])

    def pdb_group(name: str) -> str:
        return name.split("-disp")[0].split("<-")[0]

    groups: dict[str, list[int]] = {}
    for i, p in enumerate(pdbs):
        groups.setdefault(pdb_group(p), []).append(i)
    rng = random.Random(42)
    gids = list(groups.keys())
    rng.shuffle(gids)
    n_val = max(1, int(0.2 * len(gids)))
    val_groups = set(gids[:n_val])
    va_idx = [i for g in gids if g in val_groups for i in groups[g]]

    ck = torch.load(args.model, weights_only=False)
    net = EGNN(node_dim=ck["n_node_dim"], edge_dim=ck["n_edge_dim"],
               depth=6, hidden=128).to(args.device)
    net.load_state_dict(ck["model"])
    net.eval()
    b = {k: v.to(args.device) for k, v in
         collate([graphs[i] for i in va_idx]).items()}
    with torch.no_grad():
        logits = net(**b).squeeze(-1).cpu()
    y = torch.tensor(labels[va_idx], dtype=torch.float32)

    p0 = torch.sigmoid(logits).numpy()
    print(f"n={len(y)} pos={int(y.sum())}")
    print(f"before: AUROC={roc_auc_score(y, p0):.3f} "
          f"ECE={ece(p0, y.numpy()):.3f}")

    # fit temperature by minimizing BCE on logits/T
    logT = torch.zeros(1, requires_grad=True)
    opt = torch.optim.LBFGS([logT], lr=0.1, max_iter=100)
    loss_fn = nn.BCEWithLogitsLoss()

    def closure():
        opt.zero_grad()
        loss = loss_fn(logits / torch.exp(logT), y)
        loss.backward()
        return loss

    opt.step(closure)
    T = float(torch.exp(logT).item())
    p1 = torch.sigmoid(logits / T).detach().numpy()
    print(f"after : T={T:.3f} AUROC={roc_auc_score(y, p1):.3f} "
          f"ECE={ece(p1, y.numpy()):.3f}")

    with open(args.out, "w") as fh:
        json.dump({"model": str(args.model), "temperature": T,
                   "ece_before": ece(p0, y.numpy()),
                   "ece_after": ece(p1, y.numpy())}, fh, indent=2)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())