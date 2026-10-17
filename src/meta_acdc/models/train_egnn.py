"""Train a structural native-versus-synthetic-decoy scorer.

Requires an independently accepted immutable release and a new checkpoint path.
Example: python -m meta_acdc.models.train_egnn --release RELEASE_MANIFEST
         --out NEW_RUN/model.pt --epochs 30
Calibration/test folds remain unavailable to model fitting.
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import average_precision_score, roc_auc_score

from meta_acdc.models.egnn import EGNN, PEPTIDE_FLAG_INDEX

NODE_DIM = 26  # 3 phys + pLDDT + is-peptide + 20-dim identity + sidechain extent
EDGE_DIM = 16  # RBF distance (12) + atom-contact histogram (4)
NODE_DIM_V11 = 28  # + esp, burial (atom-level electrostatics)
EDGE_DIM_V11 = 18  # + coul, coul_mag


def collate(graphs: list[tuple], node_dim: int = NODE_DIM,
            edge_dim: int = EDGE_DIM) -> dict:
    """Collate a list of (node_feats, coords, edge_index, edge_feats) into one batch.

    node_dim/edge_dim slice the stored features, so one cached dataset can
    serve feature ablations (e.g. v11 28/18 vs v9.1 26/16).
    """
    from meta_acdc.structure.feature_schema import validate_node_features
    hs, xs, eis, eas, batch, masks = [], [], [], [], [], []
    offset = 0
    for nf, coords, ei, ef in graphs:
        validate_node_features(nf)
        n = len(nf)
        masks.append(torch.tensor([row[PEPTIDE_FLAG_INDEX] > .5 for row in nf], dtype=torch.bool))
        h = torch.tensor([row[:node_dim] for row in nf],
                         dtype=torch.float32)
        h = normalize_node_features(h)
        hs.append(h)
        xs.append(torch.tensor(coords, dtype=torch.float32))
        eis.append(torch.tensor(ei, dtype=torch.long).t() + offset)
        eas.append(torch.tensor([row[:edge_dim] for row in ef],
                                dtype=torch.float32))
        batch.append(torch.full((n,), len(batch), dtype=torch.long))
        offset += n
    return {
        "h": torch.cat(hs),
        "peptide_mask": torch.cat(masks),
        "x": torch.cat(xs),
        "edge_index": torch.cat(eis, dim=1),
        "edge_attr": torch.cat(eas),
        "batch": torch.cat(batch),
    }


def normalize_node_features(h: torch.Tensor) -> torch.Tensor:
    """Per-graph z-score of node features, protecting the peptide flag.

    Batch-level z-scoring made the same graph's features (and therefore its
    prediction) depend on which other graphs share the batch, and z-scored
    the binary is-peptide flag into non-binary values so the `> 0.5` peptide
    mask silently stopped matching peptide nodes. Normalize PER GRAPH and
    restore protected columns verbatim — the same transform is applied by
    score_predictions.score_graph, so train/val/scoring stay consistent.
    (AUDIT.md #3.)
    """
    h = h.clone()
    protected = [c for c in (PEPTIDE_FLAG_INDEX,) if c < h.shape[1]]
    saved = [h[:, c].clone() for c in protected]
    m = h.mean(0, keepdim=True)
    s = h.std(0, keepdim=True)
    h = (h - m) / (s + 1e-6)
    for c, v in zip(protected, saved):
        h[:, c] = v
    return h


def train_val_split(n: int, val_frac: float = 0.2, seed: int = 42) -> tuple[list, list]:
    rng = random.Random(seed)
    idx = list(range(n))
    rng.shuffle(idx)
    n_val = int(n * val_frac)
    return idx[n_val:], idx[:n_val]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=None)
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--node-dim", type=int, default=NODE_DIM,
                    help="slice node features to this width (26 baseline / 28 v11)")
    ap.add_argument("--edge-dim", type=int, default=EDGE_DIM,
                    help="slice edge features to this width (16 baseline / 18 v11)")
    ap.add_argument("--out", type=Path, required=True,
                    help="model output path (default: <data>.model.pt)")
    ap.add_argument("--seed", type=int, default=0,
                    help="torch init + batch-order seed (split seed stays 42)")
    ap.add_argument("--release", type=Path, required=True, help="validated frozen structural release manifest")
    ap.add_argument('--fixed-coordinates', action='store_true',
                    help='invariant distance/message baseline; do not update coordinates')
    args = ap.parse_args()
    from meta_acdc.data.release import load_release, frozen_indices, sha256
    manifest, approved_data = load_release(args.release, data_path=args.data)
    if args.out.exists():raise ValueError("refusing to overwrite checkpoint")
    if args.out.resolve().is_relative_to(args.release.resolve().parent):raise ValueError("checkpoint must be outside immutable release")
    args.data = approved_data
    node_dim, edge_dim = args.node_dim, args.edge_dim

    data = torch.load(args.data, weights_only=False)
    if not (25 <= node_dim <= len(data["graphs"][0][0][0]) and 1 <= edge_dim <= len(data["graphs"][0][3][0])):raise ValueError("feature width incompatible with release")
    graphs, labels = data["graphs"], torch.tensor(data["labels"], dtype=torch.float32)
    pdbs = data["pdbs"]
    print(f"dataset: {len(graphs)} graphs ({int(labels.sum())} pos / {len(labels) - int(labels.sum())} neg)",
          flush=True)
    print(f"config: node_dim={node_dim} edge_dim={edge_dim} epochs={args.epochs} "
          f"seed={args.seed} device={args.device}", flush=True)

    frozen = frozen_indices(manifest, data)
    tr_idx, va_idx = frozen['train'], frozen['validation']
    print(f"frozen split: train={len(tr_idx)} validation={len(va_idx)}; calibration/test untouched")

    va_batch = collate([graphs[i] for i in va_idx], node_dim, edge_dim)
    va_batch = {k: v.to(args.device) for k, v in va_batch.items()}
    y_va = labels[va_idx].to(args.device)

    torch.manual_seed(args.seed)
    net = EGNN(node_dim=node_dim, edge_dim=edge_dim, depth=6, hidden=128, update_coordinates=not args.fixed_coordinates).to(args.device)
    opt = torch.optim.Adam(net.parameters(), lr=3e-4, weight_decay=1e-4)
    loss_fn = nn.BCEWithLogitsLoss()

    batch_size = 32
    rng = random.Random(args.seed)
    for epoch in range(args.epochs):
        net.train()
        order = tr_idx[:]
        rng.shuffle(order)
        total_loss = 0.0
        for start in range(0, len(order), batch_size):
            idx = order[start : start + batch_size]
            batch = collate([graphs[i] for i in idx], node_dim, edge_dim)
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
        tr_batch = collate([graphs[i] for i in tr_idx], node_dim, edge_dim)
        tr_batch = {k: v.to(args.device) for k, v in tr_batch.items()}
        prob_tr = torch.sigmoid(net(**tr_batch).squeeze(-1)).cpu().numpy()
    print(f"train(subset) AUROC={roc_auc_score(labels[tr_idx].numpy(), prob_tr):.3f}",
          flush=True)

    # per-decoy-type val breakdown (displaced = geometry, graft = hard)
    for dtype in ("disp", "<-"):
        sel = [i for i in va_idx if dtype in pdbs[i]] + [i for i in va_idx if labels[i] == 1]
        if not sel:
            continue
        if len({float(labels[i]) for i in sel}) < 2:
            print("per-control AUROC unavailable: both classes required", flush=True)
            continue
        b = collate([graphs[i] for i in sel], node_dim, edge_dim)
        b = {k: v.to(args.device) for k, v in b.items()}
        with torch.no_grad():
            p = torch.sigmoid(net(**b).squeeze(-1)).cpu().numpy()
        y = labels[sel].numpy()
        name = "displaced" if dtype == "disp" else "graft"
        print(f"val [{name}] AUROC={roc_auc_score(y, p):.3f} "
              f"AUPRC={average_precision_score(y, p):.3f} (n={len(sel)})", flush=True)

    out = args.out or Path(args.data).with_suffix(".model.pt")
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"release_sha256": sha256(args.release), "task": manifest["task"], "mask_plddt": False, "update_coordinates": not args.fixed_coordinates, "model": net.state_dict(),
                "n_node_dim": node_dim, "n_edge_dim": edge_dim}, out)
    print(f"saved model to {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
