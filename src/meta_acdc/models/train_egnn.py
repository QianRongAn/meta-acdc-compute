"""EGNN prototype training on interface graphs with graft-decoy negatives.

This is the KN-4 pipeline-validation run: can the EGNN distinguish native
TCR-pMHC interfaces from grafted decoys (foreign peptide in the groove)?

Usage:
    .venv/bin/python src/meta_acdc/models/train_egnn.py \
        --data data/processed/egnn_dataset.pt --epochs 30

    # v11 feature ablation on the electrostatics-augmented dataset:
    .venv/bin/python src/meta_acdc/models/train_egnn.py \
        --data data/processed/egnn_dataset_v11.pt \
        --node-dim 28 --edge-dim 18 --out .../v11_full.model.pt
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
    hs, xs, eis, eas, batch = [], [], [], [], []
    offset = 0
    for nf, coords, ei, ef in graphs:
        n = len(nf)
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
    ap.add_argument("--data", type=Path, default=Path("data/processed/egnn_dataset.pt"))
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--node-dim", type=int, default=NODE_DIM,
                    help="slice node features to this width (26 baseline / 28 v11)")
    ap.add_argument("--edge-dim", type=int, default=EDGE_DIM,
                    help="slice edge features to this width (16 baseline / 18 v11)")
    ap.add_argument("--out", type=Path, default=None,
                    help="model output path (default: <data>.model.pt)")
    ap.add_argument("--seed", type=int, default=0,
                    help="torch init + batch-order seed (split seed stays 42)")
    ap.add_argument("--group-by", choices=["pdb", "family"], default="family",
                    help="fold grouping: pdb = same complex only; family = "
                         "also merges PDBs sharing a receptor (recommended; "
                         "needs --family-map)")
    ap.add_argument("--family-map", type=Path,
                    default=Path("data/processed/pdb_family_map.tsv"),
                    help="pdb -> family_id map (structure/family_map.py)")
    args = ap.parse_args()
    node_dim, edge_dim = args.node_dim, args.edge_dim

    data = torch.load(args.data, weights_only=False)
    graphs, labels = data["graphs"], torch.tensor(data["labels"], dtype=torch.float32)
    pdbs = data["pdbs"]
    print(f"dataset: {len(graphs)} graphs ({int(labels.sum())} pos / {len(labels) - int(labels.sum())} neg)",
          flush=True)
    print(f"config: node_dim={node_dim} edge_dim={edge_dim} epochs={args.epochs} "
          f"seed={args.seed} device={args.device}", flush=True)

    # Complex-identity grouping: every decoy and every re-prediction of the
    # same PDB must land in ONE fold. The old key only stripped "-disp"/"<-",
    # so "5brz_0", "5brz_1" and "af3-5brz_EVDPIGHLY" were three independent
    # groups and the crystal version of a complex could sit in validation
    # while its predicted copy was in training (AUDIT.md #4).
    family: dict[str, str] = {}
    if args.group_by == "family":
        fpath = args.family_map
        if not fpath.exists():
            print(f"family map {fpath} not found; falling back to pdb grouping",
                  flush=True)
            args.group_by = "pdb"
        else:
            with open(fpath, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    family[r["pdb"]] = r["family_id"]

    def pdb_group(name: str) -> str:
        n = name.split("<-")[0].split("-disp")[0]
        for pre in ("af3-", "tcr-", "inst-", "nativetcr-"):
            if n.startswith(pre):
                n = n[len(pre):]
        head = n.split("_")[0]          # peptide / replication suffix
        # PDB ids are exactly 4 characters; synthetic stems stay unique
        pdb = head if len(head) == 4 else n
        return family.get(pdb, pdb) if family else pdb

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

    va_batch = collate([graphs[i] for i in va_idx], node_dim, edge_dim)
    va_batch = {k: v.to(args.device) for k, v in va_batch.items()}
    y_va = labels[va_idx].to(args.device)

    torch.manual_seed(args.seed)
    net = EGNN(node_dim=node_dim, edge_dim=edge_dim, depth=6, hidden=128).to(args.device)
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
        tr_batch = collate([graphs[i] for i in tr_idx[:100]], node_dim, edge_dim)
        tr_batch = {k: v.to(args.device) for k, v in tr_batch.items()}
        prob_tr = torch.sigmoid(net(**tr_batch).squeeze(-1)).cpu().numpy()
    print(f"train(subset) AUROC={roc_auc_score(labels[tr_idx[:100]].numpy(), prob_tr):.3f}",
          flush=True)

    # per-decoy-type val breakdown (displaced = geometry, graft = hard)
    for dtype in ("disp", "<-"):
        sel = [i for i in va_idx if dtype in pdbs[i]] + [i for i in va_idx if labels[i] == 1]
        if not sel:
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
    torch.save({"model": net.state_dict(),
                "n_node_dim": node_dim, "n_edge_dim": edge_dim}, out)
    print(f"saved model to {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
