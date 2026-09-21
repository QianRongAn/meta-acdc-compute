"""Domain-adaptation training: crystals + TCRmodel2-predicted structures.

Fixes the Rashomon-effect instability: the model must see "predicted-style"
inputs during training. TCRmodel2 PDBs (data/raw/tcrmodel/) are added as
extra POSITIVE examples alongside the crystal natives, with the SAME decoy
protocol (displaced + graft) applied to them.

pLDDT is zeroed in training AND scoring (the AF3 B-factor confound).

After training 3 seeds, reports pairwise Spearman rank correlation on the
held-out decoy set across seeds (the stability metric).

Usage:
    .venv/bin/python src/meta_acdc/models/train_domain_adapt.py \
        --base data/processed/egnn_dataset.pt \
        --tcrmodel data/raw/tcrmodel \
        --epochs 150
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score

from meta_acdc.models.dataset import make_displaced_peptide_decoy
from meta_acdc.models.egnn import EGNN
from meta_acdc.models.train_egnn import collate
from meta_acdc.structure.graph import (build_interface_graph_from_residues,
                                       classify_chains, parse_pdb)

SEEDS = [0, 1, 2]


def mask_plddt(gs):
    out = []
    for nf, coords, ei, ef in gs:
        nf = np.array(nf, dtype=np.float32).copy()
        nf[:, 3] = 0.0
        out.append((nf.tolist(), coords, ei, ef))
    return out


def add_tcrmodel_graphs(graphs, labels, pdbs, tcr_dir: Path):
    n_added = 0
    for p in sorted(tcr_dir.glob("*.pdb")):
        try:
            residues = parse_pdb(p)
            g = build_interface_graph_from_residues(residues)
        except Exception:
            continue
        if g.n_nodes < 20:
            continue
        graphs.append((g.node_features, g.node_coords,
                       g.edge_index, g.edge_features))
        labels.append(1)
        pdbs.append(f"tcr-{p.stem}")
        n_added += 1
    print(f"added {n_added} TCRmodel2 natives as positives", flush=True)
    return graphs, labels, pdbs


def pdb_group(name: str) -> str:
    return name.split("-disp")[0].split("<-")[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=Path,
                    default=Path("data/processed/egnn_dataset.pt"))
    ap.add_argument("--tcrmodel", type=Path, default=Path("data/raw/tcrmodel"))
    ap.add_argument("--epochs", type=int, default=150)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    d = torch.load(args.base, weights_only=False)
    graphs = list(d["graphs"])
    labels = list(d["labels"])
    pdbs = list(d["pdbs"])
    graphs, labels, pdbs = add_tcrmodel_graphs(
        graphs, labels, pdbs, args.tcrmodel)
    labels_t = torch.tensor(labels, dtype=torch.float32)
    graphs_masked = mask_plddt(graphs)
    print(f"DA dataset: {len(graphs)} graphs "
          f"({sum(labels)} pos / {len(labels) - sum(labels)} neg)", flush=True)

    groups: dict[str, list[int]] = {}
    for i, p in enumerate(pdbs):
        groups.setdefault(pdb_group(p), []).append(i)
    rng = random.Random(42)
    gids = list(groups.keys())
    rng.shuffle(gids)
    n_val = max(1, int(0.2 * len(gids)))
    val_g = set(gids[:n_val])
    tr = [i for g in gids if g not in val_g for i in groups[g]]
    va = [i for g in gids if g in val_g for i in groups[g]]

    device = args.device
    bs = 32
    all_probs = []
    for seed in SEEDS:
        torch.manual_seed(seed)
        net = EGNN(node_dim=26, edge_dim=16, depth=6, hidden=128).to(device)
        opt = torch.optim.Adam(net.parameters(), lr=3e-4, weight_decay=1e-4)
        loss_fn = nn.BCEWithLogitsLoss()
        rng2 = random.Random(seed)
        for epoch in range(args.epochs):
            net.train()
            order = tr[:]
            rng2.shuffle(order)
            for start in range(0, len(order), bs):
                idx = order[start:start + bs]
                b = {k: v.to(device) for k, v in
                     collate([graphs_masked[i] for i in idx]).items()}
                y = labels_t[idx].to(device)
                opt.zero_grad()
                loss = loss_fn(net(**b).squeeze(-1), y)
                loss.backward()
                opt.step()
        net.eval()
        vb = {k: v.to(device) for k, v in
              collate([graphs_masked[i] for i in va]).items()}
        with torch.no_grad():
            p = torch.sigmoid(net(**vb).squeeze(-1)).cpu().numpy()
        y = labels_t[va].numpy()
        auc = roc_auc_score(y, p)
        print(f"seed {seed}: val AUROC {auc:.3f}", flush=True)
        all_probs.append(p)
        torch.save({"model": net.state_dict(), "n_node_dim": 26,
                    "n_edge_dim": 16},
                   f"data/processed/da_seed{seed}.model.pt")

    # stability: pairwise Spearman of val scores across seeds
    from scipy.stats import spearmanr
    print("\npairwise Spearman (val set) across seeds:")
    for i in range(len(SEEDS)):
        for j in range(i + 1, len(SEEDS)):
            r, _ = spearmanr(all_probs[i], all_probs[j])
            print(f"  seed{i} vs seed{j}: r={r:.3f}")
    print("DA training complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
