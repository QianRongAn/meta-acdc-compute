"""Does domain adaptation preserve the core crystal decoy performance?

Evaluates each DA seed (and the original v9.1 checkpoint) on the held-out
PDB-grouped crystal decoy set from egnn_dataset.pt — the same protocol as
train_egnn.py. The DA models were trained with pLDDT zeroed, so the eval
masks it too. Answers: did mixing AF3 natives into training erode the
in-distribution graft/displaced discrimination?

Usage:
    .venv/bin/python src/meta_acdc/models/eval_da_crystal.py \
        --da data/processed/da_seed0.model.pt \
             data/processed/da_seed1.model.pt \
             data/processed/da_seed2.model.pt \
        --baseline data/processed/egnn_dataset.model.pt
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import torch
from sklearn.metrics import average_precision_score, roc_auc_score

from meta_acdc.models.egnn import EGNN
from meta_acdc.models.train_egnn import collate


def mask_plddt(graphs):
    out = []
    for nf, coords, ei, ef in graphs:
        nf = np.array(nf, dtype=np.float32).copy()
        nf[:, 3] = 0.0
        out.append((nf.tolist(), coords, ei, ef))
    return out


def evaluate(ckpt_path: Path, va_idx, graphs_masked, labels, pdbs, device):
    ck = torch.load(ckpt_path, weights_only=False)
    net = EGNN(node_dim=ck["n_node_dim"], edge_dim=ck["n_edge_dim"],
               depth=6, hidden=128).to(device)
    net.load_state_dict(ck["model"])
    net.eval()
    res = {}
    for dtype, name in (("disp", "displaced"), ("<-", "graft")):
        sel = [i for i in va_idx if dtype in pdbs[i]] + \
              [i for i in va_idx if labels[i] == 1]
        b = {k: v.to(device) for k, v in
             collate([graphs_masked[i] for i in sel]).items()}
        with torch.no_grad():
            p = torch.sigmoid(net(**b).squeeze(-1)).cpu().numpy()
        y = labels[sel]
        res[name] = (roc_auc_score(y, p), average_precision_score(y, p),
                     len(sel))
    # overall val
    b = {k: v.to(device) for k, v in
         collate([graphs_masked[i] for i in va_idx]).items()}
    with torch.no_grad():
        p = torch.sigmoid(net(**b).squeeze(-1)).cpu().numpy()
    y = labels[va_idx]
    res["overall"] = (roc_auc_score(y, p), average_precision_score(y, p),
                      len(va_idx))
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path,
                    default=Path("data/processed/egnn_dataset.pt"))
    ap.add_argument("--da", type=Path, nargs="*", default=[])
    ap.add_argument("--baseline", type=Path,
                    default=Path("data/processed/egnn_dataset.model.pt"))
    ap.add_argument("--device",
                    default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    data = torch.load(args.data, weights_only=False)
    graphs = list(data["graphs"])
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

    gm = mask_plddt(graphs)
    models = [("v9.1 (baseline)", args.baseline)] + \
             [(f"DA seed {i}", p) for i, p in enumerate(args.da)]
    print(f"held-out crystal val: {len(va_idx)} graphs "
          f"({sum(labels[i] for i in va_idx)} pos)\n")
    print(f"{'model':18s} {'overall':>8s} {'displaced':>10s} {'graft':>8s}")
    for name, path in models:
        if not path.exists():
            print(f"{name:18s} (missing {path})")
            continue
        r = evaluate(path, va_idx, gm, labels, pdbs, args.device)
        print(f"{name:18s} {r['overall'][0]:8.3f} {r['displaced'][0]:10.3f} "
              f"{r['graft'][0]:8.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())