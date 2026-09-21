"""Train the hybrid sequence+graph model on the decoy task (v7 experiment).

Usage:
    .venv/bin/python src/meta_acdc/models/train_hybrid.py \
        --data data/processed/egnn_dataset.pt --epochs 200
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

from meta_acdc.models.hybrid import HybridModel, AA_IDX, PEPTIDE_MAXLEN, N_AA

NODE_DIM = 25
EDGE_DIM = 16


def collate(graphs: list[tuple]) -> tuple[dict, torch.Tensor]:
    hs, xs, eis, eas, batch = [], [], [], [], []
    seqs = []
    offset = 0
    for nf, coords, ei, ef in graphs:
        nf = torch.tensor(nf, dtype=torch.float32)
        n = len(nf)
        hs.append(nf)
        xs.append(torch.tensor(coords, dtype=torch.float32))
        eis.append(torch.tensor(ei, dtype=torch.long).t() + offset)
        eas.append(torch.tensor(ef, dtype=torch.float32))
        batch.append(torch.full((n,), len(batch), dtype=torch.long))
        offset += n
        # peptide sequence from node one-hot (nodes are in file order)
        pep = nf[nf[:, 4] > 0.5]
        if len(pep):
            seq = pep[:, 5:25].argmax(-1).tolist()
        else:
            seq = []
        onehot = torch.zeros(PEPTIDE_MAXLEN, N_AA)
        for i, aa in enumerate(seq[:PEPTIDE_MAXLEN]):
            onehot[i, aa] = 1.0
        seqs.append(onehot)
    h = torch.cat(hs)
    h = (h - h.mean(0, keepdim=True)) / (h.std(0, keepdim=True) + 1e-6)
    return {
        "h": h,
        "x": torch.cat(xs),
        "edge_index": torch.cat(eis, dim=1),
        "edge_attr": torch.cat(eas),
        "batch": torch.cat(batch),
        "seq": torch.stack(seqs),
    }


def pdb_group(name: str) -> str:
    return name.split("-disp")[0].split("<-")[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("data/processed/egnn_dataset.pt"))
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    data = torch.load(args.data, weights_only=False)
    graphs, labels = data["graphs"], torch.tensor(data["labels"], dtype=torch.float32)
    pdbs = data["pdbs"]
    print(f"dataset: {len(graphs)} graphs", flush=True)

    groups: dict[str, list[int]] = {}
    for i, p in enumerate(pdbs):
        groups.setdefault(pdb_group(p), []).append(i)
    rng = random.Random(42)
    gids = list(groups.keys())
    rng.shuffle(gids)
    n_val = max(1, int(0.2 * len(gids)))
    val_g = set(gids[:n_val])
    tr_idx = [i for g in gids if g not in val_g for i in groups[g]]
    va_idx = [i for g in gids if g in val_g for i in groups[g]]
    print(f"split: train={len(tr_idx)} val={len(va_idx)}", flush=True)

    def to_dev(b: dict) -> dict:
        return {k: v.to(args.device) for k, v in b.items()}

    va_batch = to_dev(collate([graphs[i] for i in va_idx]))
    y_va = labels[va_idx].to(args.device)

    torch.manual_seed(0)
    net = HybridModel().to(args.device)
    opt = torch.optim.Adam(net.parameters(), lr=3e-4, weight_decay=1e-4)
    loss_fn = nn.BCEWithLogitsLoss()

    bs = 32
    for epoch in range(args.epochs):
        net.train()
        order = tr_idx[:]
        rng.shuffle(order)
        total = 0.0
        for start in range(0, len(order), bs):
            idx = order[start:start + bs]
            b = to_dev(collate([graphs[i] for i in idx]))
            y = labels[idx].to(args.device)
            opt.zero_grad()
            logits = net(b["h"], b["x"], b["edge_index"], b["edge_attr"],
                         b["batch"], b["seq"]).squeeze(-1)
            loss = loss_fn(logits, y)
            loss.backward()
            opt.step()
            total += loss.item()
        if epoch % 25 == 0 or epoch == args.epochs - 1:
            net.eval()
            with torch.no_grad():
                prob = torch.sigmoid(net(**va_batch).squeeze(-1)).cpu().numpy()
                y = y_va.cpu().numpy()
            print(f"epoch {epoch:3d} loss={total / max(1, len(order) // bs):.4f} "
                  f"val AUROC={roc_auc_score(y, prob):.3f} "
                  f"AUPRC={average_precision_score(y, prob):.3f}", flush=True)

    net.eval()
    for dtype, name in (("disp", "displaced"), ("<-", "graft")):
        sel = [i for i in va_idx if dtype in pdbs[i]] + [i for i in va_idx if labels[i] == 1]
        b = to_dev(collate([graphs[i] for i in sel]))
        with torch.no_grad():
            p = torch.sigmoid(net(**b).squeeze(-1)).cpu().numpy()
        y = labels[sel].numpy()
        print(f"val [{name}] AUROC={roc_auc_score(y, p):.3f} "
              f"AUPRC={average_precision_score(y, p):.3f} (n={len(sel)})", flush=True)

    out = Path(args.data).with_suffix(".hybrid.model.pt")
    torch.save({"model": net.state_dict()}, out)
    print(f"saved to {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
