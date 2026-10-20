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

NODE_DIM = 26
EDGE_DIM = 16


def collate(graphs: list[tuple]) -> dict:
    from meta_acdc.models.train_egnn import collate as collate_graphs
    from meta_acdc.structure.feature_schema import PEPTIDE_FLAG_INDEX, AA_START, AA_STOP, AA_ORDER
    batch = collate_graphs(graphs, NODE_DIM, EDGE_DIM)
    seqs = []
    for nf, *_ in graphs:
        pep = [row for row in nf if row[PEPTIDE_FLAG_INDEX] == 1]
        if not pep or len(pep) > PEPTIDE_MAXLEN:
            raise ValueError("unsupported or missing peptide length; refusing truncation")
        seq = torch.zeros(PEPTIDE_MAXLEN, N_AA)
        for i, row in enumerate(pep):
            aa = AA_ORDER[int(np.argmax(row[AA_START:AA_STOP]))]
            seq[i, AA_IDX[aa]] = 1
        seqs.append(seq)
    batch["seq"] = torch.stack(seqs)
    return batch


def pdb_group(name: str) -> str:
    return name.split("-disp")[0].split("<-")[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=None)
    ap.add_argument("--release", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    from meta_acdc.data.release import load_release, frozen_indices, sha256
    manifest, dataset = load_release(args.release, data_path=args.data)
    if args.out.exists(): raise ValueError("refusing to overwrite checkpoint")
    if args.out.resolve().is_relative_to(args.release.resolve().parent):
        raise ValueError("checkpoint must be outside immutable release")
    data = torch.load(dataset, weights_only=False)
    graphs, labels = data["graphs"], torch.tensor(data["labels"], dtype=torch.float32)
    pdbs = data["pdbs"]
    print(f"dataset: {len(graphs)} graphs", flush=True)

    split = frozen_indices(manifest, data)
    tr_idx, va_idx = split['train'], split['validation']
    rng = random.Random(0)
    print(f"split: train={len(tr_idx)} val={len(va_idx)}", flush=True)

    def to_dev(b: dict) -> dict:
        return {k: v.to(args.device) for k, v in b.items()}

    va_batch = to_dev(collate([graphs[i] for i in va_idx]))
    y_va = labels[va_idx].to(args.device)

    torch.manual_seed(0)
    net = HybridModel(node_dim=NODE_DIM, edge_dim=EDGE_DIM).to(args.device)
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
                         b["batch"], b["seq"], peptide_mask=b["peptide_mask"]).squeeze(-1)
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
        if len({float(labels[i]) for i in sel}) < 2:
            print("per-control AUROC unavailable: both classes required", flush=True)
            continue
        b = to_dev(collate([graphs[i] for i in sel]))
        with torch.no_grad():
            p = torch.sigmoid(net(**b).squeeze(-1)).cpu().numpy()
        y = labels[sel].numpy()
        print(f"val [{name}] AUROC={roc_auc_score(y, p):.3f} "
              f"AUPRC={average_precision_score(y, p):.3f} (n={len(sel)})", flush=True)

    out = args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"model": net.state_dict(), "release_sha256": sha256(args.release),
                "task": manifest["task"], "mask_plddt": bool(manifest.get("mask_plddt_required",False)),
                "n_node_dim": NODE_DIM, "n_edge_dim": EDGE_DIM}, out)
    print(f"saved to {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
