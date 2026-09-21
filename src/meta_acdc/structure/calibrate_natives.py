"""Score all native crystal complexes -> calibration distribution.

Purpose: the EGNN was trained on native-vs-decoy discrimination; its
absolute scores on unseen natives vary by complex (e.g., 1oga crystal
0.741 vs typical 0.9+). This script scores every native crystal in the
structure set, producing the reference distribution used to interpret
AF3-predicted structure scores (e.g., z-score relative to crystal natives).

Usage:
    .venv/bin/python src/meta_acdc/structure/calibrate_natives.py \
        --struct-dir data/raw/structures \
        --model data/processed/egnn_dataset.model.pt \
        --out data/processed/native_crystal_scores.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import torch

from meta_acdc.models.egnn import EGNN
from meta_acdc.structure.graph import build_interface_graph


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct-dir", type=Path, default=Path("data/raw/structures"))
    ap.add_argument("--model", type=Path,
                    default=Path("data/processed/egnn_dataset.model.pt"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/native_crystal_scores.tsv"))
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    ckpt = torch.load(args.model, weights_only=False)
    net = EGNN(node_dim=ckpt["n_node_dim"], edge_dim=ckpt["n_edge_dim"],
               depth=6, hidden=128).to(args.device)
    net.load_state_dict(ckpt["model"])
    net.eval()

    with open(args.struct_dir / "complexes.tsv") as fh:
        complexes = [r["pdb"] for r in csv.DictReader(fh, delimiter="\t")]

    rows = []
    for pdb in complexes:
        path = args.struct_dir / f"{pdb.lower()}.pdb"
        try:
            g = build_interface_graph(path)
        except Exception:
            continue
        h = torch.tensor(g.node_features, dtype=torch.float32)
        h = (h - h.mean(0, keepdim=True)) / (h.std(0, keepdim=True) + 1e-6)
        x = torch.tensor(g.node_coords)
        ei = torch.tensor(g.edge_index, dtype=torch.long).t()
        ea = torch.tensor(g.edge_features)
        b = torch.zeros(g.n_nodes, dtype=torch.long)
        with torch.no_grad():
            s = torch.sigmoid(net(h.to(args.device), x.to(args.device),
                                  ei.to(args.device), ea.to(args.device),
                                  b.to(args.device)).squeeze(-1)).item()
        rows.append((pdb, s, g.n_nodes))

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pdb", "score", "n_nodes"])
        for r in rows:
            w.writerow(r)

    scores = [s for _, s, _ in rows]
    mean = sum(scores) / len(scores)
    sd = (sum((s - mean) ** 2 for s in scores) / max(len(scores) - 1, 1)) ** 0.5
    print(f"native crystals scored: {len(rows)}")
    print(f"score distribution: mean={mean:.3f} +/- {sd:.3f} "
          f"min={min(scores):.3f} max={max(scores):.3f}")
    print(f"saved to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
