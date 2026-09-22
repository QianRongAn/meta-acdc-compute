"""Score predicted TCR-pMHC structures with the trained EGNN (KN-5 pipeline).

Input: a directory of AF3/TCRmodel2 CIF files (one per complex).
Output: a risk table — for each job, the EGNN "compatibility" score
(sigmoid of the logit; the model was trained to output high scores for
native interfaces and low scores for incompatible ones).

Usage:
    .venv/bin/python src/meta_acdc/structure/score_predictions.py \
        --cifs data/raw/af3_predictions \
        --model data/processed/egnn_dataset.model.pt \
        --out data/processed/prediction_scores.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import torch

from meta_acdc.models.egnn import EGNN
from meta_acdc.structure.cif import parse_cif, parse_all_models
from meta_acdc.structure.graph import build_interface_graph_from_residues


def graph_to_tensors(g):
    return {
        "h": torch.tensor(g.node_features, dtype=torch.float32),
        "x": torch.tensor(g.node_coords, dtype=torch.float32),
        "edge_index": torch.tensor(g.edge_index, dtype=torch.long).t(),
        "edge_attr": torch.tensor(g.edge_features, dtype=torch.float32),
        "batch": torch.zeros(g.n_nodes, dtype=torch.long),
    }


def score_graph(net: EGNN, g, device: str, mask_plddt: bool = False) -> float:
    tensors = {k: v.to(device) for k, v in graph_to_tensors(g).items()}
    # normalize node features as in training (per-batch z-score)
    h = tensors["h"]
    if mask_plddt:
        # domain-adapted models are trained with pLDDT zeroed (col 3) —
        # scoring must apply the same masking to avoid distribution mismatch
        h = h.clone()
        h[:, 3] = 0.0
    h = (h - h.mean(0, keepdim=True)) / (h.std(0, keepdim=True) + 1e-6)
    with torch.no_grad():
        logit = net(h, tensors["x"], tensors["edge_index"],
                    tensors["edge_attr"], tensors["batch"]).squeeze(-1)
    return float(torch.sigmoid(logit).item())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cifs", type=Path, default=Path("data/raw/af3_predictions"))
    ap.add_argument("--model", type=Path, action="append", default=[],
                    help="model checkpoint; repeat for ensemble")
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/prediction_scores.tsv"))
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--mask-plddt", action="store_true",
                    help="zero the pLDDT node feature (domain-adapted "
                         "models were trained this way)")
    args = ap.parse_args()
    if not args.model:
        args.model = [Path("data/processed/egnn_dataset.model.pt")]

    nets = []
    for mp in args.model:
        ckpt = torch.load(mp, weights_only=False)
        net = EGNN(node_dim=ckpt["n_node_dim"], edge_dim=ckpt["n_edge_dim"],
                   depth=6, hidden=128).to(args.device)
        net.load_state_dict(ckpt["model"])
        net.eval()
        nets.append(net)
    print(f"ensemble: {len(nets)} models", flush=True)

    rows = []
    for cif in sorted(args.cifs.glob("*.cif")):
        try:
            residues = parse_cif(cif)
            g = build_interface_graph_from_residues(residues)
        except Exception as e:
            rows.append({"job_id": cif.stem, "n_nodes": 0, "n_edges": 0,
                         "score": "", "std": "", "error": str(e)})
            continue
        ss = [score_graph(net, g, args.device, mask_plddt=args.mask_plddt)
              for net in nets]
        mean = sum(ss) / len(ss)
        std = (sum((s - mean) ** 2 for s in ss) / max(len(ss) - 1, 1)) ** 0.5
        rows.append({"job_id": cif.stem, "n_nodes": g.n_nodes,
                     "n_edges": len(g.edge_index), "score": f"{mean:.4f}",
                     "std": f"{std:.4f}", "error": ""})
        print(f"{cif.stem}: score={mean:.4f} +/- {std:.4f} "
              f"(nodes={g.n_nodes})", flush=True)

    # aggregate per job across AF3 models (model_0..4): mean + spread
    from collections import defaultdict
    agg: dict[str, list[float]] = defaultdict(list)
    for r in rows:
        if r["score"]:
            job = r["job_id"].rsplit("_model_", 1)[0]
            agg[job].append(float(r["score"]))
    agg_rows = []
    for job, vals in sorted(agg.items()):
        mean = sum(vals) / len(vals)
        sd = (sum((v - mean) ** 2 for v in vals) / max(len(vals) - 1, 1)) ** 0.5
        agg_rows.append({"job_id": job, "n_models": len(vals),
                         "score": f"{mean:.4f}", "std": f"{sd:.4f}"})
        print(f"{job}: ensemble mean={mean:.4f} +/- {sd:.4f} over "
              f"{len(vals)} AF3 models", flush=True)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["job_id", "n_nodes", "n_edges",
                                           "score", "std", "error"],
                           delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    agg_out = Path(str(args.out).replace(".tsv", "_ensemble.tsv"))
    with open(agg_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["job_id", "n_models", "score", "std"],
                           delimiter="\t")
        w.writeheader()
        for r in agg_rows:
            w.writerow(r)
    print(f"saved {len(rows)} per-model scores to {args.out} and "
          f"{len(agg_rows)} ensemble scores to {agg_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
