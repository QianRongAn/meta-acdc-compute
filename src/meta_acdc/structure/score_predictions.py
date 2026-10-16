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
from meta_acdc.models.train_egnn import normalize_node_features
from meta_acdc.structure.cif import parse_cif, parse_all_models
from meta_acdc.structure.graph import build_interface_graph_from_residues


def graph_to_tensors(g):
    return {
        "h": torch.tensor(g.node_features, dtype=torch.float32),
        "peptide_mask": torch.tensor(g.peptide_mask, dtype=torch.bool),
        "x": torch.tensor(g.node_coords, dtype=torch.float32),
        "edge_index": torch.tensor(g.edge_index, dtype=torch.long).t(),
        "edge_attr": torch.tensor(g.edge_features, dtype=torch.float32),
        "batch": torch.zeros(g.n_nodes, dtype=torch.long),
    }


def score_graph(net: EGNN, g, device: str, mask_plddt: bool = False,
                temperature: float = 1.0) -> float:
    tensors = {k: v.to(device) for k, v in graph_to_tensors(g).items()}
    h = tensors["h"][:, :getattr(net,"n_node_dim",tensors["h"].shape[1])]
    tensors["edge_attr"]=tensors["edge_attr"][:,:getattr(net,"n_edge_dim",tensors["edge_attr"].shape[1])]
    if mask_plddt:
        # domain-adapted models are trained with pLDDT zeroed (col 3) —
        # scoring must apply the same masking to avoid distribution mismatch
        h = h.clone()
        h[:, 3] = 0.0
    # per-graph normalization, peptide flag protected — same transform as
    # training (train_egnn.normalize_node_features), so a scored graph gets
    # identical features regardless of what else is in the batch (AUDIT.md #3)
    h = normalize_node_features(h)
    with torch.no_grad():
        logit = net(h, tensors["x"], tensors["edge_index"],
                    tensors["edge_attr"], tensors["batch"], peptide_mask=tensors["peptide_mask"]).squeeze(-1)
    return float(torch.sigmoid(logit / temperature).item())


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
    ap.add_argument("--temperature", type=float, default=1.0,
                    help="divide logits by T before sigmoid (from "
                         "models/calibrate.py) for calibrated probabilities")
    ap.add_argument("--qc", type=Path,
                    default=Path("data/processed/af3_qc.tsv"),
                    help="structure-QC table from structure/af3_qc.py. Since "
                         "2026-10-15 the gate is ENFORCED: jobs without a QC "
                         "row, or with a non-empty flag, are not scored, and "
                         "the reason is written into the output instead of a "
                         "score (AUDIT.md #26)")
    ap.add_argument("--no-qc-gate", action="store_true",
                    help="score even QC-missing/failed jobs (audit/debug only; "
                         "the result then has no structure-quality basis)")
    ap.add_argument("--release", type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists() or Path(str(args.out).replace(".tsv","_ensemble.tsv")).exists():raise ValueError("refusing to overwrite scores")
    if not __import__("math").isfinite(args.temperature) or args.temperature<=0:raise ValueError("temperature must be finite and positive")
    from meta_acdc.data.release import load_release, sha256
    manifest, _ = load_release(args.release)
    import json
    identity_records = {x['complex_id']:x for x in json.loads((args.release.parent/'structure_records.json').read_text())}
    if args.no_qc_gate:
        raise ValueError("QC bypass is not allowed for accepted release scoring")
    if not args.model:
        args.model = [Path("data/processed/egnn_dataset.model.pt")]

    nets = []; masking=[];dimensions=[]
    for mp in args.model:
        ckpt = torch.load(mp, weights_only=False)
        if ckpt.get("release_sha256") != sha256(args.release):
            raise ValueError("legacy or mismatched checkpoint; retrain on the accepted release first")
        net = EGNN(node_dim=ckpt["n_node_dim"], edge_dim=ckpt["n_edge_dim"],
                   depth=6, hidden=128).to(args.device)
        net.load_state_dict(ckpt["model"])
        net.eval()
        net.n_node_dim=ckpt["n_node_dim"];net.n_edge_dim=ckpt["n_edge_dim"]
        nets.append(net);masking.append(bool(ckpt["mask_plddt"]));dimensions.append((ckpt["n_node_dim"],ckpt["n_edge_dim"]))
    if len(set(dimensions))!=1 or len(set(masking))!=1:raise ValueError("ensemble preprocessing mismatch")
    args.mask_plddt=masking[0]
    print(f"ensemble: {len(nets)} models", flush=True)

    # --- structure-quality gate -------------------------------------------
    if not args.qc.exists():
        raise ValueError("per-model QC table required")
    with args.qc.open(newline='') as fh:
        reader=csv.DictReader(fh,delimiter='\t')
        if not {'source_path','file_sha256','model_id','reference_id'} <= set(reader.fieldnames or []):
            raise ValueError("legacy job-level QC table refused; regenerate per-model QC")
        qc_status={}
        for row in reader:
            key=str(Path(row['source_path']).resolve())
            if key in qc_status:raise ValueError("duplicate per-model QC record")
            qc_status[key]=row
    rows = []
    for cif in sorted(args.cifs.glob("*.cif")):
        job_id = cif.stem
        job = job_id.rsplit("_model_", 1)[0]
        q=qc_status.get(str(cif.resolve()))
        failure='missing'
        if q:
            failure=q['flag']
            if sha256(cif)!=q['file_sha256']:failure='stale_structure_QC'
            if not q.get('confidence_path') or not Path(q['confidence_path']).is_file() or sha256(q['confidence_path'])!=q['confidence_sha256']:failure='stale_or_missing_confidence'
        if failure:
            rows.append({'job_id':job_id,'n_nodes':0,'n_edges':0,'score':'','std':'','qc':failure,'error':'per-model QC rejected'})
            continue
        try:
            residues = parse_cif(cif)
            from meta_acdc.structure.quality import annotate_prediction
            reference=identity_records[q['reference_id']]
            if q['reference_sha256']!=reference['source_sha256']:
                raise ValueError("QC reference identity version mismatch")
            annotate_prediction(residues, reference, job.split('_',1)[1])
            g = build_interface_graph_from_residues(residues,electrostatics=dimensions[0][0]>26 or dimensions[0][1]>16)
        except Exception as e:
            rows.append({"job_id": job_id, "n_nodes": 0, "n_edges": 0,
                         "score": "", "std": "",
                         "qc": "structure_revalidation_failed",
                         "error": str(e)})
            continue
        ss = [score_graph(net, g, args.device, mask_plddt=args.mask_plddt,
                          temperature=args.temperature)
              for net in nets]
        mean = sum(ss) / len(ss)
        std = (sum((s - mean) ** 2 for s in ss) / max(len(ss) - 1, 1)) ** 0.5
        rows.append({"job_id": job_id, "n_nodes": g.n_nodes,
                     "n_edges": len(g.edge_index), "score": f"{mean:.4f}",
                     "std": f"{std:.4f}",
                     "qc": "qc-pass" if not args.no_qc_gate else "not-checked",
                     "error": ""})
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
                         "score": f"{mean:.4f}", "std": f"{sd:.4f}", "qc": "qc-pass"})
        print(f"{job}: ensemble mean={mean:.4f} +/- {sd:.4f} over "
              f"{len(vals)} AF3 models", flush=True)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["job_id", "n_nodes", "n_edges",
                                           "score", "std", "qc", "error"],
                            delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    n_gated = sum(1 for r in rows if r.get("qc", "").startswith("qc-flag")
                  or r.get("qc") == "missing")
    agg_out = Path(str(args.out).replace(".tsv", "_ensemble.tsv"))
    with open(agg_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["job_id", "n_models", "score",
                                           "std", "qc"], delimiter="\t")
        w.writeheader()
        for r in agg_rows:
            w.writerow(r)
    print(f"saved {len(rows)} per-model scores to {args.out} and "
          f"{len(agg_rows)} ensemble scores to {agg_out}")
    if n_gated:
        print(f"structure-QC gate held back {n_gated} jobs (see the qc/error "
              f"columns). They are NOT part of any reported ranking.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
