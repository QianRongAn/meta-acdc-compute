"""AF3 resubmission variance: AF3 self-consistency vs EGNN scorer variance.

Question this answers: the legacy benchmark (docs/benchmarks.md, 2026-09-22)
claimed that re-submitting the same job to AlphaFold3 swings the reported
score by up to 0.94 (LLFGYPRYV 0.049 -> 0.991), and concluded "a single AF3
submission is untrustworthy". That conclusion conflated two very different
variance sources:

  1. AF3 run-to-run variance  -> measured here with AF3's OWN confidence
     metrics (ipTM/pTM), which are model-independent.
  2. EGNN scorer variance     -> measured by replaying each AF3 submission
     through the SAME (fixed) EGNN checkpoint.

Layout: data/raw/af3_instances/<instance>/{job_id}_model_{0..4}.cif,
        and confidences {job_id}_conf_0.json (B*) / _conf.json (A, legacy).

Usage:
    .venv/bin/python src/meta_acdc/structure/af3_variance.py \
        --train-root data/raw/af3_instances \
        --instances A B1 B2 B3 \
        --scorer egnn_noplddt_s10:data/processed/af3_scores_s10_{i}_ensemble.tsv \
        --scorer egnn_dataset:data/processed/af3_scores_{i}_ensemble.tsv \
        --out data/processed/af3_resubmission_variance.tsv
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


def af3_conf_metrics(inst_dir: Path, job_id: str) -> dict[str, list[float]]:
    """Collect ipTM/pTM across all AF3 model confidences for one submission."""
    out: dict[str, list[float]] = defaultdict(list)
    for f in sorted(inst_dir.glob(f"{job_id}_conf*.json")):
        with open(f) as fh:
            d = json.load(fh)
        for k in ("iptm", "ptm", "ranking_score"):
            if d.get(k) is not None:
                out[k].append(float(d[k]))
        cpm = d.get("chain_pair_iptm")
        if cpm and len(cpm) >= 5:  # chains: alpha,beta,peptide,mhc,b2m
            out["pep_mhc_iptm"].append(float(cpm[2][3]))
    return out


def load_egnn(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    d = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("score"):
                d[r["job_id"]] = float(r["score"])
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-root", type=Path,
                    default=Path("data/raw/af3_instances"))
    ap.add_argument("--instances", nargs="+", default=["A", "B1", "B2", "B3"])
    ap.add_argument("--scorer", action="append", default=[],
                    help="name:template where template contains {i} for the "
                         "instance; repeat for multiple scorers")
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/af3_resubmission_variance.tsv"))
    args = ap.parse_args()
    if not args.scorer:
        args.scorer = [
            "egnn_noplddt_s10:data/processed/af3_scores_s10_{i}_ensemble.tsv",
            "egnn_dataset:data/processed/af3_scores_{i}_ensemble.tsv",
        ]

    scorers = []
    for spec in args.scorer:
        name, tmpl = spec.split(":", 1)
        scorers.append((name, tmpl))

    # ---- AF3 self-consistency, per job per instance -----------------------
    af3: dict[str, dict[str, dict[str, list[float]]]] = {}
    for inst in args.instances:
        inst_dir = args.train_root / inst
        if not inst_dir.is_dir():
            continue
        jobs = sorted({os.path.basename(f).split("_conf")[0]
                       for f in glob.glob(str(inst_dir / "*_conf*.json"))})
        af3[inst] = {job: af3_conf_metrics(inst_dir, job) for job in jobs}

    # ---- fixed-scorer EGNN scores, per job per instance -------------------
    egnn: dict[str, dict[str, float]] = {}
    for name, tmpl in scorers:
        egnn[name] = {inst: load_egnn(Path(tmpl.format(i=inst)))
                      for inst in args.instances}

    all_jobs = sorted(set().union(*[set(d) for d in af3.values()]))
    rows = []
    for job in all_jobs:
        row = {"job_id": job}
        for inst in args.instances:
            m = af3.get(inst, {}).get(job)
            if m and "iptm" in m:
                row[f"iptm_{inst}"] = f"{np.mean(m['iptm']):.3f}"
                row[f"iptm_std_{inst}"] = f"{np.std(m['iptm']):.3f}"
            else:
                row[f"iptm_{inst}"] = ""
                row[f"iptm_std_{inst}"] = ""
        # AF3 across-submission spread (on model-mean ipTM)
        iptm_means = [np.mean(af3[i][job]["iptm"]) for i in args.instances
                      if job in af3.get(i, {}) and af3[i][job].get("iptm")]
        row["af3_iptm_range"] = (f"{max(iptm_means) - min(iptm_means):.3f}"
                                 if len(iptm_means) > 1 else "")
        for name, _ in scorers:
            vals = [egnn[name][i][job] for i in args.instances
                    if job in egnn[name].get(i, {})]
            for inst in args.instances:
                v = egnn[name].get(inst, {}).get(job)
                row[f"{name}_{inst}"] = f"{v:.4f}" if v is not None else ""
            row[f"{name}_range"] = (f"{max(vals) - min(vals):.4f}"
                                    if len(vals) > 1 else "")
        rows.append(row)

    fields = ["job_id"] + [f"{p}_{i}" for i in args.instances
                           for p in ("iptm", "iptm_std")] \
        + ["af3_iptm_range"] \
        + [f for name, _ in scorers
           for f in [f"{name}_{i}" for i in args.instances] + [f"{name}_range"]]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields,
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # ---- summary ---------------------------------------------------------
    print(f"wrote {len(rows)} jobs to {args.out}")
    ranges = [float(r["af3_iptm_range"]) for r in rows if r["af3_iptm_range"]]
    if ranges:
        print(f"\nAF3 self-consistency (model-mean ipTM across resubmissions):")
        print(f"  jobs with >1 submission : {len(ranges)}")
        print(f"  median |d ipTM|         : {np.median(ranges):.3f}")
        print(f"  max    |d ipTM|         : {max(ranges):.3f}")
    for name, _ in scorers:
        rs = [float(r[f"{name}_range"]) for r in rows if r.get(f"{name}_range")]
        if rs:
            print(f"\nEGNN scorer '{name}' resubmission range (same model):")
            print(f"  jobs with >1 submission : {len(rs)}")
            print(f"  median range            : {np.median(rs):.4f}")
            print(f"  max    range            : {max(rs):.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
