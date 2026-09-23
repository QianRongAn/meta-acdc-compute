"""Average per-seed ensemble score tables into one cross-instance ensemble.

Each input is a `*_ensemble.tsv` (job_id, n_models, score, std). The output
has one row per job with the mean score across seeds and the across-seed
standard deviation — the fixed multi-instance ensemble used for all reported
rankings (the stability gate requires it; single instances are unreliable).

Usage:
    .venv/bin/python src/meta_acdc/structure/ensemble_mean.py \
        --scores data/processed/da_scores_seed0_ensemble.tsv \
                 data/processed/da_scores_seed1_ensemble.tsv \
                 data/processed/da_scores_seed2_ensemble.tsv \
        --out data/processed/prediction_scores_ensemble.tsv
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/prediction_scores_ensemble.tsv"))
    args = ap.parse_args()

    agg: dict[str, list[float]] = {}
    for p in args.scores:
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("score"):
                    agg.setdefault(r["job_id"], []).append(float(r["score"]))

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["job_id", "score", "std"])
        for job in sorted(agg):
            v = agg[job]
            w.writerow([job, f"{statistics.mean(v):.4f}",
                        f"{statistics.pstdev(v):.4f}"])
    print(f"averaged {len(args.scores)} seeds over {len(agg)} jobs "
          f"-> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())