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


def load_seed(path: Path) -> dict[str, float]:
    d = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("score"):
                d[r["job_id"]] = float(r["score"])
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/prediction_scores_ensemble.tsv"))
    args = ap.parse_args()

    seeds = [load_seed(p) for p in args.scores]
    # degenerate-seed guard: a seed whose candidate scores are constant
    # carries no ranking signal and must not enter the reported ensemble
    good = []
    for p, d in zip(args.scores, seeds):
        vals = list(d.values())
        if len(vals) >= 2 and statistics.pstdev(vals) < 1e-6:
            print(f"excluding degenerate seed {p.name} (constant output)")
            continue
        good.append(d)
    if not good:
        print("all seeds degenerate", file=sys.stderr)
        return 1

    agg: dict[str, list[float]] = {}
    for d in good:
        for job, v in d.items():
            agg.setdefault(job, []).append(v)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["job_id", "score", "std"])
        for job in sorted(agg):
            v = agg[job]
            w.writerow([job, f"{statistics.mean(v):.4f}",
                        f"{statistics.pstdev(v):.4f}"])
    print(f"averaged {len(good)}/{len(seeds)} seeds over {len(agg)} jobs "
          f"-> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())