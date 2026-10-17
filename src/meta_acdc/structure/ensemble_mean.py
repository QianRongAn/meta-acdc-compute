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
                job = r["job_id"]
                value = float(r["score"])
                import math
                if not job or job in d or not math.isfinite(value):
                    raise ValueError("duplicate, empty, or nonfinite seed score")
                d[job] = value
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/prediction_scores_ensemble.tsv"))
    ap.add_argument("--min-coverage", type=int, default=None,
                    help="minimum number of seeds a job must be scored by to "
                         "appear in the output. Default: the full seed count "
                         "(the common job set). Output gains n_seeds and "
                         "uncertainty_known columns.")
    args = ap.parse_args()

    seeds = [load_seed(p) for p in args.scores]
    # Keep every requested seed. Constant outputs are failures of ranking,
    # not a reason to select a more agreeable survivor ensemble.
    good = seeds
    degenerate = [p.name for p,d in zip(args.scores,seeds)
                  if len(d)<2 or statistics.pstdev(list(d.values()))<1e-6]
    if degenerate:
        print(f"degenerate seeds retained; ensemble is diagnostic only: {degenerate}")
    if args.out.exists(): raise ValueError("refusing to overwrite ensemble")
    agg: dict[str, list[float]] = {}
    for d in good:
        for job, v in d.items():
            agg.setdefault(job, []).append(v)

    # A candidate scored by one seed and another by five must not both be
    # reported as an ensemble mean: the single-seed row has no uncertainty at
    # all, and a std of exactly 0 would read as "highly stable". Require the
    # common job set unless the caller opts out, and always report coverage.
    covered = [j for j, v in agg.items() if len(v) == len(good)]
    if args.min_coverage is not None:
        covered = [j for j, v in agg.items() if len(v) >= args.min_coverage]
    dropped = sorted(set(agg) - set(covered))
    if dropped:
        print(f"dropping {len(dropped)} jobs scored by fewer than "
              f"{args.min_coverage or len(good)} seeds; their 'std' would be "
              f"an artefact of missing coverage, not stability", flush=True)
    covered.sort()

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["job_id", "score", "std", "n_seeds", "uncertainty_known", "quality_eligible", "n_degenerate_seeds"])
        for job in covered:
            v = agg[job]
            known = "yes" if len(v) >= 2 else "no"   # 1 seed: std undefined, not 0
            sd = f"{statistics.pstdev(v):.4f}" if len(v) >= 2 else "NA"
            w.writerow([job, f"{statistics.mean(v):.4f}", sd, len(v), known, "no" if degenerate or len(v)!=len(seeds) else "yes", len(degenerate)])
    print(f"averaged {len(good)}/{len(seeds)} seeds over {len(covered)} jobs "
          f"({len(dropped)} dropped for insufficient coverage) -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())