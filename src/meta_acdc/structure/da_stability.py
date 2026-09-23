"""Cross-instance stability analysis for domain-adapted models (Rashomon fix).

Joins per-seed AF3 prediction scores on job_id, computes pairwise Spearman
rank correlation across model instances — the reproducibility gate that the
original v9.1 models failed (r = 0.17 / -0.24 / NaN).

Usage:
    .venv/bin/python src/meta_acdc/structure/da_stability.py \
        --scores data/processed/da_scores_seed0_ensemble.tsv \
                  data/processed/da_scores_seed1_ensemble.tsv \
                  data/processed/da_scores_seed2_ensemble.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr


def load_jobs(path: Path) -> dict[str, float]:
    jobs = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["score"]:
                jobs[row["job_id"]] = float(row["score"])
    return jobs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path, nargs="+", required=True,
                    help="per-seed ensemble score TSVs (same jobs)")
    args = ap.parse_args()

    per_seed = [load_jobs(p) for p in args.scores]
    common = set(per_seed[0])
    for j in per_seed[1:]:
        common &= set(j)
    common = sorted(common)
    if len(common) < 5:
        print(f"too few common jobs ({len(common)}) for Spearman",
              file=sys.stderr)
        return 1
    print(f"common jobs: {len(common)}", flush=True)

    M = np.array([[j[job] for job in common] for j in per_seed])
    rs = []
    for i, j in combinations(range(len(per_seed)), 2):
        r, p = spearmanr(M[i], M[j])
        rs.append(r)
        print(f"  seed{i} vs seed{j}: Spearman r={r:.3f} (P={p:.3g})",
              flush=True)
    print(f"mean pairwise r = {np.mean(rs):.3f}", flush=True)

    # bootstrap CI over jobs: how firmly is the gate passed/failed?
    rng = np.random.default_rng(0)
    boot = []
    n = len(common)
    for _ in range(2000):
        idx = rng.integers(0, n, n)
        pair_rs = []
        for i, j in combinations(range(len(per_seed)), 2):
            a, b = M[i][idx], M[j][idx]
            if np.std(a) == 0 or np.std(b) == 0:
                continue
            pair_rs.append(spearmanr(a, b)[0])
        if pair_rs:
            boot.append(np.mean(pair_rs))
    if boot:
        lo, hi = np.percentile(boot, [2.5, 97.5])
        print(f"bootstrap 95% CI on mean r: [{lo:.3f}, {hi:.3f}] "
              f"(n={n} jobs, {len(boot)} resamples)", flush=True)

    # verdict gate: mean r >= 0.5 = rankings reproducible enough to report
    verdict = "STABLE (report rankings)" if np.mean(rs) >= 0.5 else \
        "UNSTABLE (single-instance anecdotes only)"
    print(f"verdict: {verdict}", flush=True)

    # per-TCR-group breakdown (jobs are named <pdb>_<peptide>)
    groups: dict[str, list[int]] = {}
    for k, job in enumerate(common):
        groups.setdefault(job.split("_", 1)[0], []).append(k)
    print("\nper-TCR-group Spearman (>=4 jobs):")
    for g, idx in sorted(groups.items()):
        if len(idx) < 4:
            continue
        grs = [spearmanr(M[i][idx], M[j][idx])[0]
               for i, j in combinations(range(len(per_seed)), 2)]
        print(f"  {g}: mean r={np.mean(grs):.3f} over {len(idx)} jobs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
