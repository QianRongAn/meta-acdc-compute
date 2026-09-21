"""Statistical tests for the per-TCR ranking experiment.

When N native controls are scored (one per structure), test whether the
model ranks native peptides significantly better than chance:
- pooled exact test: P(rank <= observed) under random permutation of the
  candidate order, pooled over structures with the same candidate count
- per-structure binomial-style test

Currently runnable with any number of native controls; meaningful from
N >= 3-4. Rerun as user submissions arrive.

Usage:
    .venv/bin/python src/meta_acdc/structure/rank_stats.py \
        --ensemble data/processed/prediction_scores_ensemble.tsv \
        --map data/processed/structure_vdjdb_map.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ensemble", type=Path,
                    default=Path("data/processed/prediction_scores_ensemble.tsv"))
    ap.add_argument("--map", type=Path,
                    default=Path("data/processed/structure_vdjdb_map.tsv"))
    args = ap.parse_args()

    scores: dict[tuple[str, str], float] = {}
    with open(args.ensemble, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pdb, pep = row["job_id"].split("_", 1)
            scores[(pdb, pep)] = float(row["score"])

    native: dict[str, str] = {}
    with open(args.map, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            native.setdefault(row["pdb"], row["pdb_peptide"])

    by_pdb: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for (pdb, pep), s in scores.items():
        by_pdb[pdb].append((pep, s))

    rows = []
    for pdb, entries in sorted(by_pdb.items()):
        nat = native.get(pdb, "")
        if not nat:
            continue
        entries.sort(key=lambda r: -r[1])
        rank = next((i + 1 for i, (p, _) in enumerate(entries) if p == nat), None)
        if rank is None:
            continue
        n = len(entries)
        rows.append((pdb, nat, rank, n))

    print(f"native controls scored: {len(rows)}")
    for pdb, nat, rank, n in rows:
        p = rank / n
        print(f"  {pdb:8s} native={nat:12s} rank={rank}/{n} "
              f"(one-sided P={p:.3f})")
    if len(rows) >= 2:
        # pooled: product of one-sided p-values via Fisher's method
        chi2 = -2 * sum(math.log(r / n) for _, _, r, n in rows)
        df = 2 * len(rows)
        # approximate p via chi2 survival (no scipy dependency)
        from statistics import NormalDist
        z = (chi2 - df) / math.sqrt(2 * df)
        pooled_p = 1 - NormalDist().cdf(z)
        print(f"\nFisher pooled P (native ranks not better than chance): "
              f"{pooled_p:.4f} (chi2={chi2:.1f}, df={df})")
        print("interpretation: P < 0.05 => model ranks native peptides "
              "significantly above chance across structures")
    else:
        print("\nneed >= 2 native controls for the pooled test; "
              "awaiting 1qse/1qsf submissions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
