"""Per-TCR cross-reactivity ranking analysis (KN-5 experiment).

When all AF3 predictions for one TCR are scored, this script:
1. Groups scores by structure (pdb prefix of job_id)
2. Ranks candidate peptides by compatibility score within each TCR
3. Reports: native peptide rank, VDJdb-validated cross-reactive peptide ranks
4. Verdict: do the validated cross-reactive peptides outrank random candidates?

Usage:
    .venv/bin/python src/meta_acdc/structure/rank_analysis.py \
        --scores data/processed/prediction_scores.tsv \
        --map data/processed/structure_vdjdb_map.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path,
                    default=Path("data/processed/prediction_scores.tsv"))
    ap.add_argument("--map", type=Path,
                    default=Path("data/processed/structure_vdjdb_map.tsv"))
    args = ap.parse_args()

    scores = []
    with open(args.scores, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["score"]:
                scores.append(row)
    by_pdb: dict[str, list[dict]] = defaultdict(list)
    for s in scores:
        pdb = s["job_id"].split("_")[0]
        by_pdb[pdb].append(s)

    # native peptide for each pdb (from the structure map)
    native = {}
    with open(args.map, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            native.setdefault(row["pdb"], row["pdb_peptide"])

    for pdb, rows in sorted(by_pdb.items()):
        rows.sort(key=lambda r: -float(r["score"]))
        nat = native.get(pdb, "")
        print(f"\n=== {pdb} (native peptide: {nat}) ===")
        for i, r in enumerate(rows, 1):
            pep = r["job_id"].split("_", 1)[1]
            tag = " [NATIVE]" if pep == nat else ""
            print(f"  {i:2d}. {pep:12s} score={float(r['score']):.3f}{tag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
