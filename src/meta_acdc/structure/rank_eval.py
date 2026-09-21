"""Statistical evaluation of per-TCR ranking experiments (KN-5).

For each TCR with scored candidate peptides:
- recall@0.5: fraction of VDJdb-validated binders scoring above 0.5
- native rank: where the native peptide lands (should be #1-2)
- score spread: min/max for sanity

Usage:
    .venv/bin/python src/meta_acdc/structure/rank_eval.py \
        --scores data/processed/prediction_scores.tsv \
        --map data/processed/structure_vdjdb_map.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

THRESHOLD = 0.5


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path,
                    default=Path("data/processed/prediction_scores.tsv"))
    ap.add_argument("--map", type=Path,
                    default=Path("data/processed/structure_vdjdb_map.tsv"))
    args = ap.parse_args()

    scores: dict[str, float] = {}
    with open(args.scores, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["score"]:
                scores[row["job_id"].rsplit("_model", 1)[0]] = float(row["score"])

    # evidence: (pdb, peptide) -> set of cdr3s with VDJdb validation
    evidence: dict[tuple[str, str], int] = defaultdict(int)
    native: dict[str, str] = {}
    with open(args.map, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["pdb"], row["vdjdb_epitope"])
            evidence[key] += 1
            native.setdefault(row["pdb"], row["pdb_peptide"])

    by_pdb: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for job_id, s in scores.items():
        pdb, pep = job_id.split("_", 1)
        by_pdb[pdb].append((pep, s))

    print(f"{'TCR':8s} {'n':>3s} {'recall@0.5':>10s} {'native rank':>12s}  verdict")
    for pdb, rows in sorted(by_pdb.items()):
        rows.sort(key=lambda r: -r[1])
        n_val = sum(1 for p, _ in rows if evidence[(pdb, p)] > 0)
        n_high = sum(1 for p, s in rows
                     if s >= THRESHOLD and evidence[(pdb, p)] > 0)
        recall = n_high / max(n_val, 1)
        nat = native.get(pdb, "")
        nat_rank = next((i for i, (p, _) in enumerate(rows, 1) if p == nat), None)
        verdict = ""
        if nat_rank is not None:
            verdict = "native in top-2 ✓" if nat_rank <= 2 else "native ranked LOW ⚠"
        print(f"{pdb:8s} {len(rows):3d} {recall:9.1%} {str(nat_rank) if nat_rank else 'n/a':>11s}  {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
