"""Is the compatibility score just a peptide-similarity proxy?

For scored candidate peptides that have the same length as their cognate
crystal peptide, correlate the model score with sequence identity to the
cognate.

The output is a *screening* result only. A weak correlation does NOT show the
model "captures side-chain chemistry": identity to the cognate peptide is one
crude sequence proxy, and a null result can also come from structure quality,
receptor-family shifts, decoy artefacts, or noise. Supporting a mechanistic
claim needs controlled comparisons plus independent binding labels
(AUDIT.md #28). Two peptides at the same identity that score differently do
not rule out a subtler sequence shortcut either.

Usage:
    .venv/bin/python src/meta_acdc/structure/score_vs_identity.py \
        --scores data/processed/prediction_scores_ensemble.tsv \
        --map data/processed/structure_vdjdb_map.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from scipy.stats import spearmanr


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path,
                    default=Path("data/processed/prediction_scores_ensemble.tsv"))
    ap.add_argument("--map", type=Path,
                    default=Path("data/processed/structure_vdjdb_map.tsv"))
    args = ap.parse_args()

    native: dict[str, str] = {}
    with open(args.map, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            native.setdefault(r["pdb"], r["pdb_peptide"])

    scores: dict[str, float] = {}
    with open(args.scores, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("score"):
                scores[r["job_id"]] = float(r["score"])

    rows = []
    for job, s in scores.items():
        pdb, pep = job.split("_", 1)
        nat = native.get(pdb)
        if not nat or len(pep) != len(nat):
            continue
        ident = sum(a == b for a, b in zip(pep, nat)) / len(nat)
        rows.append((pdb, pep, nat, ident, s))

    print(f"same-length pairs: {len(rows)}")
    for pdb, pep, nat, ident, s in sorted(rows, key=lambda x: -x[3]):
        print(f"  {pdb:5s} {pep:13s} vs {nat:13s} ident={ident:.2f} "
              f"score={s:.3f}")
    if len(rows) >= 3:
        xs = [r[3] for r in rows]
        ys = [r[4] for r in rows]
        if len(set(xs)) > 1 and len(set(ys)) > 1:
            r, p = spearmanr(xs, ys)
            print(f"\nSpearman(identity, score) = {r:.3f} (P={p:.3g})")
            # The reading depends on the statistic, not on the null we hoped
            # for. An unconditional "no association detected" is no better
            # than the unconditional "captures chemistry" it replaced.
            if r > 0.5 and p < 0.05:
                print("Reading: a POSITIVE association with this identity "
                      "metric was detected. The score tracks similarity to "
                      "the cognate peptide on this metric.")
            elif r < -0.5 and p < 0.05:
                print("Reading: a NEGATIVE association with this identity "
                      "metric was detected.")
            elif p < 0.05:
                print("Reading: a weak association with this identity metric "
                      "was detected.")
            else:
                print("Reading: no association with THIS identity metric was "
                      "detected.")
            print("Either way this is a SCREENING result on one crude "
                  "similarity proxy. It does not show the model captures "
                  "side-chain chemistry, and it does not exclude other "
                  "sequence shortcuts, structure artefacts or noise. Do not "
                  "quote it as evidence of chemical specificity without "
                  "controlled experiments and independent binding labels.")
    return 0


if __name__ == "__main__":
    sys.exit(main())