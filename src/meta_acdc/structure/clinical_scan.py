"""KN-11 in silico clinical safety scan (structural space).

With a cross-instance-stable domain-adapted model, score each clinical
gold-standard case's cognate target and its documented off-targets, and
report the compatibility ranking. The safety interpretation is conservative:
a *fatal* off-target that scores at or above the cognate target is a
correctly-flagged hazard (the model cannot distinguish a lethal mimicry
peptide from the intended target — which is precisely the danger).

Usage:
    .venv/bin/python src/meta_acdc/structure/clinical_scan.py \
        --scores data/processed/prediction_scores_ensemble.tsv \
        --clinical data/processed/clinical_gold_standard.tsv \
        --out data/processed/kn11_clinical_scan.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path,
                    default=Path("data/processed/prediction_scores_ensemble.tsv"))
    ap.add_argument("--clinical", type=Path,
                    default=Path("data/processed/clinical_gold_standard.tsv"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/kn11_clinical_scan.tsv"))
    ap.add_argument("--margin", type=float, default=0.1,
                    help="a fatal off-target scoring within this margin of the "
                         "cognate target is flagged 'at target level' (the "
                         "model cannot separate a lethal mimicry peptide from "
                         "the intended target — a conservative hazard signal)")
    args = ap.parse_args()

    # scores keyed by (pdb, peptide); also peptide -> list of scores
    by_pdb: dict[str, dict[str, float]] = {}
    with open(args.scores, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r.get("score"):
                continue
            pdb, pep = r["job_id"].split("_", 1)
            by_pdb.setdefault(pdb, {})[pep] = float(r["score"])

    rows = []
    with open(args.clinical, newline="", encoding="utf-8") as fh:
        for c in csv.DictReader(fh, delimiter="\t"):
            tcr = c["tcr"]
            # match the clinical case to a scored structure family by target
            fam = next((p for p, d in by_pdb.items()
                        if c["target"] in d), None)
            if fam is None:
                rows.append({"case": tcr, "family": "", "peptide": c["target"],
                             "role": "target", "score": "",
                             "fatal": c["fatal"], "flag": "no structure"})
                continue
            d = by_pdb[fam]
            entries = sorted(d.items(), key=lambda kv: -kv[1])
            rank_of = {p: i + 1 for i, (p, _) in enumerate(entries)}
            target = c["target"]
            offs = [o for o in c["off_targets"].split(";") if o]
            for pep in [target] + offs:
                if pep not in d:
                    continue
                role = "target" if pep == target else "off-target"
                at_or_above = d[pep] >= d[target] - args.margin
                flag = ""
                if role == "off-target" and c["fatal"] == "yes" \
                        and at_or_above:
                    flag = "FATAL-mimicry at target level"
                elif role == "target":
                    flag = "cognate"
                rows.append({
                    "case": tcr, "family": fam, "peptide": pep, "role": role,
                    "score": f"{d[pep]:.4f}", "rank": rank_of.get(pep, ""),
                    "fatal": c["fatal"], "flag": flag,
                })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["case", "family", "peptide", "role",
                                       "score", "rank", "fatal", "flag"])
        w.writeheader()
        w.writerows(rows)

    for r in rows:
        if r["role"] == "no structure":
            continue
        print(f"{r['family']:6s} {r['peptide']:12s} {r['role']:10s} "
              f"score={r['score']:>7s} rank={r.get('rank','')} "
              f"fatal={r['fatal']:3s} {r['flag']}")
    print(f"\nwrote {len(rows)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())