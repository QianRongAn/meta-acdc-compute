"""Clinical index-case scan (structural space).

Scores each curated (case, peptide) row against the domain-adapted ensemble
and reports the ranking inside the case's structure family. Rewritten after
the 2026-10-15 audit (AUDIT.md #6, #23):

- The structure family comes from the row's `pdb` column. Selecting "the
  first family that contains this target peptide" silently mixed up
  receptors: two different TCRs recognise the same peptide, and reordering
  the score file changed which case got a FATAL flag.
- Severity is per peptide. `fatal=yes` used to be a case-level column that
  every listed candidate inherited, after which any candidate within an
  arbitrary 0.1 margin was labelled "FATAL-mimicry at target level".
- Proximity to the target is reported as a number, never as "statistically
  indistinguishable from the target". A score difference is not a hypothesis
  test.

Only rows whose peptide was actually scored are reported; the rest are listed
under `missing`.

Usage:
    .venv/bin/python src/meta_acdc/structure/clinical_scan.py \
        --scores data/processed/prediction_scores_ensemble.tsv \
        --clinical data/processed/clinical_index_cases.tsv \
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
                    default=Path("data/processed/clinical_index_cases.tsv"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/kn11_clinical_scan.tsv"))
    args = ap.parse_args()

    # (pdb, peptide) -> score, preserving every row (no dict-key overwrite)
    scores: dict[tuple[str, str], float] = {}
    with open(args.scores, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r.get("score"):
                continue
            job = r["job_id"]
            if "_" not in job:
                continue
            pdb, pep = job.split("_", 1)
            scores[(pdb, pep)] = float(r["score"])

    rows, missing = [], []
    with open(args.clinical, newline="", encoding="utf-8") as fh:
        cases = list(csv.DictReader(fh, delimiter="\t"))

    for c in cases:
        key = (c["pdb"], c["peptide"])
        target_key = (c["pdb"],
                      next((r["peptide"] for r in cases
                            if r["case"] == c["case"] and r["role"] == "target"),
                           c["peptide"]))
        if key not in scores:
            missing.append((c["case"], c["pdb"], c["peptide"]))
            continue
        s = scores[key]
        t = scores.get(target_key)
        delta = "" if (t is None or key == target_key) else f"{s - t:+.4f}"
        # Proximity is descriptive only, and only meaningful against a
        # different peptide. No margin, no "indistinguishable".
        proximity = ""
        if t is not None and key != target_key:
            proximity = ("at-or-above-target" if s >= t else "below-target")
        severity = c["severity"]
        hazard = ""
        if (severity.startswith("fatal") and t is not None
                and key != target_key and s >= t):
            hazard = "lethal-mimic scores at/above the intended target"
        rows.append({
            "case": c["case"], "tcr": c["tcr"], "pdb": c["pdb"],
            "restriction": c["restriction"], "peptide": c["peptide"],
            "role": c["role"], "severity": severity, "evidence": c["evidence"],
            "score": f"{s:.4f}", "target_score": "" if t is None else f"{t:.4f}",
            "delta_vs_target": delta, "proximity": proximity,
            "hazard_flag": hazard,
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=list(rows[0].keys()) if rows else
                           ["case", "pdb", "peptide", "score"])
        w.writeheader()
        w.writerows(rows)

    print(f"clinical index scan: {len(rows)} scored rows, "
          f"{len(missing)} not scored -> {args.out}")
    for r in rows:
        print(f"  {r['case']:26s} {r['peptide']:12s} {r['role']:22s} "
              f"{r['score']}  {r['proximity']}")
    if missing:
        print("not scored (no AF3 job for this pdb/peptide):")
        for case, pdb, pep in missing:
            print(f"  {case} {pdb}_{pep}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
