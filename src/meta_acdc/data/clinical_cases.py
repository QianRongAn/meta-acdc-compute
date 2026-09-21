"""Gold-standard clinical cross-reactivity cases (proposal Module 3).

The classic MAGE-A3/Titin case is the proposal's designated gold-standard
validation target for the safety-scanning pipeline. This module curates the
known case data from the literature (lit review 02, section 4) into a
machine-readable test set for:
- AF3/TCRmodel2 batch 2 submissions (TCR x off-target peptide pairs)
- the in silico safety-scan benchmark (KN-11): the model must rank the
  known lethal off-targets (titin ESDPIVAQY) at the top.

Usage:
    .venv/bin/python src/meta_acdc/data/clinical_cases.py \
        --out data/processed/clinical_gold_standard.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# curated from lit review 02 section 4 (all entries PMID-verified there)
CASES = [
    {
        "tcr": "A3A (affinity-enhanced MAGE-A3 TCR, HLA-A*01)",
        "target": "EVDPIGHLY",  # MAGE-A3 168-176
        "off_targets": "ESDPIVAQY;ILAKFLHWL;KVAKELVHFL",
        # titin (LETHAL, Cameron 2013); MAGE-A6; MAGE-B18
        "fatal": "yes",
        "evidence": "PMID:23999400; PMID:23999405",
        "note": "titin ESDPIVAQY: 2 deaths (cardiac), backbone RMSD 0.285A vs target (Raman 2016, PDB 5BRZ)",
    },
    {
        "tcr": "A2-restricted MAGE-A3 TCR (A118T affinity-enhanced)",
        "target": "KVAELVHFL",  # MAGE-A3 112-120, HLA-A*02:01
        "off_targets": "KMVELVHFL;KMAELVHFL",
        # MAGE-A12 (LETHAL brain, Morgan 2013); MAGE-A9
        "fatal": "yes",
        "evidence": "PMID:23470321",
        "note": "MAGE-A12: 3 neurotoxicity / 2 deaths (necrotic leukoencephalopathy)",
    },
    {
        "tcr": "Tebentafusp/Kimmtrak (gp100 TCR, HLA-A*02:01)",
        "target": "YLEPGPVTA",  # gp100 280-288
        "off_targets": "YLEPGPVTV;YLEPGPVTL",
        # on-target/off-tumor (melanocytes) — approved, manageable
        "fatal": "no",
        "evidence": "PMID:34525232",
        "note": "FDA approved 2022; skin toxicity on-target/off-tumor, no Tx-related deaths",
    },
]

FIELD_ORDER = ["tcr", "target", "off_targets", "fatal", "evidence", "note"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/clinical_gold_standard.tsv"))
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELD_ORDER, delimiter="\t")
        w.writeheader()
        for c in CASES:
            w.writerow(c)
    print(f"gold-standard clinical cases: {len(CASES)} -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
