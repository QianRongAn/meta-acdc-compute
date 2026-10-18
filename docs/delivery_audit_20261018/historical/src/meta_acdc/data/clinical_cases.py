"""Clinical cross-reactivity index cases (per-peptide evidence).

Curated from the literature. Two corrections after the 2026-10-15 audit
(AUDIT.md #6):

1. PMIDs are only written where they have been verified. The previous file
   carried PMID:23999405, which is a rare-earth electronic-structure paper and
   has nothing to do with TCR toxicity. Cameron et al. 2013 (titin, Sci Transl
   Med) is PMID 23926201. Unverified references are marked as such and must
   not be quoted as a PMID until checked.
2. Lethality is a property of a peptide, not of a case. "The patient died"
   was previously attached to the whole case, so every listed candidate
   inherited `fatal=yes` and was then flagged by a score margin. Each row now
   carries its own severity and its own evidence.

Structure identity is also explicit: 5BS0 is the MAGE-A3-reactive TCR
crystallised WITH the titin peptide (our own structure_vdjdb_map lists
5bs0 -> ESDPIVAQY), so scoring a predicted 5bs0/EVDPIGHLY complex is a
cross-structure check of the target on a related receptor, not an
independent second patient case.

Usage:
    .venv/bin/python src/meta_acdc/data/clinical_cases.py \
        --out data/processed/clinical_index_cases.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# One row per (case, peptide). Fields:
#   case, tcr, pdb, restriction, peptide, role, severity, evidence, note
CASES = [
    # --- A3A / MAGE-A3, HLA-A*01 (PDB 5BRZ, crystallised with MAGE-A3) ----
    dict(case="A3A-MAGEA3", tcr="A3A affinity-enhanced MAGE-A3 TCR",
         pdb="5brz", restriction="HLA-A*01", peptide="EVDPIGHLY",
         role="target", severity="approved-target",
         evidence="PMID:23926201 (Cameron 2013, Sci Transl Med)",
         note="MAGE-A3 168-176; affinity-enhanced TCR"),
    dict(case="A3A-MAGEA3", tcr="A3A affinity-enhanced MAGE-A3 TCR",
         pdb="5brz", restriction="HLA-A*01", peptide="ESDPIVAQY",
         role="off-target", severity="fatal",
         evidence="PMID:23926201 (Cameron 2013: 2 cardiac deaths via titin)",
         note="titin; backbone RMSD 0.285 A vs target (Raman 2016)"),
    dict(case="A3A-MAGEA3", tcr="A3A affinity-enhanced MAGE-A3 TCR",
         pdb="5brz", restriction="HLA-A*01", peptide="ILAKFLHWL",
         role="off-target", severity="not-documented",
         evidence="unverified - do not cite a PMID",
         note="listed candidate; no published fatality evidence found"),
    dict(case="A3A-MAGEA3", tcr="A3A affinity-enhanced MAGE-A3 TCR",
         pdb="5brz", restriction="HLA-A*01", peptide="KVAKELVHFL",
         role="off-target", severity="not-documented",
         evidence="unverified - do not cite a PMID",
         note="listed candidate; no published fatality evidence found"),
    # --- same peptide, related receptor crystallised WITH titin ------------
    dict(case="A3A-MAGEA3-cross-structure", tcr="MAGE-A3-reactive TCR",
         pdb="5bs0", restriction="HLA-A*01", peptide="EVDPIGHLY",
         role="target(predicted-complex)", severity="approved-target",
         evidence="cross-structure check only",
         note="5BS0 is crystallised with the titin peptide "
              "(structure_vdjdb_map: 5bs0 -> ESDPIVAQY); the MAGE-A3 target "
              "is a PREDICTED complex on this receptor, not a second case"),
    # --- A118T / MAGE-A12, HLA-A*02:01 ------------------------------------
    # NOTE: no crystal structure of the A118T receptor exists; 5BRZ is the
    # A3A receptor in HLA-A*01 (RCSB-confirmed). Scoring the A118T case with
    # 5BRZ structures would compare peptides against the WRONG receptor and
    # MHC, so pdb is left empty — clinical_scan reports "no structure" and
    # does not score it (AUDIT.md #6/#23).
    dict(case="A118T-MAGEA12", tcr="A118T affinity-enhanced MAGE-A3 TCR",
         pdb="", restriction="HLA-A*02:01", peptide="KVAELVHFL",
         role="target", severity="not-documented",
         evidence="unverified - cite as Linette 2013 (MAGE-A12 neurotoxicity)",
         note="MAGE-A3 112-120; no structure of this receptor available"),
    dict(case="A118T-MAGEA12", tcr="A118T affinity-enhanced MAGE-A3 TCR",
         pdb="", restriction="HLA-A*02:01", peptide="KMVELVHFL",
         role="off-target", severity="fatal-neurotoxicity",
         evidence="unverified - cite as Linette 2013 (MAGE-A12 neurotoxicity)",
         note="3 neurotoxicity / 2 deaths reported; no structure of this "
              "receptor available"),
    # --- Tebentafusp / Kimmtrak, gp100, HLA-A*02:01 -----------------------
    dict(case="KIMMTRAK-GP100", tcr="gp100 TCR (tebentafusp)",
         pdb="5eu6", restriction="HLA-A*02:01", peptide="YLEPGPVTA",
         role="target", severity="approved-drug",
         evidence="unverified - cite as tebentafusp pivotal trial",
         note="gp100 280-288; 5EU6 is crystallised with YLEPGPVTV"),
    dict(case="KIMMTRAK-GP100", tcr="gp100 TCR (tebentafusp)",
         pdb="5eu6", restriction="HLA-A*02:01", peptide="YLEPGPVTV",
         role="crystallised-peptide", severity="approved-drug",
         evidence="unverified - cite as tebentafusp pivotal trial",
         note="peptide present in the 5EU6 structure"),
    dict(case="KIMMTRAK-GP100", tcr="gp100 TCR (tebentafusp)",
         pdb="5eu6", restriction="HLA-A*02:01", peptide="YLEPGPVTL",
         role="off-target", severity="not-documented",
         evidence="unverified - do not cite a PMID",
         note="mimic; on-target/off-tumour skin toxicity reported, "
              "no treatment-related deaths"),
]

for _case in CASES:
    if "unverified" in _case["evidence"].lower():
        _case["pdb"] = ""

FIELD_ORDER = ["case", "tcr", "pdb", "restriction", "peptide", "role",
               "severity", "evidence", "note"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/clinical_index_cases.tsv"))
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELD_ORDER, delimiter="\t")
        w.writeheader()
        for c in CASES:
            w.writerow(c)
    print(f"clinical index cases: {len(CASES)} rows -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
