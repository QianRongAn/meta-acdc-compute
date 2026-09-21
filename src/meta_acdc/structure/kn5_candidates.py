"""Generate the KN-5 structure-prediction submission list.

From the structure-VDJdb map, build (TCR alpha/beta sequences, peptide,
MHC hint) rows for TCRmodel2 / AlphaFold submissions. Priority:
1. cross-reactivity pairs (different-epitope hits: the TCR validated against
   a peptide OTHER than the one in its solved structure)
2. same-epitope pairs (control)

The TCR full sequences come from the mapped PDB file; the peptide is the
VDJdb epitope; the MHC hint is scraped from COMPND/TITLE records.

Usage:
    .venv/bin/python src/meta_acdc/structure/kn5_candidates.py \
        --map data/processed/structure_vdjdb_map.tsv \
        --struct-dir data/raw/structures \
        --out data/processed/kn5_submission_list.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from meta_acdc.structure.graph import classify_chains, parse_pdb, TCR, AA3TO1
from meta_acdc.structure.vdjdb_structure_map import chain_sequences, aa3to1


def mhc_hint(path: Path) -> str:
    """Scrape HLA allele mention from COMPND/TITLE records."""
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(("COMPND", "TITLE")):
                m = re.search(r"HLA-?[A-Z0-9]+\*?\d*:?\d*", line, re.I)
                if m:
                    return m.group(0).upper()
    return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", type=Path,
                    default=Path("data/processed/structure_vdjdb_map.tsv"))
    ap.add_argument("--struct-dir", type=Path, default=Path("data/raw/structures"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/kn5_submission_list.tsv"))
    args = ap.parse_args()

    hits = defaultdict(list)
    with open(args.map, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            hits[row["pdb"]].append(row)

    rows = []
    for pdb, hlist in sorted(hits.items()):
        path = args.struct_dir / f"{pdb.lower()}.pdb"
        if not path.exists():
            continue
        try:
            residues = parse_pdb(path)
            roles = classify_chains(residues)
            seqs = {c: aa3to1(s) for c, s in chain_sequences(path).items()}
        except Exception:
            continue
        tcr_seqs = {c: s for c, s in seqs.items() if roles.get(c) == TCR}
        if len(tcr_seqs) < 2:
            continue
        tcr_list = sorted(tcr_seqs.values(), key=len, reverse=True)
        tcr_a, tcr_b = tcr_list[0], tcr_list[1]
        hint = mhc_hint(path)
        for h in hlist:
            rows.append({
                "pdb": pdb,
                "tcr_a": tcr_a,
                "tcr_b": tcr_b,
                "peptide": h["vdjdb_epitope"],
                "mhc_hint": hint,
                "evidence": h["match"],
                "pdb_peptide": h["pdb_peptide"],
                "cdr3": h["cdr3"],
            })

    # priority: cross-reactivity (different) first, then same-epitope controls
    rows.sort(key=lambda r: (r["evidence"] != "different", r["pdb"]))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pdb", "tcr_a", "tcr_b", "peptide",
                                           "mhc_hint", "evidence",
                                           "pdb_peptide", "cdr3"],
                           delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    n_diff = sum(1 for r in rows if r["evidence"] == "different")
    n_same = sum(1 for r in rows if r["evidence"] == "same")
    n_pdb = len({r["pdb"] for r in rows})
    print(f"KN-5 submission list: {len(rows)} pairs ({n_diff} cross-reactivity, "
          f"{n_same} same-epitope) from {n_pdb} structures")
    print(f"saved to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
