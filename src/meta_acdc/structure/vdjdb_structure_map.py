"""Map structure TCRs to VDJdb records (KN-4+).

Strategy: extract the full ATOM-record sequence of TCR chains A/B from each
complex PDB, then search every VDJdb CDR3 as an exact substring of those
chains (CDR3 is contained in the solved V domain). Hits link a structure to
real specificity data:
- match: PDB peptide == VDJdb epitope -> strong positive (same complex class)
- mismatch: TCR validated against a DIFFERENT epitope -> cross-reactivity /
  multi-specificity evidence.

Usage:
    .venv/bin/python src/meta_acdc/structure/vdjdb_structure_map.py \
        --struct-dir data/raw/structures --vdjdb data/processed/vdjdb.clean.tsv \
        --out data/processed/structure_vdjdb_map.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

from meta_acdc.structure.graph import classify_chains, parse_pdb, TCR, PEPTIDE

CDR3_MIN_LEN = 8


def chain_sequences(path: Path) -> dict[str, str]:
    """Full ATOM sequence per chain.

    Appends CA residues in FILE ORDER — crucial for CDR3 regions, where
    insertion codes (resid 100A/100B) mean consecutive residues can share a
    number. Only consecutive duplicate (resid, icode) pairs are skipped
    (alternative conformations).
    """
    last: dict[str, tuple[int, str]] = {}
    seq: dict[str, list[str]] = defaultdict(list)
    with open(path, errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21].strip()
            resid = int(line[22:26])
            icode = line[26].strip()
            if last.get(chain) == (resid, icode):
                continue
            last[chain] = (resid, icode)
            seq[chain].append(line[17:20].strip())
    return {c: "".join(s) for c, s in seq.items()}


def aa3to1(seq: str) -> str:
    table = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
        "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
        "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M",
    }
    out = []
    for i in range(0, len(seq) - 2, 3):
        out.append(table.get(seq[i:i + 3], "X"))
    return "".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct-dir", type=Path, default=Path("data/raw/structures"))
    ap.add_argument("--vdjdb", type=Path, default=Path("data/processed/vdjdb.clean.tsv"))
    ap.add_argument("--out", type=Path, default=Path("data/processed/structure_vdjdb_map.tsv"))
    args = ap.parse_args()

    # load VDJdb: cdr3 -> set of epitopes (binders only)
    cdr3_epitopes: dict[str, set[str]] = defaultdict(set)
    with open(args.vdjdb, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["label"] == "1" and row["cdr3_beta"] and row["peptide"]:
                cdr3_epitopes[row["cdr3_beta"]].add(row["peptide"])
    print(f"VDJdb positive CDR3s: {len(cdr3_epitopes)}", flush=True)

    # complexes
    with open(args.struct_dir / "complexes.tsv") as fh:
        complexes = list(csv.DictReader(fh, delimiter="\t"))
    print(f"complexes to map: {len(complexes)}", flush=True)

    matches = []
    for row in complexes:
        pdb = row["pdb"]
        path = args.struct_dir / f"{pdb.lower()}.pdb"
        residues = parse_pdb(path)
        roles = classify_chains(residues)
        # PDB peptide sequence: chain(s) with the peptide role
        pep_chain = next((c for c, r in roles.items() if r == PEPTIDE), None)
        pdb_pep = ""
        if pep_chain:
            pep_res = sorted(
                {(r.resid, r.icode): r.resname
                 for r in residues if r.chain == pep_chain}.items()
            )
            pdb_pep = aa3to1("".join(n for _, n in pep_res))

        chains = chain_sequences(path)
        tcr_seqs = {c: aa3to1(s) for c, s in chains.items()
                    if roles.get(c) == TCR}
        for cdr3, epitopes in cdr3_epitopes.items():
            if len(cdr3) < CDR3_MIN_LEN:
                continue
            for chain, seq in tcr_seqs.items():
                # locate CDR3 window: after the conserved YFC/YLC motif,
                # bounded by the FGXG motif (length cap 30)
                m = None
                for motif in ("YFC", "YLC"):
                    pos = seq.find(motif)
                    if pos != -1:
                        m = (pos, motif)
                        break
                if m is None:
                    continue
                # IMGT CDR3 starts at the conserved Cys (the C of YFC);
                # require the VDJdb CDR3 to start exactly there (idx == 0)
                # to exclude substring-overlap false positives
                start = m[0] + 2
                window = seq[start:start + 35]
                idx = window.find(cdr3)
                if idx != 0:
                    continue
                for epitope in epitopes:
                    matches.append({
                        "pdb": pdb,
                        "chain": chain,
                        "cdr3": cdr3,
                        "vdjdb_epitope": epitope,
                        "pdb_peptide": pdb_pep,
                        "match": "same" if epitope == pdb_pep else "different",
                    })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pdb", "chain", "cdr3", "vdjdb_epitope",
                                           "pdb_peptide", "match"], delimiter="\t")
        w.writeheader()
        for m in matches:
            w.writerow(m)

    n_pdb = len({m["pdb"] for m in matches})
    n_same = sum(1 for m in matches if m["match"] == "same")
    n_diff = sum(1 for m in matches if m["match"] == "different")
    print(f"mapped: {len(matches)} CDR3 hits across {n_pdb} PDBs "
          f"({n_same} same-epitope, {n_diff} different-epitope)", flush=True)
    print(f"saved to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
