"""Receptor-family map for leakage-controlled splits.

Two PDBs that share the same TCR receptor (same alpha/beta chain sequences)
are not independent data points: putting one in train and the other in
validation leaks the receptor. This builds

    data/processed/pdb_family_map.tsv   pdb -> family_id, receptor_key

where family_id groups PDBs with identical receptor sequences, and receptor_key
is a stable hash of those sequences. train_egnn.py consumes it via
`--group-by family`.

Complex identity is handled separately (same PDB, crystal vs AF3-predicted of
the same complex): pdb_group() normalises those to one group.

Usage:
    .venv/bin/python src/meta_acdc/structure/family_map.py \
        --structures data/raw/structures --out data/processed/pdb_family_map.tsv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from meta_acdc.structure.graph import parse_pdb, classify_chains  # noqa: E402
from meta_acdc.structure.graph import AA3TO1, TCR, PEPTIDE  # noqa: E402


def receptor_sequences(residues, roles) -> tuple[str, str]:
    """(alpha_seq, beta_seq) of the TCR-classified chains only, sorted.

    Including peptide/MHC/B2M chains (or chain names) made the receptor key
    change when the peptide changed — the same TCR bound to a different
    peptide became a different "family", re-opening receptor leakage. Only
    chains classified as TCR contribute. (AUDIT.md #4.)
    """
    per_chain: dict[str, list[str]] = defaultdict(list)
    for r in residues:
        if roles.get(r.chain) != TCR:
            continue
        per_chain[r.chain].append(r.resname)
    seqs = []
    for chain, names in per_chain.items():
        seqs.append("".join(AA3TO1.get(n.upper(), "X") for n in names))
    return tuple(sorted(seqs))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", type=Path, default=Path("data/raw/structures"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/pdb_family_map.tsv"))
    args = ap.parse_args()

    rows = []
    for pdb in sorted(args.structures.glob("*.pdb")):
        try:
            residues = parse_pdb(pdb)
            roles = classify_chains(residues)
        except Exception as exc:                       # noqa: BLE001
            print(f"skip {pdb.name}: {exc}", file=sys.stderr)
            continue
        if not any(v == TCR for v in roles.values()):
            continue
        if not any(v == PEPTIDE for v in roles.values()):
            continue
        seqs = receptor_sequences(residues, roles)
        key = hashlib.sha1("|".join(seqs).encode()).hexdigest()[:12]
        rows.append((pdb.stem, key))

    # PDBs sharing a receptor key are one family. Families are numbered by the
    # first PDB id so the ids are stable across runs.
    by_key: dict[str, list[str]] = defaultdict(list)
    for pdb, key in rows:
        by_key[key].append(pdb)
    fam: dict[str, str] = {}
    for key in sorted(by_key):
        members = sorted(by_key[key])
        fam_id = f"F{members[0]}"
        for m in members:
            fam[m] = fam_id

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pdb", "family_id", "receptor_key"])
        for pdb, key in sorted(rows, key=lambda r: r[0]):
            w.writerow([pdb, fam[pdb], key])

    multi = sum(1 for k, v in by_key.items() if len(v) > 1)
    print(f"wrote {len(rows)} PDBs to {args.out}")
    print(f"{len(by_key)} distinct receptors; {multi} receptors shared by >1 PDB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
