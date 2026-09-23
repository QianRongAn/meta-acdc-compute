"""Build AF3 web-submission jobs for arbitrary peptides against one TCR.

Complements `af3_batch.py` (which uses the VDJdb epitope per structure).
Here the TCR + MHC + B2M come from a chosen crystal structure and the peptide
chain is substituted with each peptide from a list — the input needed to
re-rank a proteome-scale candidate pool with the structural model
(KN-8 / KN-11 structural scan). Output is one 5-chain FASTA per peptide,
ready for copy-paste into the AF3 web UI.

Usage:
    .venv/bin/python src/meta_acdc/structure/af3_scan_batch.py \
        --pdb 5brz \
        --peptides data/processed/kn8_top50k_CASSLGRY.tsv \
        --n 200 --out data/processed/af3_scan_5brz
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from meta_acdc.structure.af3_batch import build_job, fasta_block


def read_peptides(path: Path) -> list[str]:
    """Accept a plain one-per-line list or a TSV with a 'peptide' column."""
    peps: list[str] = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
        fh.seek(0)
        if "\t" in first and "peptide" in first:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("peptide"):
                    peps.append(row["peptide"].strip())
        else:
            for ln in fh:
                ln = ln.strip()
                if ln and not ln.startswith("#"):
                    peps.append(ln)
    seen: set[str] = set()
    return [p for p in peps if not (p in seen or seen.add(p))]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True,
                    help="TCR source structure id (data/raw/structures/<pdb>.pdb)")
    ap.add_argument("--peptides", type=Path, required=True,
                    help="peptide list (one per line, or TSV with 'peptide')")
    ap.add_argument("--struct-dir", type=Path,
                    default=Path("data/raw/structures"))
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    job = build_job(args.pdb, args.struct_dir)
    if job is None:
        print(f"cannot build TCR/MHC/B2M from {args.pdb}", file=sys.stderr)
        return 1
    args.out.mkdir(parents=True, exist_ok=True)

    manifest = []
    for pep in read_peptides(args.peptides):
        if len(manifest) >= args.n:
            break
        if not pep.isalpha():
            continue
        job_id = f"{args.pdb}_{pep}"
        block = fasta_block(job_id, job["alpha"], job["beta"], pep,
                            job["mhc"], job["b2m"])
        (args.out / f"{job_id}.fasta").write_text(block)
        manifest.append({"job_id": job_id, "pdb": args.pdb, "peptide": pep})
    with open(args.out / "manifest.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["job_id", "pdb", "peptide"])
        w.writeheader()
        w.writerows(manifest)
    print(f"wrote {len(manifest)} scan jobs for {args.pdb} to {args.out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())