"""Prepare AlphaFold3 web-submission batches from the KN-5 candidate list.

For each cross-reactivity pair (pdb, vdjdb_epitope), build the full complex
input: TCR alpha/beta + peptide + MHC heavy chain + B2M — all sequences taken
from the original PDB (the peptide is the only substituted chain). Output is a
multi-chain FASTA per job, ready for copy-paste into the AF3 web UI
(alphafoldserver.com) or TCRmodel2.

Usage:
    .venv/bin/python src/meta_acdc/structure/af3_batch.py \
        --list data/processed/kn5_submission_list.tsv \
        --struct-dir data/raw/structures \
        --n 50 --out data/processed/af3_batch1
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from meta_acdc.structure.graph import classify_chains, parse_pdb, TCR, MHC, B2M, PEPTIDE
from meta_acdc.structure.vdjdb_structure_map import chain_sequences, aa3to1


def build_job(pdb: str, struct_dir: Path) -> dict[str, str] | None:
    path = struct_dir / f"{pdb.lower()}.pdb"
    if not path.exists():
        return None
    try:
        residues = parse_pdb(path)
        roles = classify_chains(residues)
        seqs = {c: aa3to1(s) for c, s in chain_sequences(path).items()}
    except Exception:
        return None
    tcr_seqs = [s for c, s in seqs.items() if roles.get(c) == TCR]
    mhc = sorted((s for c, s in seqs.items() if roles.get(c) == MHC),
                 key=len, reverse=True)
    b2m = [s for c, s in seqs.items() if roles.get(c) == B2M]
    if len(tcr_seqs) < 2 or not mhc or not b2m:
        return None
    # alpha/beta by conserved pre-CDR3 motif: YFCAS/YLCAS -> beta,
    # YFCAV/YFCAI/YFCAL -> alpha
    alpha = beta = None
    for s in tcr_seqs:
        if any(m in s for m in ("YFCAS", "YLCAS")):
            beta = s
        elif any(m in s for m in ("YFCAV", "YFCAI", "YFCAL", "YLCAV")):
            alpha = s
    if alpha is None and beta is None:
        return None
    if alpha is None:
        others = [s for s in tcr_seqs if s != beta]
        if not others:
            return None
        alpha = others[0]
    if beta is None:
        others = [s for s in tcr_seqs if s != alpha]
        if not others:
            return None
        beta = others[0]
    return {"alpha": alpha, "beta": beta, "mhc": mhc[0], "b2m": b2m[0]}


def fasta_block(job_id: str, alpha: str, beta: str, peptide: str,
                mhc: str, b2m: str) -> str:
    # NOTE: AF3 web UI rejects non-letter characters in headers (|, digits);
    # use plain chain labels only. The job identity lives in the filename.
    lines = [
        ">alpha",
        alpha,
        ">beta",
        beta,
        ">peptide",
        peptide,
        ">mhc",
        mhc,
        ">b2m",
        b2m,
    ]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", type=Path,
                    default=Path("data/processed/kn5_submission_list.tsv"))
    ap.add_argument("--struct-dir", type=Path, default=Path("data/raw/structures"))
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--evidence", default="different",
                    choices=["different", "same", "any"],
                    help="which VDJdb evidence class to submit: 'different' "
                         "= cross-reactivity candidates (default), 'same' = "
                         "native complexes (domain-adaptation positives)")
    ap.add_argument("--out", type=Path, default=Path("data/processed/af3_batch1"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    rows = []
    with open(args.list, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if args.evidence == "any" or row["evidence"] == args.evidence:
                rows.append(row)

    manifest = []
    n_written = 0
    seen_jobs: set[str] = set()
    for row in rows:
        if n_written >= args.n:
            break
        job_id = f"{row['pdb']}_{row['peptide']}"
        if job_id in seen_jobs:  # same (pdb, peptide) hit by multiple CDR3s
            continue
        job = build_job(row["pdb"], args.struct_dir)
        if job is None:
            continue
        seen_jobs.add(job_id)
        block = fasta_block(job_id, job["alpha"], job["beta"],
                            row["peptide"], job["mhc"], job["b2m"])
        (args.out / f"{job_id}.fasta").write_text(block)
        manifest.append({"job_id": job_id, "pdb": row["pdb"],
                         "peptide": row["peptide"], "cdr3": row["cdr3"],
                         "mhc_hint": row["mhc_hint"]})
        n_written += 1

    with open(args.out / "manifest.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["job_id", "pdb", "peptide", "cdr3",
                                           "mhc_hint"], delimiter="\t")
        w.writeheader()
        for m in manifest:
            w.writerow(m)
    print(f"wrote {n_written} FASTA jobs to {args.out}/ (manifest.tsv included)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
