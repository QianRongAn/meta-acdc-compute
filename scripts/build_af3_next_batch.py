"""Build the next AF3 web-submission batch from the pending-upload queue.

Protocol (established 2026-10-12 in af3_pending_upload/README.md): AF3 only
sees sequences, so "already submitted" is decided at *content* level. A job is
excluded when any of:

  1. its job_id is in `af3_pending_upload/already_uploaded.tsv`;
  2. its job_id has result files on disk (`af3_predictions/`,
     `af3_native_predictions/`, `af3_instances/{A,B1-B4}/`, quarantine dir);
  3. its 5-chain FASTA md5 equals that of any already-submitted/packaged job
     (same sequences => same AF3 prediction => duplicate submission);
  4. its job_id is an obsolete relabeled ID (`obsolete_relabeled_do_not_reupload.tsv`);
  5. it is already packaged in af3_next60_batch1 / af3_next60 / af3_next60_batch3
     / af3_next60_submitted.

Remaining content-unique jobs are copied into a new batch folder with a
`manifest.tsv` and a `SUBMIT_NOW.tsv` (submission order: grouped by PDB so the
same TCR/MHC chains are checked together). A hand-written README.md is
expected per batch (see af3_next60_batch3/README.md).

Usage:
    python3 scripts/build_af3_next_batch.py \
        --out data/processed/af3_next60_batch4 \
        --limit 60
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import os
import shutil
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PENDING_MANIFEST = "af3_pending_upload/pending_upload_manifest.tsv"
LEDGER = "af3_pending_upload/already_uploaded.tsv"
OBSOLETE = "af3_pending_upload/obsolete_relabeled_do_not_reupload.tsv"
ON_DISK_GLOBS = [
    "data/raw/af3_predictions/*_conf.json",
    "data/raw/af3_native_predictions/*_conf.json",
    "data/raw/af3_instances/*/*_conf.json",
    "data/raw/af3_predictions_quarantine_20261012/*_conf.json",
]
PACKAGED_DIRS = [
    "data/processed/af3_next60_batch1",
    "data/processed/af3_next60",
    "data/processed/af3_next60_batch3",
    "data/processed/af3_next60_submitted",
]


def md5_file(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fasta_index() -> dict[str, list[str]]:
    idx: dict[str, list[str]] = {}
    for path in glob.glob(f"{ROOT}/**/*.fasta", recursive=True):
        if "/.git/" in path:
            continue
        job_id = os.path.basename(path)[: -len(".fasta")]
        idx.setdefault(job_id, []).append(path)
    return idx


def known_jobs() -> tuple[set[str], set[str]]:
    """Return (submit-or-skip job_ids, job_ids with results on disk)."""
    known, disk = set(), set()
    with open(f"{ROOT}/{LEDGER}") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("job_id"):
                known.add(row["job_id"])
    for pattern in ON_DISK_GLOBS:
        for path in glob.glob(f"{ROOT}/{pattern}"):
            job_id = os.path.basename(path)[: -len("_conf.json")]
            disk.add(job_id)
            known.add(job_id)
    with open(f"{ROOT}/{OBSOLETE}") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("old_job_id"):
                known.add(row["old_job_id"])
    return known, disk


def categorized(rows: list[dict]) -> list[dict]:
    """Label a family as '<pep>_matrix' when it spans >=3 jobs, else single."""
    pep_n = Counter(r["peptide"] for r in rows)
    for r in rows:
        r["category"] = (
            f"{r['peptide']}_matrix" if pep_n[r["peptide"]] >= 3 else "single_point"
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="batch folder to create/populate")
    ap.add_argument("--limit", type=int, default=0, help="cap batch size (0 = all)")
    ap.add_argument(
        "--dry-run", action="store_true", help="only report what would ship"
    )
    args = ap.parse_args()

    known, _disk = known_jobs()
    packaged = {
        os.path.basename(p)[: -len(".fasta")]
        for d in PACKAGED_DIRS
        for p in glob.glob(f"{ROOT}/{d}/*.fasta")
    }
    fastas = fasta_index()

    # md5 pool of every job that is submitted, on disk or already packaged
    pool = {}
    for job_id in sorted(known | packaged):
        for path in fastas.get(job_id, []):
            pool.setdefault(md5_file(path), []).append(job_id)

    with open(f"{ROOT}/{PENDING_MANIFEST}") as fh:
        pending = list(csv.DictReader(fh, delimiter="\t"))

    survivors, dupes, no_fasta = [], [], []
    for row in pending:
        job_id = row["job_id"]
        if job_id in known or job_id in packaged:
            continue
        paths = fastas.get(job_id, [])
        if not paths:
            no_fasta.append(job_id)
            continue
        digest = md5_file(paths[0])
        if digest in pool:
            dupes.append((job_id, pool[digest][0]))
            continue
        survivors.append({**row, "source": os.path.relpath(paths[0], ROOT), "md5": digest})

    survivors = categorized(survivors)
    survivors.sort(key=lambda r: (r["pdb"], r["peptide"]))

    print(f"pending manifest jobs      : {len(pending)}")
    print(f"known/on-disk/packaged ids : {len(known | packaged)}")
    print(f"content-duplicate catches  : {len(dupes)} {dupes}")
    if no_fasta:
        print(f"WARNING no fasta found     : {len(no_fasta)} {no_fasta}")
    print(f"REMAINING content-unique   : {len(survivors)}")

    batch = survivors[: args.limit] if args.limit else survivors
    print(f"batch selection            : {len(batch)}")

    by_cat = Counter(r["category"] for r in batch)
    by_pdb = Counter(r["pdb"] for r in batch)
    print(f"categories                 : {dict(by_cat)}")
    print(f"pdbs ({len(by_pdb)})                 : {dict(by_pdb)}")

    if args.dry_run:
        return

    out_dir = f"{ROOT}/{args.out}"
    os.makedirs(out_dir, exist_ok=True)
    for stale in glob.glob(f"{out_dir}/*.fasta"):
        os.remove(stale)
    for r in batch:
        shutil.copyfile(r["source"], f"{out_dir}/{r['job_id']}.fasta")

    with open(f"{out_dir}/manifest.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["job_id", "pdb", "peptide", "category", "source", "md5"])
        for r in batch:
            w.writerow(
                [r["job_id"], r["pdb"], r["peptide"], r["category"], r["source"], r["md5"]]
            )

    with open(f"{out_dir}/SUBMIT_NOW.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["order", "job_id", "pdb", "peptide", "category"])
        for i, r in enumerate(batch, 1):
            w.writerow([i, r["job_id"], r["pdb"], r["peptide"], r["category"]])

    print(f"wrote {len(batch)} fasta + manifest.tsv + SUBMIT_NOW.tsv -> {args.out}/")


if __name__ == "__main__":
    main()
