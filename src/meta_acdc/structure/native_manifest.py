"""Build the native AF3 prediction set + manifest from submitted folds.

Reproduces the 2026-09-23 native domain-adaptation intake. The submitted AF3
jobs are matched to the *native submission FASTAs* (`af3_native*/`) by exact
chain-set fingerprint; only genuine cognate natives are copied into the
native prediction directory, and a manifest of the accepted job_ids is
written. Cross-reactivity candidate jobs are excluded from positives.

Usage:
    .venv/bin/python src/meta_acdc/structure/native_manifest.py \
        --src /mnt/c/.../folds_2026_09_23_01_55 ... \
        --out data/raw/af3_native_predictions \
        --manifest data/processed/af3_native/manifest.tsv \
        --report data/processed/af3_native/import_report.tsv
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from meta_acdc.structure.import_af3 import read_job_sequences  # noqa: E402


def read_fasta_seqs(path: Path) -> list[str]:
    return [ln.strip() for ln in open(path)
            if not ln.startswith(">") and ln.strip()]


def build_native_index(fasta_dirs: list[Path]) -> tuple[dict[frozenset, str],
                                                        set[frozenset]]:
    """Map chain-set fingerprint -> native job_id; track ambiguous sets.

    Some native complexes share an identical chain set (e.g. the same peptide
    solved with near-identical TCRs). For those the fingerprint is ambiguous:
    we keep the first id but record the set so callers can report it.
    """
    idx: dict[frozenset, str] = {}
    ambiguous: set[frozenset] = set()
    for d in fasta_dirs:
        for f in sorted(glob.glob(str(d / "*.fasta"))):
            job = os.path.basename(f)[:-6]
            key = frozenset(read_fasta_seqs(f))
            if key in idx and idx[key] != job:
                ambiguous.add(key)
            idx.setdefault(key, job)
    return idx, ambiguous


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, action="append", required=True,
                    help="folds_* download dirs (repeatable)")
    ap.add_argument("--fasta-dirs", type=Path, action="append", default=None,
                    help="native FASTA dirs (default: data/processed/af3_native*)")
    ap.add_argument("--out", type=Path,
                    default=Path("data/raw/af3_native_predictions"))
    ap.add_argument("--manifest", type=Path,
                    default=Path("data/processed/af3_native/manifest.tsv"))
    ap.add_argument("--report", type=Path,
                    default=Path("data/processed/af3_native/import_report.tsv"))
    args = ap.parse_args()

    fasta_dirs = args.fasta_dirs or sorted(
        Path("data/processed").glob("af3_native*"))
    fasta_dirs = [d for d in fasta_dirs if d.is_dir()]
    idx, ambiguous = build_native_index(fasta_dirs)
    print(f"native FASTA index: {len(idx)} complexes from "
          f"{[d.name for d in fasta_dirs]}", flush=True)
    if ambiguous:
        print(f"note: {len(ambiguous)} chain-set fingerprints are ambiguous "
              f"(shared by >1 native id); first id kept", flush=True)

    args.out.mkdir(parents=True, exist_ok=True)
    # CUMMULATIVE: union the accepted jobs with the existing manifest so a
    # re-run over a single old folds dir can never shrink the native set
    # (the 2026-09-23 cron clobber: orchestrate_da.sh reprocessed old dirs
    # and native_manifest overwrote the 89-native manifest with 1 job).
    existing: dict[str, tuple[str, str, str, str]] = {}
    if args.manifest.exists():
        with open(args.manifest, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                existing[r["job_id"]] = (r["pdb"], r["peptide"],
                                         r.get("cdr3", ""), r.get("mhc_hint", ""))

    seen: set[str] = set()
    report = []
    unmatched = 0
    for base in args.src:
        for job_dir in sorted(glob.glob(str(base) + "/*/")):
            reqs = sorted(glob.glob(job_dir + "*_job_request.json"))
            if not reqs:
                continue
            seqs = read_job_sequences(Path(reqs[0]))
            job = idx.get(frozenset(seqs)) if seqs else None
            if job is None:
                unmatched += 1
                continue
            if job in seen:
                continue
            seen.add(job)
            for cif in sorted(glob.glob(job_dir + "*_model_?.cif")):
                tag = os.path.basename(cif).split("_model_")[1].split(".")[0]
                shutil.copy2(cif, args.out / f"{job}_model_{tag}.cif")
            cf = glob.glob(job_dir + "*_summary_confidences_0.json")
            if cf:
                shutil.copy2(cf[0], args.out / f"{job}_conf.json")
            report.append((job, os.path.basename(job_dir.rstrip("/")),
                           seqs[0] if seqs else ""))

    # every native CIF currently on disk belongs in the manifest
    on_disk = {p.name.rsplit("_model_0", 1)[0]
               for p in args.out.glob("*_model_0.cif")}
    all_ids = sorted(set(existing) | {r[0] for r in report} | on_disk)
    with open(args.manifest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["job_id", "pdb", "peptide", "cdr3", "mhc_hint"])
        for job in all_ids:
            pdb, pep, cdr3, mhc = existing.get(
                job, (job.split("_", 1)[0], job.split("_", 1)[1], "", ""))
            w.writerow([job, pdb, pep, cdr3, mhc])
    with open(args.report, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["job_id", "source", "peptide"])
        w.writerows(report)
    print(f"matched {len(seen)} in this run; manifest now {len(all_ids)} "
          f"natives ({len(on_disk)} CIFs on disk); excluded {unmatched}")
    print(f"-> {args.out}, {args.manifest}, {args.report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())