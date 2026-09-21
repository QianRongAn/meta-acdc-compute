"""Import AF3 web-download results into the pipeline (KN-5).

The AF3 web UI downloads one folder per job with model_0..4 CIFs and
confidence JSONs. This script:
1. Reads each job_request.json to recover the submitted sequences
2. Matches them against kn5_submission_list.tsv (peptide + TCR alpha)
3. Copies model_0 (top-ranked) CIF + summary_confidences_0.json into
   data/raw/af3_predictions/{job_id}_model_0.cif (+ _conf.json)

Usage:
    .venv/bin/python src/meta_acdc/structure/import_af3.py \
        --src "/mnt/c/Users/Administrator/Downloads/folds_2026_09_21_11_38" \
        --list data/processed/kn5_submission_list.tsv \
        --out data/raw/af3_predictions
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path


def read_job_sequences(request_path: Path) -> list[str]:
    with open(request_path) as fh:
        req = json.load(fh)
    if isinstance(req, list):
        req = req[0]  # web downloads wrap the job in a list
    seqs = []
    for item in req.get("sequences", []):
        for kind, spec in item.items():
            if isinstance(spec, dict) and "sequence" in spec:
                seqs.append(spec["sequence"])
    return sorted(seqs, key=len)  # peptide first


def load_manifest(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, required=True)
    ap.add_argument("--list", type=Path,
                    default=Path("data/processed/kn5_submission_list.tsv"))
    ap.add_argument("--out", type=Path, default=Path("data/raw/af3_predictions"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(args.list)
    by_peptide: dict[str, list[dict]] = {}
    for row in manifest:
        by_peptide.setdefault(row["peptide"], []).append(row)

    imported = []
    unmatched = []
    for job_dir in sorted(args.src.iterdir()):
        if not job_dir.is_dir():
            continue
        reqs = sorted(job_dir.glob("*_job_request.json"))
        if not reqs:
            continue
        req = reqs[0]
        seqs = read_job_sequences(req)
        if not seqs:
            continue
        peptide = seqs[0]
        alpha = next((s for s in seqs if s != peptide and ("YFC" in s or "YLC" in s)),
                     None)
        # match manifest: same peptide + same TCR alpha
        candidates = by_peptide.get(peptide, [])
        match = None
        for c in candidates:
            if alpha and c["tcr_a"] and alpha == c["tcr_a"]:
                match = c
                break
        if match is None and candidates:
            match = candidates[0]
        if match is None:
            unmatched.append((job_dir.name, peptide))
            continue
        job_id = f"{match['pdb']}_{peptide}"
        cif_srcs = sorted(job_dir.glob("*_model_0.cif"))
        conf_srcs = sorted(job_dir.glob("*_summary_confidences_0.json"))
        cif_src = cif_srcs[0] if cif_srcs else None
        conf_src = conf_srcs[0] if conf_srcs else None
        cif_dst = args.out / f"{job_id}_model_0.cif"
        conf_dst = args.out / f"{job_id}_conf.json"
        if cif_src:
            shutil.copy2(cif_src, cif_dst)
        if conf_src:
            shutil.copy2(conf_src, conf_dst)
        imported.append((job_dir.name, job_id, peptide, match["evidence"]))

    print(f"imported {len(imported)} jobs:")
    for src, job_id, pep, ev in imported:
        print(f"  {job_id:30s} <- {src}  [{ev}]")
    if unmatched:
        print(f"unmatched {len(unmatched)}: {unmatched}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
