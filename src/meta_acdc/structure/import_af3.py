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


def tcr_pair(seqs: list[str], peptide: str) -> tuple[str | None, str | None]:
    """Alpha/beta by conserved pre-CDR3 motif (order-independent)."""
    rest = [s for s in seqs if s != peptide]
    beta = next((s for s in rest if "YFCAS" in s or "YLCAS" in s), None)
    alpha = next((s for s in rest
                  if any(m in s for m in ("YFCAV", "YFCAI", "YFCAL", "YLCAV"))),
                 None)
    return alpha, beta


def matches_chain_pair(alpha, beta, cand: dict) -> bool:
    """True iff the job's {alpha,beta} equals the candidate's {tcr_a,tcr_b}.

    kn5_submission_list.tsv stores tcr_a = beta and tcr_b = alpha (columns
    are swapped relative to the AF3 job FASTA), so an order-dependent
    alpha == tcr_a test never fires — hence the set comparison."""
    if not alpha or not beta:
        return False
    return {alpha, beta} == {cand.get("tcr_a"), cand.get("tcr_b")}


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
        alpha, beta = tcr_pair(seqs, peptide)
        # clinical gold-standard set (MAG-IC3/5brz TCR): detect FIRST by the
        # exact alpha-chain sequence, before manifest matching
        CLINICAL_ALPHA = (
            "AQEVTQIPAALSVPEGENLVLNCSFTDSAIYNLQWFRQDPGKGLTSLLYVRPYQREQTSGRLNASLDKS"
            "SGRSTLYIAASQPGDSATYLCAVRPGGAGPFFVVFGKGTKLSVIPNIQNPDPAVYQLRDSKSSDKSVCL"
            "FTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIP")
        if alpha == CLINICAL_ALPHA:
            match = {"pdb": "5brz", "evidence": "clinical"}
        else:
            # match manifest: same peptide + same {alpha,beta} chain pair
            candidates = by_peptide.get(peptide, [])
            match = next((c for c in candidates
                          if matches_chain_pair(alpha, beta, c)), None)
            if match is None and candidates:
                match = candidates[0]
            if match is None:
                unmatched.append((job_dir.name, peptide))
                continue
        job_id = f"{match['pdb']}_{peptide}"
        # copy ALL five ranked models (model_0..4) for structural ensembling
        n_copied = 0
        for cif_src in sorted(job_dir.glob("*_model_?.cif")):
            model_tag = cif_src.name.split("_model_")[1].split(".")[0]
            shutil.copy2(cif_src, args.out / f"{job_id}_model_{model_tag}.cif")
            n_copied += 1
        conf_srcs = sorted(job_dir.glob("*_summary_confidences_0.json"))
        if conf_srcs:
            shutil.copy2(conf_srcs[0], args.out / f"{job_id}_conf.json")
        imported.append((job_dir.name, job_id, peptide,
                         match["evidence"], n_copied))

    print(f"imported {len(imported)} jobs:")
    for src, job_id, pep, ev, n in imported:
        print(f"  {job_id:30s} <- {src}  [{ev}] {n} models")
    if unmatched:
        print(f"unmatched {len(unmatched)}: {unmatched}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
