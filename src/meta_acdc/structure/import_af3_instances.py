"""Replicate-aware import of AF3 web results (KN-5 resubmission variance).

Unlike import_af3.py (which overwrites {job_id}_model_0..4 on every import),
this script keeps every independent AF3 submission as its own *instance*
directory so that resubmission variance can be quantified:

    data/raw/af3_instances/<instance>/{job_id}_model_{0..4}.cif (+ _conf.json)

- Existing already-imported predictions (data/raw/af3_predictions) are copied
  in as instance "A".
- Each new AF3 web job folder is assigned to B1, B2, B3 ... = the 1st/2nd/3rd
  *independent* submission seen for that job_id inside --src.

Usage:
    .venv/bin/python src/meta_acdc/structure/import_af3_instances.py \
        --src "/mnt/c/Users/Administrator/Documents/folds_2026_09_22_18_41" \
        --list data/processed/kn5_submission_list.tsv \
        --existing data/raw/af3_predictions \
        --out data/raw/af3_instances \
        --manifest data/processed/af3_instance_manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from collections import defaultdict
from pathlib import Path

from meta_acdc.structure.import_af3 import load_manifest, read_job_sequences

CLINICAL_ALPHA = (
    "AQEVTQIPAALSVPEGENLVLNCSFTDSAIYNLQWFRQDPGKGLTSLLYVRPYQREQTSGRLNASLDKS"
    "SGRSTLYIAASQPGDSATYLCAVRPGGAGPFFVVFGKGTKLSVIPNIQNPDPAVYQLRDSKSSDKSVCL"
    "FTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIP"
)


def resolve_job(job_dir: Path, by_peptide: dict[str, list[dict]]) -> dict | None:
    reqs = sorted(job_dir.glob("*_job_request.json"))
    if not reqs:
        return None
    seqs = read_job_sequences(reqs[0])
    if not seqs:
        return None
    peptide = seqs[0]
    alpha = next((s for s in seqs if s != peptide and ("YFC" in s or "YLC" in s)),
                 None)
    if alpha == CLINICAL_ALPHA:
        return {"pdb": "5brz", "peptide": peptide, "evidence": "clinical"}
    candidates = by_peptide.get(peptide, [])
    match = None
    for c in candidates:
        if alpha and c["tcr_a"] and alpha == c["tcr_a"]:
            match = c
            break
    if match is None and candidates:
        match = candidates[0]
    if match is None:
        return None
    return {"pdb": match["pdb"], "peptide": peptide,
            "evidence": match["evidence"]}


def copy_instance_files(src: Path, dest: Path, job_id: str) -> int:
    dest.mkdir(parents=True, exist_ok=True)
    n = 0
    for cif in sorted(src.glob("*_model_?.cif")):
        tag = cif.name.split("_model_")[1].split(".")[0]
        shutil.copy2(cif, dest / f"{job_id}_model_{tag}.cif")
        n += 1
    confs = sorted(src.glob("*_summary_confidences_?.json"))
    if confs:
        for conf in confs:
            tag = conf.name.split("_summary_confidences_")[1].split(".")[0]
            shutil.copy2(conf, dest / f"{job_id}_conf_{tag}.json")
    else:
        conf0 = sorted(src.glob("*_conf*.json"))
        if conf0:
            shutil.copy2(conf0[0], dest / f"{job_id}_conf.json")
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, required=True)
    ap.add_argument("--list", type=Path,
                    default=Path("data/processed/kn5_submission_list.tsv"))
    ap.add_argument("--existing", type=Path,
                    default=Path("data/raw/af3_predictions"))
    ap.add_argument("--out", type=Path, default=Path("data/raw/af3_instances"))
    ap.add_argument("--manifest", type=Path,
                    default=Path("data/processed/af3_instance_manifest.tsv"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    by_peptide: dict[str, list[dict]] = {}
    for row in load_manifest(args.list):
        by_peptide.setdefault(row["peptide"], []).append(row)

    manifest: list[dict] = []

    # 1. seed instance A from the previously-imported predictions
    if args.existing.is_dir():
        dest = args.out / "A"
        seen_a: set[str] = set()
        for conf in sorted(args.existing.glob("*_conf*.json")):
            job_id = conf.name.split("_conf")[0]
            if job_id in seen_a:
                continue
            seen_a.add(job_id)
            n = sum(1 for _ in args.existing.glob(f"{job_id}_model_?.cif"))
            for cif in sorted(args.existing.glob(f"{job_id}_model_?.cif")):
                tag = cif.name.split("_model_")[1].split(".")[0]
                dest.mkdir(parents=True, exist_ok=True)
                shutil.copy2(cif, dest / f"{job_id}_model_{tag}.cif")
            shutil.copy2(conf, dest / f"{job_id}_conf.json")
            manifest.append({"instance": "A", "job_id": job_id,
                             "src_folder": args.existing.name,
                             "n_models": n, "evidence": ""})

    # 2. assign each new folder to B1, B2, B3 ... per job_id (chronological)
    occ: dict[str, int] = defaultdict(int)
    for job_dir in sorted(args.src.iterdir()):
        if not job_dir.is_dir():
            continue
        info = resolve_job(job_dir, by_peptide)
        if info is None:
            print(f"  UNMATCHED {job_dir.name}", flush=True)
            continue
        job_id = f"{info['pdb']}_{info['peptide']}"
        occ[job_id] += 1
        instance = f"B{occ[job_id]}"
        n = copy_instance_files(job_dir, args.out / instance, job_id)
        manifest.append({"instance": instance, "job_id": job_id,
                         "src_folder": job_dir.name, "n_models": n,
                         "evidence": info["evidence"]})
        print(f"  {instance:3s} {job_id:22s} <- {job_dir.name} "
              f"({n} models, {info['evidence']})", flush=True)

    with open(args.manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=["instance", "job_id",
                         "src_folder", "n_models", "evidence"])
        w.writeheader()
        for m in manifest:
            w.writerow(m)
    print(f"wrote {len(manifest)} entries to {args.manifest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
