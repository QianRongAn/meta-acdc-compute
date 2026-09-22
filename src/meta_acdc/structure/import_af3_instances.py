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
import hashlib
import shutil
import sys
from collections import defaultdict
from pathlib import Path

from meta_acdc.structure.import_af3 import (load_manifest, matches_chain_pair,
                                            read_job_sequences, tcr_pair)

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
    alpha, beta = tcr_pair(seqs, peptide)
    if alpha == CLINICAL_ALPHA:
        return {"pdb": "5brz", "peptide": peptide, "evidence": "clinical"}
    candidates = by_peptide.get(peptide, [])
    match = next((c for c in candidates
                  if matches_chain_pair(alpha, beta, c)), None)
    if match is None and candidates:
        match = candidates[0]
    if match is None:
        return None
    return {"pdb": match["pdb"], "peptide": peptide,
            "evidence": match["evidence"]}


def file_md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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

    # 2. assign each new folder to B1, B2, B3 ... per job_id (chronological).
    # Seed the per-job counter from existing B* instance dirs so successive
    # import runs continue numbering (B4, B5, ...) instead of colliding.
    occ: dict[str, int] = defaultdict(int)
    for inst_dir in sorted(args.out.iterdir()):
        if not inst_dir.is_dir() or not inst_dir.name.startswith("B"):
            continue
        try:
            idx = int(inst_dir.name[1:])
        except ValueError:
            continue
        for f in inst_dir.glob("*_model_0.cif"):
            job = f.name.rsplit("_model_0", 1)[0]
            occ[job] = max(occ[job], idx)
    # content signatures of everything already imported -> idempotent re-runs
    existing_sigs: set[str] = set()
    for d in args.out.iterdir():
        if d.is_dir():
            for c in d.glob("*_model_0.cif"):
                existing_sigs.add(file_md5(c))

    for job_dir in sorted(args.src.iterdir()):
        if not job_dir.is_dir():
            continue
        m0 = sorted(job_dir.glob("*_model_0.cif"))
        if m0 and file_md5(m0[0]) in existing_sigs:
            print(f"  skip {job_dir.name} (already imported)", flush=True)
            continue
        info = resolve_job(job_dir, by_peptide)
        if info is None:
            print(f"  UNMATCHED {job_dir.name}", flush=True)
            continue
        job_id = f"{info['pdb']}_{info['peptide']}"
        occ[job_id] += 1
        instance = f"B{occ[job_id]}"
        n = copy_instance_files(job_dir, args.out / instance, job_id)
        if m0:
            existing_sigs.add(file_md5(m0[0]))
        print(f"  {instance:3s} {job_id:22s} <- {job_dir.name} "
              f"({n} models, {info['evidence']})", flush=True)

    # rebuild the FULL manifest from disk so successive runs stay cumulative
    # (not just this run's additions)
    manifest = []
    for inst_dir in sorted(args.out.iterdir()):
        if not inst_dir.is_dir():
            continue
        inst = inst_dir.name
        for cif in sorted(inst_dir.glob("*_model_0.cif")):
            job_id = cif.name.rsplit("_model_0", 1)[0]
            n = sum(1 for _ in inst_dir.glob(f"{job_id}_model_?.cif"))
            peptide = job_id.split("_", 1)[1]
            cands = by_peptide.get(peptide, [])
            ev = cands[0]["evidence"] if cands else ""
            manifest.append({"instance": inst, "job_id": job_id,
                             "src_folder": "", "n_models": n, "evidence": ev})
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
