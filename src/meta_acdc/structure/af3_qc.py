"""Strict per-model AF3 structural QC.

Requires --release, repeatable --root, and a new --out TSV. Confidence must
belong to the same model. Parses/identities/peptide integrity must pass even
when confidence is high. ipTM>=0.75 is a project QC filter, not proof of binding.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

from meta_acdc.structure.cif import parse_cif
from meta_acdc.structure.graph import PEPTIDE, TCR, classify_chains

IPMT_LOW = 0.75       # below this, resubmission instability is observed
TCR_CA_MIN = 400      # complete alpha+beta is ~439-443 CA in this set


def conf_values(inst_dir: Path, job: str) -> tuple[list[float], list[float]]:
    iptm, ptm = [], []
    for f in sorted(inst_dir.glob(f"{job}_conf*.json")):
        try:
            d = json.load(open(f))
        except Exception:  # noqa: BLE001
            continue
        if d.get("iptm") is not None:
            iptm.append(float(d["iptm"]))
        if d.get("ptm") is not None:
            ptm.append(float(d["ptm"]))
    return iptm, ptm


def chain_ca_counts(cif: Path) -> tuple[int, int]:
    residues = parse_cif(cif, model="1")
    roles = classify_chains(residues)
    from collections import Counter
    n = Counter(r.chain for r in residues)
    tcr = sum(c for ch, c in n.items() if roles.get(ch) == TCR)
    pep = sum(c for ch, c in n.items() if roles.get(ch) == PEPTIDE)
    return tcr, pep


def main() -> int:
    ap=argparse.ArgumentParser(description="Per-model QC with source hashes and reference identity")
    ap.add_argument('--root',type=Path,action='append',required=True)
    ap.add_argument('--release',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    from meta_acdc.data.release import load_release
    from meta_acdc.structure.quality import build_qc
    load_release(args.release,require_training=False)
    references=json.loads((args.release.parent/'structure_records.json').read_text())
    rows=build_qc(args.root,references,args.out)
    print(f"per-model QC rows={len(rows)} pass={sum(not x['flag'] for x in rows)}")
    return 0


def _qc_one_dir(inst_dir: Path, inst: str, rows: list[dict]) -> None:
    jobs = sorted({os.path.basename(f).rsplit("_model_0", 1)[0]
                   for f in glob.glob(str(inst_dir / "*_model_0.cif"))})
    for job in jobs:
        iptm, ptm = conf_values(inst_dir, job)
        tcr_ca = pep_ca = -1
        try:
            tcr_ca, pep_ca = chain_ca_counts(inst_dir / f"{job}_model_0.cif")
        except Exception:  # noqa: BLE001
            tcr_ca = pep_ca = -1
        m = float(np.mean(iptm)) if iptm else None
        flags = []
        if tcr_ca < 0 or pep_ca < 0:
            flags.append("parse_failed")
        elif not tcr_ca or not pep_ca:
            flags.append("missing_tcr_or_peptide")
        if not iptm:
            # no confidence basis at all: explicit rejection flag, NOT an
            # empty flag that downstream would read as "passed" (AUDIT.md #26)
            flags.append("no_confidence")
        if m is not None and m < IPMT_LOW:
            flags.append("low_iptm")
        # NB: tcr_ca is reported but NOT flagged — some crystals (e.g.
        # 1ao7, 324 CA) are legitimately truncated constructs, so a low
        # CA count is not by itself an AF3 defect (corrected 2026-09-23).
        rows.append({
            "instance": inst, "job_id": job, "n_conf": len(iptm),
            "iptm_mean": "" if m is None else f"{m:.3f}",
            "iptm_std": "" if len(iptm) < 2 else f"{np.std(iptm):.3f}",
            "ptm_mean": "" if not ptm else f"{np.mean(ptm):.3f}",
            "tcr_ca": tcr_ca, "pep_ca": pep_ca,
            "flag": ";".join(flags),
            "basis": "conf-file" if iptm else "none",
        })


if __name__ == "__main__":
    sys.exit(main())
