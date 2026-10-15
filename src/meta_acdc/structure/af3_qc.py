"""Structural QC / reliability gate for AF3 predictions.

Derived from the 2026-09-23 resubmission analysis (docs/benchmarks.md):
- high-confidence predictions (ipTM >= ~0.85) reproduce across independent
  AF3 runs (ipTM Delta <= 0.07, peptide register RMSD <= 0.43 A);
- low-confidence complexes (ipTM ~0.5, e.g. JM22) genuinely vary across runs;
- some first submissions were defective (truncated TCR chains, 324 vs 443 CA).

This script flags, per (instance, job) model set, the metrics that predict
whether a structure is trustworthy — an ipTM gate plus a chain-completeness
check — so downstream scoring/ranking can exclude or annotate bad inputs.

Usage:
    .venv/bin/python src/meta_acdc/structure/af3_qc.py \
        --root data/raw/af3_instances \
        --out data/processed/af3_qc.tsv
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, action="append", default=[],
                    help="structure dir to QC (repeatable); defaults to "
                         "data/raw/af3_instances and data/raw/af3_predictions "
                         "so the QC table covers both the instance set and "
                         "the scoring panel")
    ap.add_argument("--out", type=Path, default=Path("data/processed/af3_qc.tsv"))
    args = ap.parse_args()
    roots = args.root or [Path("data/raw/af3_instances"),
                          Path("data/raw/af3_predictions")]

    rows = []
    for root in roots:
        if not root.exists():
            print(f"QC root not found: {root}", file=sys.stderr)
            continue
        if any(root.glob("*_model_0.cif")):
            _qc_one_dir(root, root.name, rows)
        else:
            for inst_dir in sorted(root.iterdir()):
                if not inst_dir.is_dir():
                    continue
                _qc_one_dir(inst_dir, inst_dir.name, rows)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "instance", "job_id", "n_conf", "iptm_mean", "iptm_std",
            "ptm_mean", "tcr_ca", "pep_ca", "flag", "basis"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    n = len(rows)
    bad = [r for r in rows if r["flag"]]
    n_low = sum(1 for r in rows if "low_iptm" in r["flag"])
    n_nc = sum(1 for r in rows if "no_confidence" in r["flag"])
    print(f"wrote {n} structure records to {args.out}")
    print(f"  flagged: {len(bad)} (low_iptm={n_low}, no_confidence={n_nc}); "
          f"gate: ipTM<{IPMT_LOW}")
    for r in bad:
        print(f"  {r['instance']:3s} {r['job_id']:22s} "
              f"ipTM={r['iptm_mean']} tcr_ca={r['tcr_ca']} [{r['flag']}]")
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
            pass
        m = float(np.mean(iptm)) if iptm else None
        flags = []
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
