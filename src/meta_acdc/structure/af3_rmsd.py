"""Model-independent structural agreement between AF3 resubmissions.

Complements af3_variance.py: AF3's own confidence (ipTM) is stable across
resubmissions, but confidence is not geometry. This script superposes the
receptor (TCR + MHC + B2M) C-alpha atoms of two AF3 submissions with Kabsch
and reports the peptide C-alpha RMSD (and interface RMSD) between them.

If the peptide geometry is essentially identical across resubmissions while
the EGNN score swings, the instability is unambiguously model-side.

Usage:
    .venv/bin/python src/meta_acdc/structure/af3_rmsd.py \
        --root data/raw/af3_instances \
        --instances A B1 B2 B3 \
        --out data/processed/af3_resubmission_rmsd.tsv
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from meta_acdc.structure.cif import parse_cif
from meta_acdc.structure.graph import (B2M, MHC, PEPTIDE, TCR,
                                       classify_chains)


def chain_ca(residues) -> dict[str, dict[int, np.ndarray]]:
    out: dict[str, dict[int, np.ndarray]] = defaultdict(dict)
    for r in residues:
        ca = next((a for a in r.atoms if a[0] == "CA"), None)
        if ca is not None:
            out[r.chain][r.resid] = np.array(ca[1:], dtype=float)
    return out


def load(path: Path) -> tuple[dict[str, dict[int, np.ndarray]], dict[str, str]]:
    residues = parse_cif(path)
    roles = classify_chains(residues)
    return chain_ca(residues), roles


def role_chains(roles: dict[str, str]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    for ch, role in roles.items():
        out[role].append(ch)
    return out


def matched(coords_a, roles_a, coords_b, roles_b, kinds: set[str]):
    """Matched C-alpha arrays, pairing chains by ROLE (unambiguous for
    mhc/b2m/peptide; TCR alpha/beta are skipped — their chain IDs are not
    guaranteed consistent across AF3 jobs)."""
    ra, rb = role_chains(roles_a), role_chains(roles_b)
    pa, pb = [], []
    for role in kinds:
        ca, cb = ra.get(role, []), rb.get(role, [])
        if len(ca) != 1 or len(cb) != 1:
            continue
        for resid, xyz in coords_a[ca[0]].items():
            if resid in coords_b[cb[0]]:
                pa.append(xyz)
                pb.append(coords_b[cb[0]][resid])
    if len(pa) < 3:
        return None, None
    return np.asarray(pa), np.asarray(pb)


def rmsd_after_superpose(ref, coords, roles, ref_coords, ref_roles,
                         kinds_target, kinds_fit):
    """RMSD over kinds_target after fitting on kinds_fit (ref -> coords)."""
    Pf, Qf = matched(ref_coords, ref_roles, coords, roles, kinds_fit)
    if Pf is None:
        return None
    Pt, Qt = matched(ref_coords, ref_roles, coords, roles, kinds_target)
    if Pt is None:
        return None
    # solve fit on kinds_fit
    Pc = Pf - Pf.mean(0)
    Qc = Qf - Qf.mean(0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = Qf.mean(0) - R @ Pf.mean(0)
    Qt_pred = (R @ Pt.T).T + t
    return float(np.sqrt(np.mean(np.sum((Qt_pred - Qt) ** 2, axis=1))))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/raw/af3_instances"))
    ap.add_argument("--instances", nargs="+", default=None)
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/af3_resubmission_rmsd.tsv"))
    args = ap.parse_args()
    if args.instances is None:
        args.instances = sorted(
            (d.name for d in args.root.iterdir() if d.is_dir()),
            key=lambda s: (s != "A", s))

    receptor = {TCR, MHC, B2M}

    def ca_count(coords, roles, kinds):
        return sum(len(m) for ch, m in coords.items() if roles.get(ch) in kinds)

    jobs: dict[str, dict[str, tuple]] = defaultdict(dict)
    for inst in args.instances:
        d = args.root / inst
        if not d.is_dir():
            continue
        for f in sorted(d.glob("*_model_0.cif")):
            job = os.path.basename(f).rsplit("_model_0", 1)[0]
            try:
                jobs[job][inst] = load(f)
            except Exception as e:  # noqa: BLE001
                print(f"  skip {inst}/{job}: {e}", flush=True)

    rows = []
    for job, per in sorted(jobs.items()):
        insts = [i for i in args.instances if i in per]
        # reference = most complete receptor (a defective first submission
        # with truncated TCR chains must not be used as the frame)
        ref_i = max(insts, key=lambda i: ca_count(per[i][0], per[i][1], receptor))
        ref_coords, ref_roles = per[ref_i]
        entry = {"job_id": job, "ref": ref_i}
        for inst in insts:
            coords, roles = per[inst]
            entry[f"tcr_ca_{inst}"] = str(ca_count(coords, roles, {TCR}))

        def add(col, inst, kinds_target, kinds_fit):
            coords, roles = per[inst]
            v = rmsd_after_superpose(ref_coords, coords, roles, ref_coords,
                                     ref_roles, kinds_target, kinds_fit)
            entry[f"{col}_{inst}"] = "" if v is None else f"{v:.3f}"

        for inst in insts:
            if inst == ref_i:
                entry[f"pep_rmsd_{inst}"] = "0.000"
                entry[f"mhc_rmsd_{inst}"] = "0.000"
                continue
            add("pep_rmsd", inst, {PEPTIDE}, {MHC, B2M})
            add("mhc_rmsd", inst, {MHC}, {MHC, B2M})
        rows.append(entry)

    fields = ["job_id", "ref"] \
        + [f"tcr_ca_{i}" for i in args.instances] \
        + [f"{k}_{i}" for i in args.instances
           for k in ("pep_rmsd", "mhc_rmsd")]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields,
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # summary split: pure resubmission (B-B) vs first-vs-resubmission (A-B)
    def collect(pair_filter):
        vals = []
        for r in rows:
            for inst in args.instances:
                k = f"pep_rmsd_{inst}"
                if not r.get(k, "") or inst == r["ref"]:
                    continue
                if pair_filter(r["ref"], inst):
                    vals.append(float(r[k]))
        return vals

    for label, filt in [("B-B (pure resubmission)",
                         lambda a, b: a.startswith("B") and b.startswith("B")),
                        ("A-B (first vs resubmission)",
                         lambda a, b: "A" in (a, b))]:
        vals = collect(filt)
        if vals:
            print(f"{label}: n={len(vals)} peptide CA RMSD median "
                  f"{np.median(vals):.3f} A, max {max(vals):.3f} A")
    return 0


if __name__ == "__main__":
    sys.exit(main())
