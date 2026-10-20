"""Synthetic structural controls; their class labels are not experimental binding labels.

Grafting aligns donor N/CA/C and transfers its side chain without relaxation.
This can create clashes and local geometry artifacts. Success on this task
only establishes discrimination of this construction protocol. Use the
immutable release builder rather than this legacy exploratory dataset CLI.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from pathlib import Path

from meta_acdc.structure.graph import (
    Residue,
    build_interface_graph_from_residues,
    classify_chains,
    parse_pdb,
    AA3TO1,
    TCR,
    MHC,
    PEPTIDE,
)

PEPTIDE_LEN = (5, 20)


def _split_sides(residues: list[Residue]) -> tuple[list[Residue], list[Residue]]:
    roles = classify_chains(residues)
    tcr = [r for r in residues if roles[r.chain] == TCR]
    pmhc = [r for r in residues if roles[r.chain] != TCR]
    return tcr, pmhc


def _peptide_chains(residues: list[Residue]) -> dict[str, list[Residue]]:
    """Group residues by peptide-role chain."""
    roles = classify_chains(residues)
    chains: dict[str, list[Residue]] = {}
    for r in residues:
        if roles[r.chain] == PEPTIDE:
            chains.setdefault(r.chain, []).append(r)
    return chains


import numpy as np


def _kabsch(src: np.ndarray, dst: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Optimal rotation+translation mapping src onto dst (both (n,3))."""
    c_src = src.mean(axis=0)
    c_dst = dst.mean(axis=0)
    h = (src - c_src).T @ (dst - c_dst)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    t = c_dst - r @ c_src
    return r, t


def _backbone_coords(res: Residue) -> np.ndarray | None:
    """N, CA, C coordinates as (3,3); None if any is missing."""
    m = {a[0]: np.array([a[1], a[2], a[3]]) for a in res.atoms}
    if all(k in m for k in ("N", "CA", "C")):
        return np.vstack([m["N"], m["CA"], m["C"]])
    return None


def _transplanted_atoms(native: Residue, foreign: Residue) -> list[tuple]:
    """Foreign side-chain atoms mapped into the native residue's backbone frame.

    Returns native backbone atoms (kept verbatim) + transformed foreign
    side-chain atoms — the residue's geometry is untouched at C-alpha level,
    only the side-chain chemistry changes.
    """
    src = _backbone_coords(foreign)
    dst = _backbone_coords(native)
    if src is None or dst is None:
        raise ValueError("graft requires complete N/CA/C backbone in both residues")
    r, t = _kabsch(src, dst)
    out = [a for a in native.atoms if a[0] in ("N", "CA", "C", "O", "OXT")]
    for name, x, y, z in foreign.atoms:
        if name in ("N", "CA", "C", "O", "OXT"):
            continue
        p = r @ np.array([x, y, z]) + t
        out.append((name, float(p[0]), float(p[1]), float(p[2])))
    return out


def make_decoy(
    residues: list[Residue],
    peptide: list[Residue],
) -> list[Residue]:
    """Graft a foreign peptide chain in place of the native one.

    Synthetic construction control (not an experimentally validated nonbinder): the foreign peptide's side-chain chemistry
    is transplanted onto the NATIVE peptide's backbone via Kabsch alignment.
    Backbone coordinates, node set and chain labels are untouched; only the
    side-chain atoms (and residue identities) change — the C-alpha-invisible
    signal. This does not establish realistic conformations or nonbinding.
    Peptide lengths must match; no prefix swaps are allowed.
    """
    roles = classify_chains(residues)
    native_pep_chains = [c for c, r in roles.items() if r == PEPTIDE]
    if len(native_pep_chains) != 1:
        raise ValueError("graft requires exactly one peptide chain")
    native_pep_chain = native_pep_chains[0]
    native_pep = [r for r in residues if r.chain == native_pep_chain]
    if len(native_pep) != len(peptide):
        raise ValueError("graft peptide lengths must match; prefix swaps forbidden")
    swap = {id(r): f for r, f in zip(native_pep, peptide)}
    out = []
    for r in residues:
        f = swap.get(id(r))
        if f is not None:
            atoms = _transplanted_atoms(r, f)
            from dataclasses import replace
            out.append(replace(r, resname=f.resname, atoms=atoms))
        else:
            out.append(r)
    return out


def _centroid(rs: list[Residue], atom_names: tuple[str, ...] = ("CA",)) -> np.ndarray:
    """Centroid of the named atoms over a residue list (np.array of 3)."""
    pts = [np.array(a[1:]) for r in rs for a in r.atoms if a[0] in atom_names]
    if not pts:
        raise ValueError("centroid requested over residues without the named atoms")
    return np.mean(np.vstack(pts), axis=0)


def make_displaced_peptide_decoy(
    residues: list[Residue],
    shift: float = 6.0,
) -> list[Residue]:
    """Synthetic direction-aware translation control, not an experimental nonbinder.

    The peptide is translated along the peptide-centroid -> MHC-groove
    direction (i.e. away from the groove in the deposited geometry), not
    along a fixed global axis. The direction is data-dependent, so a
    displacement that would move the peptide sideways (or toward the TCR)
    for one complex is not silently produced for another.

    Still a construction rule, not physics: a 6 A shift can introduce
    peptide/environment clashes rather than remove the peptide from the
    groove, and neither overlap-free geometry nor physical plausibility
    follows. Node membership is pinned by the accepted builder. Apply the
    same artifact screen as to native and graft controls and record
    exclusions.
    """
    if not math.isfinite(shift) or shift <= 0:
        raise ValueError("displacement shift must be finite and positive")
    roles = classify_chains(residues)
    pep = [r for r in residues if roles[r.chain] == PEPTIDE]
    mhc = [r for r in residues if roles[r.chain] == MHC]
    if not pep:
        raise ValueError("displacement requires peptide-role residues")
    if not mhc:
        raise ValueError(
            "direction-aware displacement requires MHC-role residues "
            "(peptide->groove direction is undefined without them)")
    p_c = _centroid(pep)
    g_c = _centroid(mhc)
    d = p_c - g_c
    norm = float(np.linalg.norm(d))
    if norm < 1e-6:
        raise ValueError("degenerate peptide/groove centroids; direction undefined")
    u = d / norm
    out = []
    for r in residues:
        if roles[r.chain] == PEPTIDE:
            shifted_atoms = [(n, x + shift * u[0], y + shift * u[1],
                              z + shift * u[2]) for n, x, y, z in r.atoms]
            from dataclasses import replace
            out.append(replace(r, x=r.x + shift * u[0], y=r.y + shift * u[1],
                               z=r.z + shift * u[2], atoms=shifted_atoms))
        else:
            out.append(r)
    return out


def build_dataset(*args,**kwargs):
    """Historical builder retired; helpers above remain for incident reconstruction."""
    raise ValueError('legacy structural dataset builder retired; use meta_acdc.data.build_release --previous-release-dir')


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct-dir", type=Path, default=Path("data/raw/structures"))
    ap.add_argument("--out", type=Path, default=Path("data/processed/egnn_dataset.pt"))
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--with-electrostatics", action="store_true",
                    help="append v11 atom-level electrostatics/desolvation "
                         "features (node 26->28, edge 16->18)")
    args = ap.parse_args()
    stats = build_dataset(args.struct_dir, args.out, args.seed,
                          with_electrostatics=args.with_electrostatics)
    print(f"dataset saved to {args.out}: {stats}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
