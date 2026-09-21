"""EGNN dataset: interface graphs + graft-decoy negatives (KN-4 prep).

Positive examples: real TCR-pMHC complex interface graphs (from STCRDab /
RCSB, 289 complexes curated by fetch_structures.py).

Negative examples ("graft decoys"): the peptide chain of complex X is
replaced by the peptide of complex Y (different antigen), then the interface
graph is rebuilt with the same 10A rule. The TCR and MHC geometry stay
real; only the antigen changes — the model must learn physical compatibility
of the interface, not sequence shortcuts. This mirrors the proposal's
hard-negative strategy at the structural level.

Usage:
    .venv/bin/python src/meta_acdc/models/dataset.py \
        --struct-dir data/raw/structures --out data/processed/egnn_dataset.pt
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

from meta_acdc.structure.graph import (
    Residue,
    build_interface_graph_from_residues,
    classify_chains,
    parse_pdb,
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
        return list(native.atoms)  # fallback: keep native atoms
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

    TRUE structural-mimicry decoy: the foreign peptide's side-chain chemistry
    is transplanted onto the NATIVE peptide's backbone via Kabsch alignment.
    Backbone coordinates, node set and chain labels are untouched; only the
    side-chain atoms (and residue identities) change — the C-alpha-invisible
    signal that real cross-reactivity (e.g., titin vs MAGE-A3) exhibits.
    (Peptide lengths may differ: swap the aligned prefix.)
    """
    roles = classify_chains(residues)
    native_pep_chains = [c for c, r in roles.items() if r == PEPTIDE]
    if not native_pep_chains:
        return residues
    native_pep_chain = native_pep_chains[0]
    native_pep = [r for r in residues if r.chain == native_pep_chain]
    swap = {id(r): f for r, f in zip(native_pep, peptide)}
    out = []
    for r in residues:
        f = swap.get(id(r))
        if f is not None:
            atoms = _transplanted_atoms(r, f)
            out.append(Residue(r.chain, f.resname, r.resid, r.x, r.y, r.z,
                               r.plddt, r.icode, atoms))
        else:
            out.append(r)
    return out


def make_displaced_peptide_decoy(
    residues: list[Residue],
    shift: float = 6.0,
) -> list[Residue]:
    """Negative decoy: the native peptide chain rigidly translated out of the
    groove (along +z). Shift is 6A: the peptide leaves the binding groove but
    STAYS WITHIN the 10A interface selection, so peptide nodes remain in the
    graph — the is-peptide flag does not leak the label. Geometry is disrupted,
    not removed. (The graft decoy is the weaker-signal variant.)
    """
    roles = classify_chains(residues)
    out = []
    for r in residues:
        if roles[r.chain] == PEPTIDE:
            shifted_atoms = [(n, x, y, z + shift) for n, x, y, z in r.atoms]
            out.append(Residue(r.chain, r.resname, r.resid, r.x, r.y, r.z + shift,
                               r.plddt, r.icode, shifted_atoms))
        else:
            out.append(r)
    return out


def build_dataset(struct_dir: Path, out_path: Path, seed: int = 42) -> dict:
    rng = random.Random(seed)
    with open(struct_dir / "complexes.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    complexes: dict[str, list[Residue]] = {}
    skipped = 0
    n_split = 0
    for row in rows:
        pdb = row["pdb"]
        path = struct_dir / f"{pdb.lower()}.pdb"
        try:
            residues = parse_pdb(path)
            roles = classify_chains(residues)
        except Exception:
            continue
        n_tcr = sum(1 for v in roles.values() if v == TCR)
        n_peptide = sum(1 for v in roles.values() if v == PEPTIDE)
        n_mhc = sum(1 for v in roles.values() if v == MHC)
        if n_tcr >= 1 and n_peptide == 1 and n_mhc >= 1:
            complexes[pdb] = residues
        elif (n_tcr > 2 or n_peptide > 1) and n_mhc >= 1:
            # multi-complex asymmetric unit: split into sub-complexes
            try:
                from meta_acdc.structure.split_complexes import split_complex
                subs = split_complex(residues)
            except Exception:
                subs = []
            for k, sub in enumerate(subs):
                sub_roles = classify_chains(sub)
                if (sum(1 for v in sub_roles.values() if v == TCR) >= 1
                        and sum(1 for v in sub_roles.values() if v == PEPTIDE) == 1
                        and sum(1 for v in sub_roles.values() if v == MHC) >= 1):
                    complexes[f"{pdb}_{k}"] = sub
                    n_split += 1
        else:
            skipped += 1
    if skipped:
        print(f"skipped {skipped} ambiguous entries", flush=True)
    if n_split:
        print(f"split {n_split} sub-complexes from multi-complex files", flush=True)

    graphs = []
    labels = []
    pdbs = []

    def add(residues: list[Residue], label: int, pdb: str,
            forced_nodes: set | None = None):
        try:
            g = build_interface_graph_from_residues(residues,
                                                     forced_nodes=forced_nodes)
        except ValueError:
            return
        if g.n_nodes < 20:
            return
        graphs.append((g.node_features, g.node_coords, g.edge_index, g.edge_features))
        labels.append(label)
        pdbs.append(pdb)
        return g

    # positives: native complexes (capture node sets for decoy pinning)
    n_pos = 0
    native_keys: dict[str, set] = {}
    for pdb, residues in complexes.items():
        g = add(residues, 1, pdb)
        if g is not None:
            native_keys[pdb] = set(g.node_keys)
            n_pos += 1

    # negatives: displaced-peptide decoys — SAME node set as native,
    # only peptide coordinates move; no node-set/flag leakage possible
    n_neg = 0
    for pdb, residues in complexes.items():
        if pdb not in native_keys:
            continue
        decoy = make_displaced_peptide_decoy(residues)
        add(decoy, 0, f"{pdb}-disp", forced_nodes=native_keys[pdb])
        n_neg += 1

    # extra negatives: graft decoys (foreign peptide, weak signal — hard task)
    n_graft = 0
    pdb_list = list(complexes.keys())
    for pdb, residues in complexes.items():
        if pdb not in native_keys:
            continue
        peptides = _peptide_chains(residues)
        if not peptides:
            continue
        native_seq = "".join(r.resname[:1] for r in next(iter(peptides.values())))
        # pick a donor peptide with a DIFFERENT sequence (else decoy == native)
        graft = None
        donor_pdb = None
        for _ in range(20):
            donor_pdb = rng.choice(pdb_list)
            donor_peptides = _peptide_chains(complexes[donor_pdb])
            if not donor_peptides:
                continue
            candidate = rng.choice(list(donor_peptides.values()))
            if "".join(r.resname[:1] for r in candidate) != native_seq:
                graft = candidate
                break
        if graft is None:
            continue
        decoy = make_decoy(residues, graft)
        add(decoy, 0, f"{pdb}<-{donor_pdb}", forced_nodes=native_keys[pdb])
        n_graft += 1

    import torch  # lazy: graph logic above is torch-free

    torch.save(
        {"graphs": graphs, "labels": labels, "pdbs": pdbs},
        out_path,
    )
    return {"positives": n_pos, "negatives": n_neg + n_graft,
            "displaced": n_neg, "graft": n_graft, "total": len(labels)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct-dir", type=Path, default=Path("data/raw/structures"))
    ap.add_argument("--out", type=Path, default=Path("data/processed/egnn_dataset.pt"))
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    stats = build_dataset(args.struct_dir, args.out, args.seed)
    print(f"dataset saved to {args.out}: {stats}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
