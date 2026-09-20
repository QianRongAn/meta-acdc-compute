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
    parse_pdb,
    is_tcr_chain,
)

PEPTIDE_LEN = (5, 20)


def _split_sides(residues: list[Residue]) -> tuple[list[Residue], list[Residue]]:
    tcr = [r for r in residues if is_tcr_chain(r.chain)]
    pmhc = [r for r in residues if not is_tcr_chain(r.chain)]
    return tcr, pmhc


def _peptide_chains(residues: list[Residue]) -> dict[str, list[Residue]]:
    """Group pMHC-side residues by chain; peptide chains have 5-20 residues."""
    chains: dict[str, list[Residue]] = {}
    for r in residues:
        if not is_tcr_chain(r.chain):
            chains.setdefault(r.chain, []).append(r)
    return {c: rs for c, rs in chains.items() if PEPTIDE_LEN[0] <= len(rs) <= PEPTIDE_LEN[1]}


def make_decoy(
    residues: list[Residue],
    peptide: list[Residue],
) -> list[Residue]:
    """Graft a foreign peptide chain in place of the native one.

    Keeps the peptide's own chain label ('P') but its coordinates/identity
    come from the foreign complex; geometry near the groove is approximate —
    acceptable for decoy negatives at v0.
    """
    out = [r for r in residues if r.chain != "P"]  # drop native peptide chain
    # translate foreign peptide to sit near the native groove (use native B2M/MHC centroid)
    mhc = [r for r in out if not is_tcr_chain(r.chain)]
    if not mhc:
        return out + peptide
    cx = sum(r.x for r in mhc) / len(mhc)
    cy = sum(r.y for r in mhc) / len(mhc)
    cz = sum(r.z for r in mhc) / len(mhc)
    for r in peptide:
        out.append(Residue("P", r.resname, r.resid, r.x - cx, r.y - cy, r.z - cz + 2.0,
                           r.plddt))
    return out


def make_displaced_peptide_decoy(
    residues: list[Residue],
    shift: float = 15.0,
) -> list[Residue]:
    """Negative decoy: the native peptide chain rigidly translated out of the
    groove (along +z). Real interface geometry destroyed — a strong structural
    signal for pipeline validation. (The graft decoy is the weak-signal variant
    we aim to learn later with better readouts.)
    """
    chain_sizes: dict[str, int] = {}
    for r in residues:
        chain_sizes[r.chain] = chain_sizes.get(r.chain, 0) + 1
    peptide_chains = {
        c for c, n in chain_sizes.items()
        if not is_tcr_chain(c) and PEPTIDE_LEN[0] <= n <= PEPTIDE_LEN[1]
    }
    out = []
    for r in residues:
        if r.chain in peptide_chains:
            out.append(Residue(r.chain, r.resname, r.resid, r.x, r.y, r.z + shift, r.plddt))
        else:
            out.append(r)
    return out


def build_dataset(struct_dir: Path, out_path: Path, seed: int = 42) -> dict:
    rng = random.Random(seed)
    with open(struct_dir / "complexes.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    complexes: dict[str, list[Residue]] = {}
    skipped = 0
    for row in rows:
        pdb = row["pdb"]
        path = struct_dir / f"{pdb.lower()}.pdb"
        try:
            residues = parse_pdb(path)
        except Exception:
            continue
        # keep only clean single complexes: exactly 2 TCR chains (A,B),
        # exactly one peptide chain (5-20 res) and >=1 MHC-like chain (>100 res)
        chain_sizes: dict[str, int] = {}
        for r in residues:
            chain_sizes[r.chain] = chain_sizes.get(r.chain, 0) + 1
        tcr_chains = [c for c in ("A", "B") if c in chain_sizes]
        other = {c: n for c, n in chain_sizes.items() if c not in ("A", "B")}
        n_peptide = sum(1 for n in other.values() if 5 <= n <= 20)
        n_mhc = sum(1 for n in other.values() if n > 100)
        if len(tcr_chains) < 2 or n_peptide != 1 or n_mhc < 1:
            skipped += 1
            continue
        complexes[pdb] = residues
    if skipped:
        print(f"skipped {skipped} multi-complex / ambiguous entries", flush=True)

    graphs = []
    labels = []
    pdbs = []

    def add(residues: list[Residue], label: int, pdb: str):
        try:
            g = build_interface_graph_from_residues(residues)
        except ValueError:
            return
        if g.n_nodes < 20:
            return
        graphs.append((g.node_features, g.node_coords, g.edge_index, g.edge_features))
        labels.append(label)
        pdbs.append(pdb)

    # positives: native complexes
    n_pos = 0
    for pdb, residues in complexes.items():
        add(residues, 1, pdb)
        n_pos += 1

    # negatives: displaced-peptide decoys (strong signal, pipeline validation)
    n_neg = 0
    for pdb, residues in complexes.items():
        decoy = make_displaced_peptide_decoy(residues)
        add(decoy, 0, f"{pdb}-disp")
        n_neg += 1

    # extra negatives: graft decoys (foreign peptide, weak signal — hard task)
    n_graft = 0
    pdb_list = list(complexes.keys())
    for pdb, residues in complexes.items():
        peptides = _peptide_chains(residues)
        if not peptides:
            continue
        donor = rng.choice(pdb_list)
        donor_peptides = _peptide_chains(complexes[donor])
        if not donor_peptides:
            continue
        graft = rng.choice(list(donor_peptides.values()))
        decoy = make_decoy(residues, graft)
        add(decoy, 0, f"{pdb}<-{donor}")
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
