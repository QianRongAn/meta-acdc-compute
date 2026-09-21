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


def make_decoy(
    residues: list[Residue],
    peptide: list[Residue],
) -> list[Residue]:
    """Graft a foreign peptide chain in place of the native one.

    TRUE structural-mimicry decoy: the foreign peptide's CHEMISTRY is written
    onto the NATIVE peptide's coordinates. Geometry, node set and chain labels
    are untouched — only peptide node features change. The model must detect
    chemical incompatibility from features + edge interactions, with no
    geometric shortcut. (Peptide lengths may differ: swap the aligned prefix.)
    """
    roles = classify_chains(residues)
    native_pep_chains = [c for c, r in roles.items() if r == PEPTIDE]
    if not native_pep_chains:
        return residues
    native_pep_chain = native_pep_chains[0]
    native_pep = [r for r in residues if r.chain == native_pep_chain]
    # pair by position: foreign chemistry onto native coordinates
    swap = {id(r): f for r, f in zip(native_pep, peptide)}
    out = []
    for r in residues:
        f = swap.get(id(r))
        if f is not None:
            # keep native position/chain/resid, take foreign chemistry
            out.append(Residue(r.chain, f.resname, r.resid, r.x, r.y, r.z,
                               r.plddt, r.icode))
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
            out.append(Residue(r.chain, r.resname, r.resid, r.x, r.y, r.z + shift,
                               r.plddt, r.icode))
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
            roles = classify_chains(residues)
        except Exception:
            continue
        # keep only clean single complexes: exactly 2 TCR chains,
        # exactly one peptide chain, >=1 MHC chain
        n_tcr = sum(1 for v in roles.values() if v == TCR)
        n_peptide = sum(1 for v in roles.values() if v == PEPTIDE)
        n_mhc = sum(1 for v in roles.values() if v == MHC)
        if n_tcr != 2 or n_peptide != 1 or n_mhc < 1:
            skipped += 1
            continue
        complexes[pdb] = residues
    if skipped:
        print(f"skipped {skipped} multi-complex / ambiguous entries", flush=True)

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
        donor = rng.choice(pdb_list)
        donor_peptides = _peptide_chains(complexes[donor])
        if not donor_peptides:
            continue
        graft = rng.choice(list(donor_peptides.values()))
        decoy = make_decoy(residues, graft)
        add(decoy, 0, f"{pdb}<-{donor}", forced_nodes=native_keys[pdb])
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
