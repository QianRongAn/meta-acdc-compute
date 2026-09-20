"""Interface-graph construction from TCR-pMHC structures (proposal Module 1).

Plan-of-record decisions (lit review 02):
- Parse PDB files with a minimal stdlib parser (no Biopython dependency at v0).
- Node = amino-acid residue (per chain). Node features: physicochemical
  properties (charge, hydrophobicity, side-chain volume) + AlphaFold-derived
  confidence (pLDDT) when available.
- Interface extraction: residues within INTERFACE_RADIUS (10 A) of ANY atom of
  the opposite complex side (TCR chains vs pMHC chains), per the proposal.
- Edges: Euclidean distance between residue centroids < EDGE_CUTOFF (8 A).
  Edge features: RBF-encoded distance (GVP-style) + sequence separation.

Complex side assignment:
- TCR: chains whose labels are alpha/beta (common PDB labels: A/B, D/E, ...)
  identified by naming convention (TRA/TRB in chain ids or COMPND records).
- pMHC: MHC chain (HLA-A/B/C...) + peptide chain (usually chain C or P).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

# Physicochemical features per standard amino acid (indexed by 1-letter code)
# order: [charge, hydrophobicity (Kyte-Doolittle), side-chain volume (A^3)]
AA_PHYSICOCHEM = {
    "A": [0.0, 1.8, 88.6], "R": [1.0, -4.5, 173.4], "N": [0.0, -3.5, 114.1],
    "D": [-1.0, -3.5, 111.1], "C": [0.0, 2.5, 108.5], "Q": [0.0, -3.5, 143.8],
    "E": [-1.0, -3.5, 138.4], "G": [0.0, -0.4, 60.1], "H": [0.1, -3.2, 153.2],
    "I": [0.0, 4.5, 166.7], "L": [0.0, 3.8, 166.7], "K": [1.0, -3.9, 168.6],
    "M": [0.0, 1.9, 162.9], "F": [0.0, 2.8, 189.9], "P": [0.0, -1.6, 112.7],
    "S": [0.0, -0.8, 89.0], "T": [0.0, -0.7, 116.1], "W": [0.0, -0.9, 227.8],
    "Y": [0.0, -1.3, 193.6], "V": [0.0, 4.2, 140.0], "X": [0.0, 0.0, 150.0],
}

INTERFACE_RADIUS = 10.0  # A
EDGE_CUTOFF = 8.0       # A
RBF_CENTERS = 12        # GVP-style distance encoding


@dataclass
class Residue:
    chain: str
    resname: str
    resid: int
    x: float
    y: float
    z: float
    plddt: float | None = None  # from B-factor column when available


@dataclass
class InterfaceGraph:
    """One TCR-pMHC interface as node/edge tensors (lists at v0; torch later)."""

    n_nodes: int = 0
    node_features: list[list[float]] = field(default_factory=list)  # phys + pLDDT
    edge_index: list[tuple[int, int]] = field(default_factory=list)
    edge_features: list[list[float]] = field(default_factory=list)  # RBF distances
    residue_info: list[tuple[str, str, int]] = field(default_factory=list)

    def summary(self) -> str:
        n_tcr = sum(1 for c, _, _ in self.residue_info if is_tcr_chain(c))
        n_pmhc = self.n_nodes - n_tcr
        return f"nodes={self.n_nodes} (TCR={n_tcr}, pMHC={n_pmhc}), edges={len(self.edge_index)}"


# --- chain classification -------------------------------------------------

TCR_CHAIN_HINTS = {"A", "B"}  # PDB convention: TCR alpha/beta chains are A, B
PEPTIDE_CHAIN_HINTS = {"C", "P"}
MHC_CHAIN_HINTS = {"A", "B", "C", "D", "E", "F", "G", "H"}


def is_tcr_chain(chain: str) -> bool:
    return chain in TCR_CHAIN_HINTS


def is_pmhc_chain(chain: str) -> bool:
    return chain not in TCR_CHAIN_HINTS


# --- PDB parsing (minimal, stdlib) ----------------------------------------

def parse_pdb(path: Path) -> list[Residue]:
    """Extract one record per residue (centroid of ATOM lines; CA preferred)."""
    atoms: dict[tuple[str, int], list[tuple[float, float, float]]] = {}
    plddt: dict[tuple[str, int], float | None] = {}
    resname: dict[tuple[str, int], str] = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21].strip() or line[72:76].strip()
            resid = int(line[22:26])
            name = line[12:16].strip()
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            bfac = float(line[60:66]) if len(line) > 66 else 0.0
            key = (chain, resid)
            atoms.setdefault(key, []).append((x, y, z))
            resname.setdefault(key, line[17:20].strip())
            # treat B-factor as pLDDT when it looks like a confidence score
            if bfac > 0:
                plddt[key] = min(bfac, 100.0)

    residues = []
    for (chain, resid), coords in atoms.items():
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        residues.append(Residue(chain, resname[(chain, resid)], resid, cx, cy, cz,
                                plddt.get((chain, resid))))
    return residues


# --- graph construction ----------------------------------------------------

def _rbf(dist: float) -> list[float]:
    return [math.exp(-((dist - (i * EDGE_CUTOFF / RBF_CENTERS)) ** 2) / 2.0)
            for i in range(RBF_CENTERS)]


def _node_features(res: Residue) -> list[float]:
    aa = res.resname[:1].upper()
    props = AA_PHYSICOCHEM.get(aa, AA_PHYSICOCHEM["X"])
    conf = res.plddt / 100.0 if res.plddt is not None else 0.5
    return [*props, conf]


def build_interface_graph(
    pdb_path: Path,
    interface_radius: float = INTERFACE_RADIUS,
    edge_cutoff: float = EDGE_CUTOFF,
) -> InterfaceGraph:
    """Build an interface graph from a TCR-pMHC PDB file."""
    residues = parse_pdb(pdb_path)
    return build_interface_graph_from_residues(
        residues, pdb_path, interface_radius, edge_cutoff
    )


def build_interface_graph_from_residues(
    residues: list[Residue],
    pdb_path: Path | None = None,
    interface_radius: float = INTERFACE_RADIUS,
    edge_cutoff: float = EDGE_CUTOFF,
) -> InterfaceGraph:
    """Build an interface graph from an explicit residue list.

    Chain labels drive the side split (A/B = TCR, rest = pMHC). Allows
    grafting chains (e.g., decoy peptides) before building the graph.
    """
    if len(residues) < 10:
        raise ValueError(f"{pdb_path or '<residues>'}: too few residues ({len(residues)})")

    # split complex sides by chain
    tcr = [r for r in residues if is_tcr_chain(r.chain)]
    pmhc = [r for r in residues if is_pmhc_chain(r.chain)]
    if not tcr or not pmhc:
        raise ValueError(f"{pdb_path or '<residues>'}: cannot split TCR/pMHC chains "
                         f"(tcr={len(tcr)}, pmhc={len(pmhc)})")

    def dist(a: Residue, b: Residue) -> float:
        return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)

    # interface selection: within interface_radius of the opposite side
    keep: set[int] = set()
    for i, r in enumerate(residues):
        side = tcr if is_tcr_chain(r.chain) else pmhc
        other = pmhc if is_tcr_chain(r.chain) else tcr
        if any(dist(r, o) <= interface_radius for o in other):
            keep.add(i)
    if len(keep) < 5:
        # fall back to the full complex if the interface is too small
        keep = set(range(len(residues)))

    idx = sorted(keep)
    pos = {old: new for new, old in enumerate(idx)}
    graph = InterfaceGraph(n_nodes=len(idx))
    graph.residue_info = [(residues[i].chain, residues[i].resname, residues[i].resid)
                          for i in idx]
    graph.node_features = [_node_features(residues[i]) for i in idx]

    for new_i, old_i in enumerate(idx):
        for new_j, old_j in enumerate(idx):
            if new_j <= new_i:
                continue
            d = dist(residues[old_i], residues[old_j])
            if d <= edge_cutoff:
                graph.edge_index.append((new_i, new_j))
                graph.edge_index.append((new_j, new_i))  # undirected
                graph.edge_features.append(_rbf(d))
                graph.edge_features.append(_rbf(d))
    return graph
