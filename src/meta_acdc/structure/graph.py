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
    icode: str = ""  # insertion code (CDR3 regions: resid 100A/100B...)
    atoms: list[tuple[str, float, float, float]] = field(default_factory=list)
    # (atom name, x, y, z) — ALL heavy atoms (backbone N/CA/C/O + side chains)


@dataclass
class InterfaceGraph:
    """One TCR-pMHC interface as node/edge tensors (lists at v0; torch later)."""

    n_nodes: int = 0
    node_features: list[list[float]] = field(default_factory=list)  # phys + pLDDT
    node_coords: list[tuple[float, float, float]] = field(default_factory=list)
    edge_index: list[tuple[int, int]] = field(default_factory=list)
    edge_features: list[list[float]] = field(default_factory=list)  # RBF distances
    residue_info: list[tuple[str, str, int]] = field(default_factory=list)
    node_keys: list[tuple[str, int, str]] = field(default_factory=list)  # (chain, resid, icode)

    def summary(self) -> str:
        n_tcr = sum(1 for c, _, _ in self.residue_info if is_tcr_chain(c))
        n_pmhc = self.n_nodes - n_tcr
        return f"nodes={self.n_nodes} (TCR={n_tcr}, pMHC={n_pmhc}), edges={len(self.edge_index)}"


# --- chain classification -------------------------------------------------
# Chain LETTERS are unreliable across PDB entries (1AO7: TCR on D/E, MHC on A;
# 1BD2: TCR on A/B). Classify by sequence content instead:
#   - TCR chains carry the conserved "YFC" motif immediately before CDR3
#     (beta: YFCAS/YLCAS..., alpha: YFCAV/YFCAL...)
#   - MHC class I heavy chain ~275 residues; class II ~180-230
#   - beta-2 microglobulin ~99-100 residues
#   - peptide 5-20 residues

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M",
}

# roles
TCR = "tcr"
MHC = "mhc"
B2M = "b2m"
PEPTIDE = "peptide"
OTHER = "other"


def is_tcr_chain(chain: str) -> bool:
    """Legacy chain-letter heuristic (kept for backward compat; prefer
    classify_chains for correctness)."""
    return chain in {"A", "B"}


def is_pmhc_chain(chain: str) -> bool:
    return not is_tcr_chain(chain)


def classify_chains(residues: list[Residue]) -> dict[str, str]:
    """Content-based chain role classification."""
    # per-chain sequence in file order (insertions included)
    last: dict[str, tuple[int, str]] = {}
    seqs: dict[str, list[str]] = {}
    for r in residues:
        key = (r.resid, r.icode)
        if last.get(r.chain) == key:
            continue
        last[r.chain] = key
        seqs.setdefault(r.chain, []).append(AA3TO1.get(r.resname[:3].upper(), "X"))
    seqs = {c: "".join(s) for c, s in seqs.items()}

    roles: dict[str, str] = {}
    for chain, seq in seqs.items():
        n = len(seq)
        if "YFC" in seq or "YLC" in seq:
            roles[chain] = TCR
        elif n > 250:
            roles[chain] = MHC
        elif 90 <= n <= 105:
            roles[chain] = B2M
        elif 5 <= n <= 20:
            roles[chain] = PEPTIDE
        elif 150 <= n <= 250:
            roles[chain] = MHC  # class II
        else:
            roles[chain] = OTHER
    return roles


# --- PDB parsing (minimal, stdlib) ----------------------------------------

def parse_pdb(path: Path) -> list[Residue]:
    """Extract one record per residue (centroid of ATOM lines; CA preferred)."""
    atoms: dict[tuple[str, int, str], list[tuple[float, float, float]]] = {}
    plddt: dict[tuple[str, int, str], float | None] = {}
    resname: dict[tuple[str, int, str], str] = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resn = line[17:20].strip()
            if resn in ("HOH", "WAT"):  # skip water (HETATM pollution)
                continue
            if line.startswith("HETATM"):  # skip ligands/sugars entirely
                continue
            altloc = line[16].strip()
            if altloc not in ("", "A"):  # take only the A alternate conformation
                continue
            chain = line[21].strip() or line[72:76].strip()
            resid = int(line[22:26])
            icode = line[26].strip()
            name = line[12:16].strip()
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            bfac = float(line[60:66]) if len(line) > 66 else 0.0
            key = (chain, resid, icode)
            atoms.setdefault(key, []).append((name, x, y, z))
            resname.setdefault(key, resn)
            # treat B-factor as pLDDT when it looks like a confidence score
            if bfac > 0:
                plddt[key] = min(bfac, 100.0)

    residues = []
    for (chain, resid, icode), coords in atoms.items():
        cx = sum(c[1] for c in coords) / len(coords)
        cy = sum(c[2] for c in coords) / len(coords)
        cz = sum(c[3] for c in coords) / len(coords)
        # all heavy atoms (H atoms filtered at collection time below)
        heavy = [(n, x, y, z) for n, x, y, z in coords if not n.startswith("H")]
        residues.append(Residue(chain, resname[(chain, resid, icode)], resid, cx, cy, cz,
                                plddt.get((chain, resid, icode)), icode, heavy))
    return residues


# --- graph construction ----------------------------------------------------

def _rbf(dist: float) -> list[float]:
    return [math.exp(-((dist - (i * EDGE_CUTOFF / RBF_CENTERS)) ** 2) / 2.0)
            for i in range(RBF_CENTERS)]


CONTACT_CUTOFF = 4.5  # A, van-der-Waals contact


def _is_carbon(name: str) -> bool:
    return name.upper().startswith("C")


BACKBONE_NAMES = ("N", "CA", "C", "O", "OXT")


def _sidechain_atoms(atoms: list[tuple]) -> list[tuple]:
    return [a for a in atoms if a[0] not in BACKBONE_NAMES]


def _contact_hist(a_atoms: list[tuple], b_atoms: list[tuple]) -> list[float]:
    """Atom-contact histogram between two residues' SIDE CHAINS.

    Bins: [C-C, C-hetero, hetero-hetero, total] counts of atom pairs within
    CONTACT_CUTOFF. Captures side-chain packing — the C-alpha-invisible
    signal that graft decoys alter.
    """
    a_sc = _sidechain_atoms(a_atoms)
    b_sc = _sidechain_atoms(b_atoms)
    cc = ch = hh = 0
    for na, ax, ay, az in a_sc:
        for nb, bx, by, bz in b_sc:
            if (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2 > CONTACT_CUTOFF ** 2:
                continue
            ca, cb = _is_carbon(na), _is_carbon(nb)
            if ca and cb:
                cc += 1
            elif ca != cb:
                ch += 1
            else:
                hh += 1
    return [float(cc), float(ch), float(hh), float(cc + ch + hh)]


AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"  # 20 canonical residues


def _sidechain_extent(res: Residue) -> float:
    """Distance from CA centroid to side-chain centroid (0 for GLY)."""
    sc = _sidechain_atoms(res.atoms)
    if not sc:
        return 0.0
    cx = sum(a[1] for a in sc) / len(sc)
    cy = sum(a[2] for a in sc) / len(sc)
    cz = sum(a[3] for a in sc) / len(sc)
    return math.sqrt((cx - res.x) ** 2 + (cy - res.y) ** 2 + (cz - res.z) ** 2)


def _node_features(res: Residue) -> list[float]:
    aa = res.resname[:1].upper()
    props = AA_PHYSICOCHEM.get(aa, AA_PHYSICOCHEM["X"])
    conf = res.plddt / 100.0 if res.plddt is not None else 0.5
    onehot = [1.0 if c == aa else 0.0 for c in AA_ORDER]
    return [*props, conf, *onehot]  # 3 phys + pLDDT + 20-dim identity


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
    forced_nodes: set[tuple[str, int, str]] | None = None,
) -> InterfaceGraph:
    """Build an interface graph from an explicit residue list.

    Chain roles drive the side split (content-based classification). If
    forced_nodes is given (decoy generation), the node set is pinned to those
    (chain, resid, icode) keys — only coordinates/edges change, so node
    features cannot leak the label via node-set differences.
    """
    if len(residues) < 10:
        raise ValueError(f"{pdb_path or '<residues>'}: too few residues ({len(residues)})")

    # split complex sides by CONTENT-based chain roles
    roles = classify_chains(residues)
    tcr = [r for r in residues if roles[r.chain] == TCR]
    pmhc = [r for r in residues if roles[r.chain] != TCR]
    if not tcr or not pmhc:
        raise ValueError(f"{pdb_path or '<residues>'}: cannot split TCR/pMHC chains "
                         f"(tcr={len(tcr)}, pmhc={len(pmhc)}); roles={roles}")

    def dist(a: Residue, b: Residue) -> float:
        return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)

    if forced_nodes is not None:
        by_key = {(r.chain, r.resid, r.icode): i for i, r in enumerate(residues)}
        keep = {by_key[k] for k in forced_nodes if k in by_key}
        if len(keep) < 5:
            keep = set(range(len(residues)))
    else:
        # interface selection: within interface_radius of the opposite side
        keep = set()
        for i, r in enumerate(residues):
            is_t = roles[r.chain] == TCR
            other = pmhc if is_t else tcr
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
    graph.node_keys = [(residues[i].chain, residues[i].resid, residues[i].icode)
                       for i in idx]
    # peptide-chain flag from role classification
    peptide_chains = {c for c, r in roles.items() if r == PEPTIDE}
    graph.node_features = [
        [*_node_features(residues[i]),
         1.0 if residues[i].chain in peptide_chains else 0.0,
         _sidechain_extent(residues[i])]
        for i in idx
    ]
    graph.node_coords = [(residues[i].x, residues[i].y, residues[i].z) for i in idx]

    for new_i, old_i in enumerate(idx):
        for new_j, old_j in enumerate(idx):
            if new_j <= new_i:
                continue
            d = dist(residues[old_i], residues[old_j])
            if d <= edge_cutoff:
                contact = _contact_hist(residues[old_i].atoms,
                                        residues[old_j].atoms)
                feats = [*_rbf(d), *contact]  # 12 RBF + 4 contact = 16
                graph.edge_index.append((new_i, new_j))
                graph.edge_index.append((new_j, new_i))  # undirected
                graph.edge_features.append(feats)
                graph.edge_features.append(feats)
    return graph
