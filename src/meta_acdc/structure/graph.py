"""Interface-graph construction from TCR-pMHC structures (proposal Module 1).

Plan-of-record decisions (lit review 02):
- Parse PDB files with a minimal stdlib parser (no Biopython dependency at v0).
- Node = amino-acid residue (per chain). Node features: physicochemical
  properties (charge, hydrophobicity, side-chain volume) + AlphaFold-derived
  confidence (pLDDT) when available.
- Interface extraction: residues within INTERFACE_RADIUS (10 A), measured between heavy-atom centroids, of
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

import numpy as np
from .feature_schema import AA_ORDER, validate_node_features

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
    plddt: float | None = None  # only for an explicitly identified predicted structure
    icode: str = ""  # insertion code (CDR3 regions: resid 100A/100B...)
    atoms: list[tuple[str, float, float, float]] = field(default_factory=list)
    role: str | None = None
    chain_kind: str | None = None
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
    peptide_mask: list[bool] = field(default_factory=list)
    node_roles: list[str] = field(default_factory=list)
    # Diagnostic mode may retain incomplete nodes; accepted builds reject these.
    missing_atom_nodes: list[dict] = field(default_factory=list)
    node_keys: list[tuple[str, int, str]] = field(default_factory=list)  # (chain, resid, icode)

    def summary(self) -> str:
        if len(self.node_roles) != self.n_nodes:
            raise ValueError("explicit node roles required for graph summary")
        n_tcr = self.node_roles.count(TCR)
        n_pmhc = sum(r in {MHC, B2M, PEPTIDE} for r in self.node_roles)
        if n_tcr + n_pmhc != self.n_nodes:
            raise ValueError("unknown node role")
        return f"nodes={self.n_nodes} (TCR={n_tcr}, pMHC={n_pmhc}), edges={len(self.edge_index)}"


# --- chain classification -------------------------------------------------
# Historical heuristics removed; chain letters and motifs do not prove roles.
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


def is_tcr_chain(chain: str, roles: dict[str, str] | None = None) -> bool:
    """Require verified roles; chain letters carry no biological identity."""
    if roles is None or chain not in roles:
        raise ValueError("explicit chain roles required")
    return roles[chain] == TCR


def is_pmhc_chain(chain: str, roles: dict[str, str] | None = None) -> bool:
    if roles is None or chain not in roles:
        raise ValueError("explicit chain roles required")
    return roles[chain] in {MHC, B2M, PEPTIDE}


def classify_chains(residues: list[Residue]) -> dict[str, str]:
    """Use depositor/reference annotations; motifs and lengths cannot prove identity."""
    roles = {}
    for chain in {r.chain for r in residues}:
        explicit = {r.role for r in residues if r.chain == chain and r.role is not None}
        if len(explicit) > 1:
            raise ValueError("conflicting chain role annotations")
        role = next(iter(explicit), OTHER)
        if role not in {TCR, MHC, B2M, PEPTIDE, OTHER}:
            raise ValueError("invalid chain role")
        roles[chain] = role
    return roles


# --- PDB parsing (minimal, stdlib) ----------------------------------------

def parse_pdb(path: Path) -> list[Residue]:
    """Extract one record per residue; experimental B factors are not pLDDT."""
    atoms: dict[tuple[str, int, str], list[tuple[float, float, float]]] = {}
    plddt: dict[tuple[str, int, str], float | None] = {}
    resname: dict[tuple[str, int, str], str] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"): break
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
            # Crystallographic B factors are not pLDDT.
            plddt[key] = None

    residues = []
    for (chain, resid, icode), coords in atoms.items():
        cx = sum(c[1] for c in coords) / len(coords)
        cy = sum(c[2] for c in coords) / len(coords)
        cz = sum(c[3] for c in coords) / len(coords)
        # all heavy atoms (H atoms filtered at collection time below)
        heavy = [(n, x, y, z) for n, x, y, z in coords if not n.startswith("H")]
        residues.append(Residue(chain, resname[(chain, resid, icode)], resid, cx, cy, cz,
                                plddt.get((chain, resid, icode)), icode, heavy))
    from .identity import annotate
    annotate(residues, path)
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




def _sidechain_extent(res: Residue) -> float:
    """Distance from the CA atom to the side-chain centroid (0 for GLY).

    The CA reference is taken from the residue's CURRENT atoms: grafted
    residues keep the receiver's backbone but gain foreign side chains, so
    the stored (x, y, z) centroid from parse_pdb no longer matches the atoms
    and would put the extent on a stale coordinate frame (AUDIT.md #5).
    """
    sc = _sidechain_atoms(res.atoms)
    if not sc:
        return 0.0
    ca = [a for a in res.atoms if a[0].upper() == "CA"]
    if ca:
        bx = sum(a[1] for a in ca) / len(ca)
        by = sum(a[2] for a in ca) / len(ca)
        bz = sum(a[3] for a in ca) / len(ca)
    else:
        bx, by, bz = res.x, res.y, res.z
    cx = sum(a[1] for a in sc) / len(sc)
    cy = sum(a[2] for a in sc) / len(sc)
    cz = sum(a[3] for a in sc) / len(sc)
    return math.sqrt((cx - bx) ** 2 + (cy - by) ** 2 + (cz - bz) ** 2)


def _node_features(res: Residue) -> list[float]:
    # Map the 3-letter name through AA3TO1. The previous `res.resname[:1]`
    # collapsed ARG/ASN/ASP onto A, GLU onto G, LYS onto L, PHE onto P and
    # TRP/TYR onto T, so charged and aromatic residues inherited wrong
    # physicochemical values and identity codes (LYS looked like LEU). That is
    # an input error, not a display bug. See AUDIT.md #2.
    aa = AA3TO1.get(res.resname.upper(), "X")
    props = AA_PHYSICOCHEM.get(aa, AA_PHYSICOCHEM["X"])
    conf = res.plddt / 100.0 if res.plddt is not None else 0.5
    onehot = [1.0 if c == aa else 0.0 for c in AA_ORDER]
    return [*props, conf, *onehot]  # 3 phys + pLDDT + 20-dim identity


def build_interface_graph(
    pdb_path: Path,
    interface_radius: float = INTERFACE_RADIUS,
    edge_cutoff: float = EDGE_CUTOFF,
    electrostatics: bool = False,
) -> InterfaceGraph:
    """Build an interface graph from a TCR-pMHC PDB file."""
    residues = parse_pdb(pdb_path)
    return build_interface_graph_from_residues(
        residues, pdb_path, interface_radius, edge_cutoff,
        electrostatics=electrostatics,
    )


def residue_centroid(res: Residue) -> tuple[float, float, float]:
    """Centroid of the residue's CURRENT heavy atoms.

    parse_pdb stores the full-atom centroid as (x, y, z), but a grafted
    residue keeps the receiver's backbone and gains foreign side chains, so
    its stored coordinates no longer match its atoms. Keeping the stale
    centroid meant a "side-chain-only" change also moved the node by up to
    several Å, which invalidates the claim that geometry is pinned, and the
    edge distances disagreed with the reported node coordinates.
    (AUDIT.md #5.)
    """
    heavy = [a for a in res.atoms if not a[0].startswith("H")]
    if not heavy:
        return (res.x, res.y, res.z)
    n = len(heavy)
    return (sum(a[1] for a in heavy) / n,
            sum(a[2] for a in heavy) / n,
            sum(a[3] for a in heavy) / n)


def residue_dist(a: Residue, b: Residue) -> float:
    ca = residue_centroid(a)
    cb = residue_centroid(b)
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(ca, cb)))


def build_interface_graph_from_residues(
    residues: list[Residue],
    pdb_path: Path | None = None,
    interface_radius: float = INTERFACE_RADIUS,
    edge_cutoff: float = EDGE_CUTOFF,
    forced_nodes: set[tuple[str, int, str]] | None = None,
    electrostatics: bool = False,
    require_complete_atoms: bool = True,
) -> InterfaceGraph:
    """Build an interface graph from an explicit residue list.

    Chain roles drive the side split (content-based classification). If
    forced_nodes is given (decoy generation), the node set is pinned to those
    (chain, resid, icode) keys — only coordinates/edges change, so node
    features cannot leak the label via node-set differences.

    electrostatics=True appends the v11 atom-level features (see
    electrostatics.py): node += [esp, burial], edge += [coul, coul_mag].
    """
    if len(residues) < 10:
        raise ValueError(f"{pdb_path or '<residues>'}: too few residues ({len(residues)})")

    residues = sorted(residues, key=lambda r: (r.chain, r.resid, r.icode))
    identities = [(r.chain, r.resid, r.icode) for r in residues]
    if len(set(identities)) != len(identities):
        raise ValueError("duplicate residue identity")
    # split complex sides by explicit chain roles
    roles = classify_chains(residues)
    if OTHER in roles.values():
        raise ValueError("unverified chain identity: isolate an explicitly annotated complex first")
    tcr = [r for r in residues if roles[r.chain] == TCR]
    pmhc = [r for r in residues if roles[r.chain] != TCR]
    if not tcr or not pmhc:
        raise ValueError(f"{pdb_path or '<residues>'}: cannot split TCR/pMHC chains "
                         f"(tcr={len(tcr)}, pmhc={len(pmhc)}); roles={roles}")

    centers = np.asarray([residue_centroid(r) for r in residues])
    distances = np.linalg.norm(centers[:, None, :] - centers[None, :, :], axis=-1)

    def dist(a: Residue, b: Residue) -> float:
        # current-atom centroids: stale stored coords would disagree with
        # the node coordinates reported in the output graph (AUDIT.md #5)
        return residue_dist(a, b)

    by_key = {(r.chain, r.resid, r.icode): i for i, r in enumerate(residues)}
    if forced_nodes is not None:
        if not forced_nodes or not forced_nodes <= set(by_key):
            raise ValueError("pinned node identities missing from perturbed structure")
        keep = {by_key[k] for k in forced_nodes}
    else:
        is_t = np.asarray([roles[r.chain] == TCR for r in residues])
        keep = set(np.where(np.where(is_t, distances[:, ~is_t].min(1), distances[:, is_t].min(1)) <= interface_radius)[0].tolist())
        if len(keep) < 5:
            raise ValueError("no valid interface; full-complex fallback forbidden")
    if not any(roles[residues[i].chain] == PEPTIDE for i in keep):
        raise ValueError("interface contains no peptide nodes")

    idx = sorted(keep)
    from .identity import missing_heavy_atoms
    incomplete = [(residues[i].chain, residues[i].resid, missing_heavy_atoms(residues[i]))
                  for i in idx if missing_heavy_atoms(residues[i])]
    if incomplete and require_complete_atoms:
        raise ValueError(f"incomplete interface heavy atoms: {incomplete[:5]}")
    pos = {old: new for new, old in enumerate(idx)}
    graph = InterfaceGraph(n_nodes=len(idx))
    graph.residue_info = [(residues[i].chain, residues[i].resname, residues[i].resid)
                          for i in idx]
    graph.missing_atom_nodes = [{"chain": c, "resid": r, "missing": m} for c, r, m in incomplete]
    graph.node_keys = [(residues[i].chain, residues[i].resid, residues[i].icode)
                       for i in idx]
    graph.node_roles = [roles[residues[i].chain] for i in idx]
    graph.peptide_mask = [x == PEPTIDE for x in graph.node_roles]
    # peptide-chain flag from role classification
    peptide_chains = {c for c, r in roles.items() if r == PEPTIDE}

    def _centroid(res: Residue) -> tuple[float, float, float]:
        """Centroid of the residue's CURRENT heavy atoms (see residue_centroid)."""
        return residue_centroid(res)

    graph.node_features = [
        [*_node_features(residues[i]),
         1.0 if residues[i].chain in peptide_chains else 0.0,
         _sidechain_extent(residues[i])]
        for i in idx
    ]
    validate_node_features(graph.node_features)
    graph.node_coords = [_centroid(residues[i]) for i in idx]

    # --- v11 atom-level electrostatics (opt-in, appended features) --------
    es = None
    if electrostatics:
        from meta_acdc.structure.electrostatics import (
            burial_feature, edge_coulomb, esp_feature, residue_charged_atoms)
        # first residue of each chain (file order) -> N-terminal charges
        seen_chains: set[str] = set()
        first_keys: set[tuple[str, int, str]] = set()
        for r in residues:
            if r.chain not in seen_chains:
                seen_chains.add(r.chain)
                first_keys.add((r.chain, r.resid, r.icode))
        charge_rows = residue_charged_atoms(residues, first_keys)
        # stacked source atoms + per-residue slices (for self-exclusion)
        slices, offset = [], 0
        for rows in charge_rows:
            slices.append((offset, offset + rows.shape[0]))
            offset += rows.shape[0]
        sources = (np.vstack(charge_rows)
                   if charge_rows and offset else
                   np.empty((0, 4)))
        # opposite-side atom pools for the desolvation proxy
        tcr_idx = [i for i, r in enumerate(residues) if roles[r.chain] == TCR]
        pmhc_idx = [i for i, r in enumerate(residues) if roles[r.chain] != TCR]
        tcr_atoms = np.vstack([charge_rows[i] for i in tcr_idx]) \
            if tcr_idx else np.empty((0, 4))
        pmhc_atoms = np.vstack([charge_rows[i] for i in pmhc_idx]) \
            if pmhc_idx else np.empty((0, 4))
        es = {"sources": sources, "slices": slices,
              "tcr_atoms": tcr_atoms, "pmhc_atoms": pmhc_atoms,
              "charge_rows": charge_rows}
        for new_i, old_i in enumerate(idx):
            is_t = roles[residues[old_i].chain] == TCR
            opposite = pmhc_atoms if is_t else tcr_atoms
            esp = esp_feature(graph.node_coords[new_i], sources, slices[old_i])
            bur = burial_feature(charge_rows[old_i], opposite)
            graph.node_features[new_i].extend([esp, bur])

    for new_i, old_i in enumerate(idx):
        for new_j, old_j in enumerate(idx):
            if new_j <= new_i:
                continue
            d = float(distances[old_i, old_j])
            if d <= edge_cutoff:
                contact = _contact_hist(residues[old_i].atoms,
                                        residues[old_j].atoms)
                feats = [*_rbf(d), *contact]  # 12 RBF + 4 contact = 16
                if es is not None:
                    coul, mag = edge_coulomb(es["charge_rows"][old_i],
                                             es["charge_rows"][old_j])
                    feats.extend([coul, mag])  # + 2 electrostatics = 18
                graph.edge_index.append((new_i, new_j))
                graph.edge_index.append((new_j, new_i))  # undirected
                graph.edge_features.append(feats)
                graph.edge_features.append(feats)
    return graph
