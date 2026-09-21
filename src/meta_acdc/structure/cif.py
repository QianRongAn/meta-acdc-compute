"""Minimal mmCIF parser for AlphaFold3 output structures.

AF3 (Server and local) writes structures in mmCIF format. We map the
_atom_site loop onto the same Residue list used by the PDB pipeline
(graph.build_interface_graph_from_residues), so AF3 predictions feed the
whole downstream stack unchanged (interface graphs, atom contacts, EGNN).

Usage:
    from meta_acdc.structure.cif import parse_cif
    residues = parse_cif(Path("1ao7_LGYGFVNYI_model_0.cif"))
"""

from __future__ import annotations

from pathlib import Path

from meta_acdc.structure.graph import Residue

# _atom_site column names as they appear in AF3 CIFs
COL_GROUP = "_atom_site.group_PDB"
COL_ATOM = "_atom_site.label_atom_id"
COL_CHAIN = "_atom_site.auth_asym_id"      # fallback: label_asym_id
COL_CHAIN2 = "_atom_site.label_asym_id"
COL_RESN = "_atom_site.auth_comp_id"       # fallback: label_comp_id
COL_RESN2 = "_atom_site.label_comp_id"
COL_RESID = "_atom_site.auth_seq_id"       # fallback: label_seq_id
COL_RESID2 = "_atom_site.label_seq_id"
COL_X = "_atom_site.Cartn_x"
COL_Y = "_atom_site.Cartn_y"
COL_Z = "_atom_site.Cartn_z"
COL_B = "_atom_site.B_iso_or_equiv"
COL_MODEL = "_atom_site.pdbx_PDB_model_num"


def parse_cif(path: Path, model: str = "1") -> list[Residue]:
    """Parse one model of an mmCIF file into Residue records.

    Returns residues keyed by (chain, resid) with heavy-atom lists, matching
    the PDB-pipeline Residue layout. AF3 confidence (pLDDT) lands in the
    B-factor slot when present.
    """
    with open(path, errors="replace") as fh:
        lines = fh.readlines()

    # locate the atom_site loop
    start = None
    for i, ln in enumerate(lines):
        if ln.startswith("_atom_site."):
            start = i
            break
    if start is None:
        raise ValueError(f"{path}: no _atom_site loop found")

    header = []
    i = start
    while i < len(lines) and lines[i].startswith("_atom_site."):
        header.append(lines[i].split()[0])
        i += 1
    # skip the '#' loop delimiter line
    while i < len(lines) and not lines[i].startswith(("_", "data_", "loop_")):
        if lines[i].strip() and not lines[i].startswith("#"):
            break
        i += 1
    data_start = i

    def col(name: str, alt: str) -> int | None:
        if name in header:
            return header.index(name)
        if alt in header:
            return header.index(alt)
        return None

    i_grp = col(COL_GROUP, COL_GROUP)
    i_atom = col(COL_ATOM, COL_ATOM)
    i_chain = col(COL_CHAIN, COL_CHAIN2)
    i_resn = col(COL_RESN, COL_RESN2)
    i_resid = col(COL_RESID, COL_RESID2)
    i_x = col(COL_X, COL_X)
    i_y = col(COL_Y, COL_Y)
    i_z = col(COL_Z, COL_Z)
    i_b = col(COL_B, COL_B)
    i_model = col(COL_MODEL, COL_MODEL)
    if None in (i_atom, i_chain, i_resn, i_resid, i_x, i_y, i_z):
        raise ValueError(f"{path}: missing required _atom_site columns "
                         f"(header={header})")

    def get(row: list[str], idx: int | None) -> str:
        if idx is None or idx >= len(row):
            return ""
        return row[idx]

    atoms: dict[tuple[str, int], list[tuple[str, float, float, float]]] = {}
    resname: dict[tuple[str, int], str] = {}
    plddt: dict[tuple[str, int], float | None] = {}

    for ln in lines[data_start:]:
        row = ln.split()
        if not row or row[0].startswith(("#", "_", "loop_", "data_", "stop_")):
            continue
        if i_grp is not None and get(row, i_grp) not in ("ATOM", ""):
            continue
        if i_model is not None and get(row, i_model) not in (model, ""):
            continue
        name = get(row, i_atom)
        chain = get(row, i_chain) or get(row, i_chain)
        if not chain:
            continue
        try:
            resid = int(get(row, i_resid))
            x = float(get(row, i_x))
            y = float(get(row, i_y))
            z = float(get(row, i_z))
        except ValueError:
            continue
        if name.startswith("H"):
            continue
        resn = get(row, i_resn)
        bfac = 0.0
        if i_b is not None:
            try:
                bfac = float(get(row, i_b))
            except ValueError:
                pass
        key = (chain, resid)
        atoms.setdefault(key, []).append((name, x, y, z))
        resname.setdefault(key, resn)
        if bfac > 0:
            plddt[key] = min(bfac, 100.0)

    residues = []
    for (chain, resid), coords in atoms.items():
        cx = sum(c[1] for c in coords) / len(coords)
        cy = sum(c[2] for c in coords) / len(coords)
        cz = sum(c[3] for c in coords) / len(coords)
        residues.append(Residue(chain, resname[(chain, resid)], resid, cx, cy, cz,
                                plddt.get((chain, resid)), "", coords))
    return residues


def parse_all_models(path: Path) -> dict[str, list[Residue]]:
    """Parse every model in an AF3 CIF (keyed by model number)."""
    with open(path, errors="replace") as fh:
        content = fh.read()
    models = set()
    import re
    for m in re.finditer(r"_atom_site\.pdbx_PDB_model_num\s+(\d+)", content):
        pass  # per-row model numbers discovered during parse
    # simpler: parse once and group
    residues = parse_cif(path, model="")
    by_model: dict[str, list[Residue]] = {}
    # parse_cif filters to one model; rerun per model by scanning unique values
    with open(path, errors="replace") as fh:
        text = fh.read()
    for m in re.findall(r"^\s*(\d+)\s+ATOM", text, re.M):
        if m not in models:
            models.add(m)
    for m in sorted(models):
        try:
            by_model[m] = parse_cif(path, model=m)
        except Exception:
            continue
    return by_model
