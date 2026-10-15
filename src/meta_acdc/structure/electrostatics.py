"""Atom-level electrostatic / desolvation interface features (v11).

Motivation (README 6.2.1): v10's RESIDUE-level Coulomb term (q_i*q_j/d from
net residue charges) was a negative result -- graft AUROC 0.933 -> 0.917.
Net residue charge is too coarse: it cannot see how a charged group is
placed relative to the partner's atoms. v11 uses FORCE-FIELD PARTIAL
CHARGES (ff14SB, united-atom; see charges.py) at the atom level:

- edge  += [coul, coul_mag]     pair Coulomb with distance-dependent
                                dielectric (eps = 4*r, kcal/mol) and its
                                unsigned magnitude (charge-density contact),
                                summed over atom pairs within 10 A.
- node  += [esp, burial]        electrostatic potential at the residue
                                centroid (all atoms within 10 A), and a
                                desolvation proxy: burial of charged atoms
                                (|q| > 0.3 e) by OPPOSITE-side heavy atoms
                                within 5 A.

All four features are squashed with tanh (bounded, so batch/per-graph
normalization stays well-behaved). They are appended AFTER the v9.1
node/edge feature blocks, so indices 0..25 (node) / 0..15 (edge) keep
their meaning -- the is-peptide flag stays at node index 4.

Physics notes:
- Distance-dependent dielectric eps(r) = 4*r is the standard implicit-solvent
  Coulomb screening (GROMACS/GROMOS convention).
- Heavy-atom-only structures (PDB) -> hydrogen charges are folded onto their
  bonded heavy atoms (united-atom), which preserves group charges exactly.
- Contributions are truncated at 10 A (with eps = 4*r they decay as ~1/r^2,
  so the truncation is mild).
"""

from __future__ import annotations

import math

import numpy as np

from meta_acdc.structure.charges import unified_charges

COULOMB_KCAL = 332.0637   # kcal/mol per (e^2 / A)
DIELECTRIC = 4.0          # eps(r) = DIELECTRIC * r
ATOM_CUTOFF = 10.0        # A, atom-pair / potential truncation
BURIAL_RADIUS = 5.0       # A, desolvation neighbor shell
CHARGED_Q = 0.3           # e, |q| threshold for "charged atom"
COUL_SQUASH = 5.0         # kcal/mol -> tanh scale (signed coulomb)
MAG_SQUASH = 10.0         # kcal/mol -> tanh scale (unsigned magnitude)
ESP_SQUASH = 5.0          # kcal/mol/e -> tanh scale (potential)
BURIAL_SQUASH = 5.0       # contacts -> tanh scale (burial count)


def residue_charged_atoms(
    residues: list,
    first_keys: set[tuple[str, int, str]],
) -> list[np.ndarray]:
    """United-atom (x, y, z, q) rows for each residue, in list order.

    Residues without any chargeable heavy atom yield an empty (0, 4) array.
    first_keys marks the first residue of each chain (N-terminal charges).
    """
    out: list[np.ndarray] = []
    for r in residues:
        key = (r.chain, r.resid, r.icode)
        qtab = unified_charges(
            r.resname,
            has_oxt=any(a[0] == "OXT" for a in r.atoms),
            first_in_chain=key in first_keys,
        )
        rows = [(x, y, z, qtab[name]) for name, x, y, z in r.atoms
                if name in qtab]
        out.append(np.asarray(rows, dtype=np.float64).reshape(-1, 4)
                   if rows else np.empty((0, 4)))
    return out


def _r2_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Squared distances between rows of (n,3) and (m,3) -> (n, m)."""
    d = a[:, None, :3] - b[None, :, :3]
    return np.einsum("nmk,nmk->nm", d, d)


def esp_feature(centroid: tuple[float, float, float],
                sources: np.ndarray, own_slice: tuple[int, int]) -> float:
    """Electrostatic potential (kcal/mol/e) at a residue centroid.

    Sources: all charged atoms of the complex (own residue's atoms excluded
    via own_slice). eps(r) = DIELECTRIC * r, so the screened potential is
    q / (eps(r) * r) ~ 1/r^2; truncation 10 A, tanh-squashed.

    The node term used q/(4r) until 2026-10-15 while the edge term and the
    docstring used 1/r^2, i.e. node and edge ran different shielding models
    and doubling the distance only halved the node potential. Both now use
    eps(r) = 4r consistently (AUDIT.md #16).
    """
    if sources.shape[0] == 0:
        return 0.0
    c = np.asarray(centroid, dtype=np.float64)[None, :]
    r2 = _r2_matrix(sources, c)[:, 0]
    mask = (r2 > 0.0) & (r2 < ATOM_CUTOFF ** 2)
    if own_slice[0] < own_slice[1]:
        mask[own_slice[0]:own_slice[1]] = False
    if not mask.any():
        return 0.0
    q = sources[mask, 3]
    r2m = r2[mask]
    phi = COULOMB_KCAL * float(np.sum(q / (DIELECTRIC * r2m)))
    return math.tanh(phi / ESP_SQUASH)


def burial_feature(charged_atoms: np.ndarray,
                   opposite_xyz: np.ndarray) -> float:
    """Desolvation proxy: burial of this residue's charged atoms.

    For each charged atom (|q| > CHARGED_Q): weight |q| by the number of
    OPPOSITE-side heavy atoms within BURIAL_RADIUS. Sum over charged atoms,
    tanh-squashed. A charged group buried in the partner's surface costs
    desolvation free energy -- the signal the v10 net-charge term misses.
    """
    if charged_atoms.shape[0] == 0 or opposite_xyz.shape[0] == 0:
        return 0.0
    charged = charged_atoms[np.abs(charged_atoms[:, 3]) > CHARGED_Q]
    if charged.shape[0] == 0:
        return 0.0
    r2 = _r2_matrix(charged, opposite_xyz)          # (n_charged, n_opp)
    contacts = (r2 < BURIAL_RADIUS ** 2).sum(axis=1)
    burial = float(np.sum(np.abs(charged[:, 3]) * contacts))
    return math.tanh(burial / BURIAL_SQUASH)


def edge_coulomb(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Atom-level pair Coulomb features between two residues.

    Returns (coul, coul_mag), both tanh-squashed:
    - coul     = signed 332.06 * sum q_i q_j / (4 r^2)  over pairs < 10 A
    - coul_mag = unsigned version (|q_i q_j|) -- charge-density contact,
                 a proxy for the electrostatic cost of bringing the two
                 side-chain chemistries into contact.
    """
    if a.shape[0] == 0 or b.shape[0] == 0:
        return 0.0, 0.0
    r2 = _r2_matrix(a, b)                            # (n_a, n_b)
    mask = (r2 > 0.0) & (r2 < ATOM_CUTOFF ** 2)
    if not mask.any():
        return 0.0, 0.0
    qq = a[:, 3, None] * b[None, :, 3]
    inv = COULOMB_KCAL / (DIELECTRIC * r2[mask])
    e = float(np.sum(qq[mask] * inv))
    m = float(np.sum(np.abs(qq[mask]) * inv))
    return math.tanh(e / COUL_SQUASH), math.tanh(m / MAG_SQUASH)
