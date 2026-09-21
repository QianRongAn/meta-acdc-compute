"""Hard-negative mining (proposal Module 1, Support Set Construction).

The proposal's triple-criterion filter for "difficult negatives":
  1. Sequence similarity:  Levenshtein distance <= 3 to the target peptide,
     or identical anchor residues.
  2. Structural mimicry:   Pearson correlation > 0.8 of surface electrostatic
     potential between candidate and target pMHC.
  3. Model uncertainty:    prediction probability in [0.4, 0.6] with maximal
     mutual information.

v0 notes:
- Criterion 2 uses a physicochemical-profile proxy (per-position charge +
  hydrophobicity vectors correlated across an aligned sliding window), since
  true electrostatic surfaces require APBS/Poisson-Boltzmann computation.
  The proxy is documented as such; the module exposes the same interface so
  the real computation can be swapped in later (KN-6 upgrade path).
- Criterion 3 reuses acquisition.hard_negative_filter.

Usage example in the meta-learning task builder (models/meta_tasks.py).
"""

from __future__ import annotations

from itertools import product

import numpy as np

from meta_acdc.active_learning.acquisition import hard_negative_filter

AA_PHYS = {
    # [charge, Kyte-Doolittle hydrophobicity]
    "A": [0.0, 1.8], "R": [1.0, -4.5], "N": [0.0, -3.5], "D": [-1.0, -3.5],
    "C": [0.0, 2.5], "Q": [0.0, -3.5], "E": [-1.0, -3.5], "G": [0.0, -0.4],
    "H": [0.1, -3.2], "I": [0.0, 4.5], "L": [0.0, 3.8], "K": [1.0, -3.9],
    "M": [0.0, 1.9], "F": [0.0, 2.8], "P": [0.0, -1.6], "S": [0.0, -0.8],
    "T": [0.0, -0.7], "W": [0.0, -0.9], "Y": [0.0, -1.3], "V": [0.0, 4.2],
    "X": [0.0, 0.0],
}


# --- criterion 1: sequence similarity --------------------------------------

def levenshtein(a: str, b: str) -> int:
    """Classic Levenshtein distance (O(n*m), stdlib-only)."""
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(
                prev[j] + 1,          # deletion
                cur[j - 1] + 1,       # insertion
                prev[j - 1] + (ca != cb),  # substitution
            ))
        prev = cur
    return prev[-1]


def anchor_residues_match(a: str, b: str, n_anchor: int = 2) -> bool:
    """Anchor positions: first 2 + last 1 positions (typical MHC-I anchors
    P1/P2/PC). True if all anchors identical."""
    if len(a) != len(b) or len(a) < n_anchor + 1:
        return False
    a_anchor = [a[0], a[1], a[-1]]
    b_anchor = [b[0], b[1], b[-1]]
    return a_anchor == b_anchor


def sequence_hard_negatives(
    target: str,
    candidates: list[str],
    max_distance: int = 3,
) -> list[str]:
    """Candidates within Levenshtein <= max_distance OR identical anchors."""
    out = []
    for c in candidates:
        if c == target:
            continue
        if levenshtein(target, c) <= max_distance:
            out.append(c)
        elif anchor_residues_match(target, c):
            out.append(c)
    return out


# --- criterion 2: structural mimicry proxy ---------------------------------

def surface_proxy_vector(peptide: str) -> np.ndarray:
    """Per-position [charge, hydrophobicity] profile, length-padded to 15."""
    vec = [AA_PHYS.get(aa, AA_PHYS["X"]) for aa in peptide]
    vec += [[0.0, 0.0]] * (15 - len(vec))
    return np.asarray(vec[:15], dtype=float).ravel()


def structural_mimicry_score(target: str, candidate: str) -> float:
    """Pearson correlation of physicochemical profiles (electrostatics proxy).

    Matches the proposal's >0.8 threshold semantics; replace with true
    surface-potential correlation once APBS structures are computed.
    """
    a = surface_proxy_vector(target)
    b = surface_proxy_vector(candidate)
    a_c = a - a.mean()
    b_c = b - b.mean()
    denom = np.sqrt((a_c ** 2).sum() * (b_c ** 2).sum())
    if denom == 0:
        return 0.0
    return float((a_c * b_c).sum() / denom)


def structural_mimics(
    target: str,
    candidates: list[str],
    threshold: float = 0.8,
) -> list[str]:
    """Candidates with proxy surface correlation > threshold."""
    return [c for c in candidates if c != target
            and structural_mimicry_score(target, c) > threshold]


# --- combined triple-criterion filter --------------------------------------

def hard_negative_candidates(
    target: str,
    candidates: list[str],
    mean_probs: np.ndarray | None = None,
    max_distance: int = 3,
    mimicry_threshold: float = 0.8,
    uncertainty_band: tuple[float, float] = (0.4, 0.6),
) -> tuple[list[str], dict]:
    """Union of the three criteria (proposal: dynamically selected 200).

    Returns (selected, stats). mean_probs aligned with candidates; when None,
    the uncertainty criterion is skipped.
    """
    selected: set[str] = set()
    seq_hits = sequence_hard_negatives(target, candidates, max_distance)
    selected.update(seq_hits)
    struct_hits = structural_mimics(target, candidates, mimicry_threshold)
    selected.update(struct_hits)
    n_unc = 0
    if mean_probs is not None:
        mask = hard_negative_filter(np.asarray(mean_probs),
                                    uncertainty_band[0], uncertainty_band[1])
        unc_hits = [c for c, m in zip(candidates, mask) if m]
        n_unc = len(unc_hits)
        selected.update(unc_hits)
    selected.discard(target)
    return sorted(selected), {
        "sequence": len(seq_hits),
        "structural": len(struct_hits),
        "uncertainty": n_unc,
        "union": len(selected),
    }


def mutate_all_singles(peptide: str) -> list[str]:
    """All single-position substitutions to the 20 standard amino acids."""
    out = []
    for i, orig in enumerate(peptide):
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            if aa != orig:
                out.append(peptide[:i] + aa + peptide[i + 1:])
    return out
