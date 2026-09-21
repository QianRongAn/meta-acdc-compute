"""Split multi-complex PDB files into individual TCR-pMHC complexes.

Some asymmetric units contain 2+ TCR-pMHC complexes (e.g., 1d9k: two TCRs,
two MHCs, two peptides). Grouping rule:
- Pair the TCR-role chains (alpha+beta) into TCR units.
- For each TCR unit, its peptide = the peptide-role chain whose centroid is
  closest to the TCR unit; its MHC unit = the MHC-role chain(s) closest to
  that peptide's centroid.
- Emit one residue list per (TCR, peptide, MHC) grouping.

Usage:
    from meta_acdc.structure.split_complexes import split_complex
    sub_complexes = split_complex(residues)  # list[list[Residue]]
"""

from __future__ import annotations

from meta_acdc.structure.graph import Residue, classify_chains, TCR, MHC, PEPTIDE, B2M


def _centroid(rs: list[Residue]) -> tuple[float, float, float]:
    n = len(rs)
    return (sum(r.x for r in rs) / n, sum(r.y for r in rs) / n,
            sum(r.z for r in rs) / n)


def _dist2(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def split_complex(residues: list[Residue]) -> list[list[Residue]]:
    """Split into sub-complexes; returns [residues...] (one per TCR unit)."""
    roles = classify_chains(residues)
    tcr_chains = sorted(c for c, r in roles.items() if r == TCR)
    pep_chains = sorted(c for c, r in roles.items() if r == PEPTIDE)
    mhc_chains = sorted(c for c, r in roles.items() if r in (MHC, B2M))

    if len(tcr_chains) < 2 or len(pep_chains) < 1:
        return []  # caller keeps its own single-complex handling

    def chain_res(c: str) -> list[Residue]:
        return [r for r in residues if r.chain == c]

    # pair TCR chains: greedy nearest pairing of chain centroids
    centroids = {c: _centroid(chain_res(c)) for c in tcr_chains}
    tcr_units: list[tuple[str, str]] = []
    remaining = list(tcr_chains)
    while len(remaining) >= 2:
        a = remaining.pop(0)
        b = min(remaining, key=lambda c: _dist2(centroids[a], centroids[c]))
        remaining.remove(b)
        tcr_units.append((a, b))

    out: list[list[Residue]] = []
    for a, b in tcr_units:
        unit_cent = _centroid(chain_res(a) + chain_res(b))
        pep = min(pep_chains, key=lambda c: _dist2(unit_cent, centroids.get(c, _centroid(chain_res(c)))))
        pep_cent = _centroid(chain_res(pep))
        mhc = min(mhc_chains, key=lambda c: _dist2(pep_cent, _centroid(chain_res(c))))
        keep_chains = {a, b, pep, mhc}
        sub = [r for r in residues if r.chain in keep_chains]
        out.append(sub)
    return out
