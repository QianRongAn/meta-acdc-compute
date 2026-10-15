"""Split multi-complex PDB entries into individual TCR-pMHC complexes.

Some asymmetric units contain 2+ TCR-pMHC complexes (e.g. 1d9k: two TCRs, two
MHCs, two peptides).

Revision after the 2026-10-15 audit (AUDIT.md #27). The previous rule paired
TCR chains by centroid proximity, then picked ONE MHC-role chain by distance to
the peptide. Three problems:

1. centroid proximity is not contact — two receptors can sit far apart in
   space yet be unrelated, or be a genuine alpha/beta pair;
2. a pMHC needs its heavy chain AND beta-2 microglobulin; keeping one chain
   silently truncated the complex;
3. nothing stopped two sub-complexes from claiming the same peptide, so the
   423 extracted structures were not demonstrably 423 independent complexes.

The assignment is now contact-based, and the split reports its own conflicts:

- a TCR unit is an alpha/beta pair whose residues genuinely contact each other
- a peptide belongs to the unit it contacts most strongly
- the MHC unit is every MHC-role chain contacted by that peptide (heavy chain
  plus beta-2 microglobulin)
- if the assignment is ambiguous, the split is reported as UNRESOLVED instead
  of silently guessing

Usage:
    from meta_acdc.structure.split_complexes import split_complex
    sub_complexes, warnings = split_complex(residues)
"""

from __future__ import annotations

from meta_acdc.structure.graph import Residue, classify_chains, TCR, MHC, PEPTIDE, B2M

CONTACT_CUTOFF = 6.0        # A, any-atom contact between two residues


def _res_contact(rs1: list[Residue], rs2: list[Residue],
                 cutoff: float = CONTACT_CUTOFF) -> int:
    """Number of residue pairs within `cutoff` (A) between two chain sets."""
    hit = 0
    c2 = cutoff * cutoff
    for r1 in rs1:
        for r2 in rs2:
            if ((r1.x - r2.x) ** 2 + (r1.y - r2.y) ** 2
                    + (r1.z - r2.z) ** 2) <= c2:
                hit += 1
                break
    return hit


def _centroid(rs: list[Residue]) -> tuple[float, float, float]:
    n = len(rs)
    return (sum(r.x for r in rs) / n, sum(r.y for r in rs) / n,
            sum(r.z for r in rs) / n)


def _dist2(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def split_complex(residues: list[Residue],
                  min_contacts: int = 3
                  ) -> tuple[list[list[Residue]], list[str]]:
    """Split into sub-complexes.

    Returns (complexes, warnings). `warnings` describes every ambiguous
    assignment; callers must surface them instead of treating the output as
    ground truth. An empty list means "single complex" (caller's own path)
    or "unresolved" (see warnings).
    """
    warnings: list[str] = []
    roles = classify_chains(residues)
    tcr_chains = sorted(c for c, r in roles.items() if r == TCR)
    pep_chains = sorted(c for c, r in roles.items() if r == PEPTIDE)
    mhc_chains = sorted(c for c, r in roles.items() if r in (MHC, B2M))

    if len(tcr_chains) < 2 or len(pep_chains) < 1:
        return [], ["single-complex: caller's own handling"]

    def chain_res(c: str) -> list[Residue]:
        return [r for r in residues if r.chain == c]

    # --- TCR units: contact-based alpha/beta pairing -----------------------
    contacts: dict[tuple[str, str], int] = {}
    for i, a in enumerate(tcr_chains):
        for b in tcr_chains[i + 1:]:
            n = _res_contact(chain_res(a), chain_res(b))
            contacts[(a, b)] = n
            contacts[(b, a)] = n

    tcr_units: list[tuple[str, str]] = []
    used: set[str] = set()
    for (a, b), n in sorted(contacts.items(), key=lambda kv: -kv[1]):
        if a in used or b in used or n < min_contacts:
            continue
        # two chains that contact each other are a candidate receptor; when a
        # chain is contacted by two partners we keep the strongest and log it
        rivals = [c for c in tcr_chains
                  if c not in (a, b) and contacts.get((a, c), 0) >= n]
        if rivals:
            warnings.append(
                f"chain {a} contacts {b} ({n}) and {rivals[0]} "
                f"({contacts[(a, rivals[0])]}) equally; kept {b}")
        tcr_units.append((a, b))
        used.update((a, b))
    for c in tcr_chains:
        if c not in used:
            warnings.append(f"TCR-role chain {c} had no contact partner; "
                            f"left out of every sub-complex")

    # --- peptide and MHC assignment by contact strength --------------------
    out: list[list[Residue]] = []
    assigned_peptides: dict[str, tuple[str, int]] = {}
    for a, b in tcr_units:
        unit_res = chain_res(a) + chain_res(b)
        scored = sorted(
            ((_res_contact(unit_res, chain_res(p)), p) for p in pep_chains),
            reverse=True)
        best_score, pep = scored[0]
        if best_score == 0:
            warnings.append(f"TCR unit {a}/{b} contacts no peptide chain; "
                            f"no complex emitted")
            continue
        if len(scored) > 1 and scored[1][0] > 0:
            warnings.append(
                f"TCR unit {a}/{b} contacts peptides {pep} ({best_score}) and "
                f"{scored[1][1]} ({scored[1][0]}); kept the stronger contact")
        if pep in assigned_peptides:
            other, other_score = assigned_peptides[pep]
            if other_score > best_score:
                warnings.append(
                    f"peptide {pep} also claimed by TCR {other} "
                    f"({other_score} > {best_score}); this complex dropped")
                continue
            warnings.append(
                f"peptide {pep} was claimed by TCR {other} "
                f"({other_score}); {a}/{b} has stronger contact, previous "
                f"complex dropped")
            out = [c for c in out if not _uses_chain(c, pep)]
        assigned_peptides[pep] = (f"{a}/{b}", best_score)

        # pMHC = heavy chain(s) contacted by the peptide + B2M chains
        # contacted by those heavy chains. B2M usually contacts the heavy
        # chain, NOT the peptide directly, so peptide-B2M contact is the
        # wrong test; and an incomplete pMHC (missing heavy or B2M) must be
        # dropped, not emitted with a warning (AUDIT.md #27).
        pep_res = chain_res(pep)
        heavy = [c for c in mhc_chains if roles[c] == MHC
                 and _res_contact(pep_res, chain_res(c)) > 0]
        b2m = [c for c in mhc_chains if roles[c] == B2M
               and any(_res_contact(chain_res(c), chain_res(h)) > 0
                       for h in heavy)]
        mhc_unit = heavy + b2m
        if not heavy or not b2m:
            warnings.append(
                f"peptide {pep}: incomplete pMHC (heavy={heavy}, b2m={b2m}); "
                f"a pMHC needs the heavy chain and beta-2 microglobulin — "
                f"no complex emitted")
            continue
        keep = {a, b, pep} | set(mhc_unit)
        out.append([r for r in residues if r.chain in keep])

    if warnings:
        warnings.insert(0, f"{len(warnings)} assignment warnings")
    return out, warnings


def _uses_chain(sub: list[Residue], chain: str) -> bool:
    return any(r.chain == chain for r in sub)
