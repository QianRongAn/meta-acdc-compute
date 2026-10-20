"""Statistical tests for the per-TCR ranking experiment.

When N native controls are scored (one per structure), test whether the
model ranks native peptides significantly better than chance:
- pooled exact test: P(rank <= observed) under random permutation of the
  candidate order, pooled over structures with the same candidate count
- per-structure binomial-style test
- family-stratified Fisher pooling (RECOMMENDED): structures that share a
  receptor family are pooled first (within-family Fisher), and the
  family-level p-values are combined across families. This treats receptors
  from the same family as one unit instead of independent structures.

Currently runnable with any number of native controls; meaningful from
N >= 3-4. Rerun as user submissions arrive.

Usage:
    .venv/bin/python src/meta_acdc/structure/rank_stats.py \
        --ensemble data/processed/prediction_scores_ensemble.tsv \
        --map data/processed/structure_vdjdb_map.tsv \
        --families data/processed/structure_family_map.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

from scipy.stats import chi2 as _chi2_dist      # noqa: E402

chi2_sf = _chi2_dist.sf

# p-values are rank/n with rank >= 1 and n >= 1, so p is never 0; the floor
# only guards against underflow through the double log.
_P_FLOOR = 1e-300


def fisher_pool(pvalues: list[float]) -> tuple[float, float]:
    """Fisher-combine p-values; returns (chi2, combined_p).

    Uses the chi2 *survival function*, not a normal approximation of it.
    With k=2 the normal form is anti-conservative (two rank-1-of-10
    controls give chi2=9.21, df=4: exact P=0.0561, normal approx 0.0327
    — it turns "not significant" into "significant").
    """
    if not pvalues:
        raise ValueError("fisher_pool requires at least one p-value")
    if any(not (0.0 < p <= 1.0) for p in pvalues):
        raise ValueError("p-values must be in (0, 1]")
    chi2 = -2.0 * sum(math.log(max(p, _P_FLOOR)) for p in pvalues)
    df = 2 * len(pvalues)
    return chi2, float(chi2_sf(chi2, df))


def family_stratified_fisher(
    rows: list[tuple[str, str, int, int]],
    family_of: dict[str, str],
) -> tuple[float, float, dict[str, float]]:
    """Two-level Fisher pooling.

    `rows` is (pdb, native_peptide, rank, n) — one per scored structure.
    Structures sharing a `family_id` are pooled first (their p-values are
    Fisher-combined into one family-level p), then the family-level
    p-values are Fisher-combined across families. A family with a single
    structure passes its p-value through unchanged.

    Returns (chi2, stratified_p, family_p_values).
    """
    if not rows:
        raise ValueError("family_stratified_fisher requires scored rows")
    by_family: dict[str, list[float]] = defaultdict(list)
    for pdb, _nat, rank, n in rows:
        by_family[family_of[pdb]].append(rank / n)
    fam_ps: dict[str, float] = {}
    for fam, ps in sorted(by_family.items()):
        if len(ps) == 1:
            fam_ps[fam] = ps[0]
        else:
            fam_ps[fam] = fisher_pool(ps)[1]
    chi2, pooled = fisher_pool(list(fam_ps.values()))
    return chi2, pooled, fam_ps


def load_family_map(path: Path) -> dict[str, str]:
    """TSV with header `pdb<TAB>family`; one family id per PDB."""
    families: dict[str, str] = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = set(reader.fieldnames or [])
        if "pdb" not in cols or not ({"family", "family_id"} & cols):
            raise ValueError("families file needs 'pdb' and 'family' columns")
        fam_col = "family" if "family" in cols else "family_id"
        for row in reader:
            pdb, fam = row["pdb"].strip(), row[fam_col].strip()
            if not pdb or not fam:
                continue
            if pdb in families and families[pdb] != fam:
                raise ValueError(f"conflicting family id for {pdb}")
            families[pdb] = fam
    return families


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ensemble", type=Path,
                    default=Path("data/processed/prediction_scores_ensemble.tsv"))
    ap.add_argument("--map", type=Path,
                    default=Path("data/processed/structure_vdjdb_map.tsv"))
    ap.add_argument("--families", type=Path, default=None,
                    help="TSV: pdb<TAB>family. When given, report the "
                         "family-stratified Fisher pooled test in addition "
                         "to the naive (independent-structure) pooled test. "
                         "Structures that share a receptor family are pooled "
                         "within family first; receptors from the same "
                         "family are never treated as independent units.")
    args = ap.parse_args()

    scores: dict[tuple[str, str], float] = {}
    with open(args.ensemble, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pdb, pep = row["job_id"].split("_", 1)
            scores[(pdb, pep)] = float(row["score"])

    native: dict[str, str] = {}
    with open(args.map, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            native.setdefault(row["pdb"], row["pdb_peptide"])
    # NOTE (2026-09-23): the old hardcoded "A6/B7 share chains" alias
    # (native["1qrn"] = "LLFGYPVYV") has been removed — chain fingerprints of
    # the 1ao7 and 1qrn crystals differ, and the alias was based on a
    # mislabeled pdb attribution. Native peptides now come solely from the
    # structure->VDJdb map.

    by_pdb: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for (pdb, pep), s in scores.items():
        by_pdb[pdb].append((pep, s))

    rows = []
    for pdb, entries in sorted(by_pdb.items()):
        nat = native.get(pdb, "")
        if not nat:
            continue
        entries.sort(key=lambda r: -r[1])
        rank = next((i + 1 for i, (p, _) in enumerate(entries) if p == nat), None)
        if rank is None:
            continue
        n = len(entries)
        rows.append((pdb, nat, rank, n))

    print(f"native controls scored: {len(rows)}")
    for pdb, nat, rank, n in rows:
        p = rank / n
        print(f"  {pdb:8s} native={nat:12s} rank={rank}/{n} "
              f"(one-sided P={p:.3f})")
    if len(rows) >= 2:
        # Pooled by Fisher's method: -2*sum(log p) ~ chi2(2k). The p-value is
        # the chi2 *survival function*, not a normal approximation of it.
        chi2 = -2 * sum(math.log(r / n) for _, _, r, n in rows)
        df = 2 * len(rows)
        pooled_p = float(chi2_sf(chi2, df))
        print(f"\nFisher pooled P (structures treated as independent): "
              f"{pooled_p:.4f} (chi2={chi2:.1f}, df={df})")

        if args.families is not None:
            family_of = load_family_map(args.families)
            missing = {pdb for pdb, *_ in rows if pdb not in family_of}
            if missing:
                raise ValueError(
                    "families file lacks entries for scored structures: "
                    + ", ".join(sorted(missing)))
            s_chi2, s_p, fam_ps = family_stratified_fisher(rows, family_of)
            print(f"family-stratified Fisher pooled P: {s_p:.4f} "
                  f"(chi2={s_chi2:.1f}, df={2 * len(fam_ps)}; "
                  f"families={len(fam_ps)})")
            for fam in sorted(fam_ps):
                members = sum(1 for pdb, *_ in rows if family_of[pdb] == fam)
                print(f"  family {fam}: {members} structure(s), "
                      f"family-level P={fam_ps[fam]:.4f}")
            print("interpretation: the stratified test treats receptors from "
                  "the same family as one unit; the naive test treats every "
                  "structure as independent and is an UPPER BOUND on "
                  "significance when families contain several structures. "
                  "P < 0.05 => model ranks native peptides significantly "
                  "above chance under the corresponding pooling.")
        else:
            print("interpretation: P < 0.05 => model ranks native peptides "
                  "significantly above chance across structures. NOTE: this "
                  "treats the listed structures as independent; receptors "
                  "from the same family are not, so treat this as an upper "
                  "bound on significance — pass --families for the "
                  "family-stratified version.")
    else:
        print("\nneed >= 2 native controls for the pooled test; "
              "awaiting 1qse/1qsf submissions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
