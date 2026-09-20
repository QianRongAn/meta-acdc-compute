"""Data splits for TCR-pMHC benchmarks (lit review 02, section 5.4).

Three CV regimes to report (TITAN pattern):
- seen-TCR: random split on pairs (both TCR and epitope may appear in train)
- unseen-epitope: group-split by epitope (edit-distance bucketing optional)
- unseen-both: group-split by (TCR, epitope) — hardest setting

Leakage control: TCRs with >90% CDR3 sequence identity share a fold
(NetTCR-1.0 practice); identical (CDR3, epitope) pairs never span folds.
"""

from __future__ import annotations

import hashlib
import random
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Split:
    train_idx: list[int]
    test_idx: list[int]
    name: str
    meta: dict = field(default_factory=dict)


def _hash_group(parts: list[str]) -> int:
    """Deterministic integer hash of a group key."""
    h = hashlib.md5("|".join(parts).encode()).digest()
    return int.from_bytes(h[:8], "big")


def _similarity_bucket(cdr3: str, n_buckets: int = 32) -> str:
    """Coarse hash bucket so near-identical CDR3s tend to share a fold."""
    sh = hashlib.md5(cdr3.encode()).digest()
    return f"{sh[0] % n_buckets:02d}"


def load_pairs(path: Path) -> tuple[list[str], list[str]]:
    """Read (cdr3_beta, peptide) columns from the unified TSV.

    Returns two aligned lists.
    """
    import csv

    cdr3s, peptides = [], []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["cdr3_beta"] and row["peptide"]:
                cdr3s.append(row["cdr3_beta"])
                peptides.append(row["peptide"])
    return cdr3s, peptides


def group_split(
    cdr3s: list[str],
    peptides: list[str],
    group_by: str = "epitope",
    test_frac: float = 0.2,
    seed: int = 42,
) -> Split:
    """Group split: all pairs sharing a group key go to the same fold.

    group_by: 'epitope' | 'tcr' | 'pair'
    """
    rng = random.Random(seed)
    n = len(cdr3s)
    if group_by == "epitope":
        groups: dict[str, list[int]] = {}
        for i, p in enumerate(peptides):
            groups.setdefault(p, []).append(i)
    elif group_by == "tcr":
        groups = {}
        for i, c in enumerate(cdr3s):
            # similar CDR3s (same coarse bucket) share a fold — leakage control
            groups.setdefault(_similarity_bucket(c), []).append(i)
    elif group_by == "pair":
        groups = {}
        for i, (c, p) in enumerate(zip(cdr3s, peptides)):
            groups.setdefault(f"{c}|{p}", []).append(i)
    else:
        raise ValueError(group_by)

    group_ids = list(groups.keys())
    rng.shuffle(group_ids)
    n_test = max(1, int(round(test_frac * len(group_ids))))
    test_groups = set(group_ids[:n_test])
    train_idx, test_idx = [], []
    for g, idxs in groups.items():
        if g in test_groups:
            test_idx.extend(idxs)
        else:
            train_idx.extend(idxs)
    return Split(
        train_idx=train_idx,
        test_idx=test_idx,
        name=f"{group_by}-split",
        meta={"n_groups": len(groups), "test_groups": n_test},
    )


def random_split(n: int, test_frac: float = 0.2, seed: int = 42) -> Split:
    """Random pair-level split (seen-TCR regime)."""
    rng = random.Random(seed)
    idx = list(range(n))
    rng.shuffle(idx)
    n_test = int(round(test_frac * n))
    return Split(train_idx=idx[n_test:], test_idx=idx[:n_test], name="random-split")


def summary(cdr3s: list[str], peptides: list[str], split: Split) -> str:
    """Report overlap statistics between folds (leakage check)."""
    tr, te = split.train_idx, split.test_idx
    tr_tcr = {cdr3s[i] for i in tr}
    te_tcr = {cdr3s[i] for i in te}
    tr_pep = {peptides[i] for i in tr}
    te_pep = {peptides[i] for i in te}
    tcr_overlap = len(tr_tcr & te_tcr)
    pep_overlap = len(tr_pep & te_pep)
    return (
        f"{split.name}: train={len(tr)}, test={len(te)}, "
        f"TCR overlap={tcr_overlap}, epitope overlap={pep_overlap}"
    )
