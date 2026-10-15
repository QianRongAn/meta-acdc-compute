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


def _kmers(seq: str, k: int = 3) -> set[str]:
    if len(seq) < k:
        return {seq} if seq else set()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def cluster_sequences(seqs: list[str], k: int = 3,
                      threshold: float = 0.6) -> list[int]:
    """Single-linkage sequence clusters by k-mer Jaccard similarity.

    Replaces the MD5 `_similarity_bucket` used until 2026-10-15: hashing does
    not preserve similarity, so three CDR3s differing in their last residue
    landed in buckets 26, 31 and 22 and were split across folds (AUDIT.md #13).
    Returns a cluster id per input sequence.
    """
    n = len(seqs)
    parent = list(range(n))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    sets = [_kmers(s, k) if s else set() for s in seqs]
    index: dict[str, list[int]] = {}
    for i, ks in enumerate(sets):
        for km in ks:
            index.setdefault(km, []).append(i)

    for i, ks in enumerate(sets):
        if not ks:
            continue
        counts: dict[int, int] = {}
        for km in ks:
            for j in index[km]:
                if j != i:
                    counts[j] = counts.get(j, 0) + 1
        for j, shared in counts.items():
            other = sets[j]
            denom = len(ks) + len(other) - shared
            if denom and shared / denom >= threshold:
                union(i, j)

    roots = {}
    out = []
    for i in range(n):
        r = find(i)
        if r not in roots:
            roots[r] = len(roots)
        out.append(roots[r])
    return out


def _cluster_key(cdr3: str, clusters: dict[str, int]) -> str:
    return f"T{clusters.get(cdr3, hash(cdr3) % 10**9)}"


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

    group_by:
      'epitope'    — group by peptide
      'tcr'        — group by CDR3 similarity CLUSTER (not by hash bucket)
      'pair'       — group by (TCR, peptide); only the exact pair is isolated
      'unseen-both'— test TCR clusters and test peptides are both disjoint
                     from train: neither the receptor nor the epitope was seen
    """
    rng = random.Random(seed)
    n = len(cdr3s)
    if group_by == "epitope":
        groups: dict[str, list[int]] = {}
        for i, p in enumerate(peptides):
            groups.setdefault(p, []).append(i)
    elif group_by == "tcr":
        clusters = {c: i for i, c in
                    enumerate(cluster_sequences(sorted(set(cdr3s))))}
        groups = {}
        for i, c in enumerate(cdr3s):
            groups.setdefault(_cluster_key(c, clusters), []).append(i)
    elif group_by == "pair":
        groups = {}
        for i, (c, p) in enumerate(zip(cdr3s, peptides)):
            groups.setdefault(f"{c}|{p}", []).append(i)
    elif group_by == "unseen-both":
        # Build the split from the epitope/TCR-cluster incidence, so the same
        # entity never appears on both sides. The old 'pair' split only
        # isolated the exact pair: train and test still shared every TCR and
        # every epitope.
        clusters = {c: i for i, c in
                    enumerate(cluster_sequences(sorted(set(cdr3s))))}
        by_pep: dict[str, set[str]] = {}
        by_tcr: dict[str, set[str]] = {}
        for c, p in zip(cdr3s, peptides):
            tp = _cluster_key(c, clusters)
            by_pep.setdefault(p, set()).add(tp)
            by_tcr.setdefault(tp, set()).add(p)
        pep_ids = sorted(by_pep)
        rng.shuffle(pep_ids)
        test_peps: set[str] = set()
        train_tcr: set[str] = set()
        target = int(test_frac * len(pep_ids))
        for p in pep_ids:
            tcrs = by_pep[p]
            if len(test_peps) < target and not (tcrs & train_tcr):
                test_peps.add(p)
                train_tcr |= {
                    t for t in by_tcr if p in by_tcr[t] and t not in tcrs}
            else:
                for t in tcrs:
                    train_tcr.add(t)
        groups = {"__test__": [], "__train__": []}
        for i, (c, p) in enumerate(zip(cdr3s, peptides)):
            side = "__test__" if p in test_peps else "__train__"
            groups[side].append(i)
        meta = {
            "n_groups": len(test_peps),
            "test_groups": len(test_peps),
            "test_tcr_clusters": sorted({_cluster_key(c, clusters)
                                         for c, p in zip(cdr3s, peptides)
                                         if p in test_peps}),
            "train_tcr_clusters": sorted(train_tcr),
        }
        return Split(
            train_idx=groups["__train__"], test_idx=groups["__test__"],
            name="unseen-both-split", meta=meta)
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
