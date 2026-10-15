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
    # cluster map: {cdr3 -> cluster id} built from the SORTED unique set so
    # the ids line up with cluster_sequences()'s return order. The previous
    # `{c: i for i, c in enumerate(cluster_sequences(...))}` iterated the
    # CLUSTER ID list, produced an identity map on ints, and every CDR3 then
    # fell back to a random hash bucket — similar sequences were split again
    # (AUDIT.md #4/#13).
    def build_clusters() -> dict[str, int]:
        uniq = sorted(set(cdr3s))
        return {c: cid for c, cid in zip(uniq, cluster_sequences(uniq))}

    if group_by == "epitope":
        groups: dict[str, list[int]] = {}
        for i, p in enumerate(peptides):
            groups.setdefault(p, []).append(i)
    elif group_by == "tcr":
        clusters = build_clusters()
        groups = {}
        for i, c in enumerate(cdr3s):
            groups.setdefault(_cluster_key(c, clusters), []).append(i)
    elif group_by == "pair":
        groups = {}
        for i, (c, p) in enumerate(zip(cdr3s, peptides)):
            groups.setdefault(f"{c}|{p}", []).append(i)
    elif group_by == "unseen-both":
        # Test peptides are chosen as whole connected components of the
        # peptide × TCR-cluster incidence graph: neither a TCR cluster nor an
        # epitope then appears on both sides. The old greedy-by-peptide
        # version left the exclusion set empty, so one TCR still appeared in
        # both train and test (AUDIT.md #13).
        clusters = build_clusters()
        by_pep: dict[str, set[str]] = {}
        by_tcr: dict[str, set[str]] = {}
        for c, p in zip(cdr3s, peptides):
            tp = _cluster_key(c, clusters)
            by_pep.setdefault(p, set()).add(tp)
            by_tcr.setdefault(tp, set()).add(p)
        parent: dict[str, str] = {}

        def find(x: str) -> str:
            parent.setdefault(x, x)
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a: str, b: str) -> None:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        for p, ts in by_pep.items():
            for t in ts:
                for q in by_tcr[t]:
                    union(("p", p), ("p", q))
        comps: dict[str, set[str]] = {}
        for p in by_pep:
            comps.setdefault(find(("p", p)), set()).add(p)
        comp_list = sorted(comps.values(), key=len)
        rng.shuffle(comp_list)
        target = max(1, int(test_frac * len(by_pep)))
        test_peps: set[str] = set()
        for comp in comp_list:
            if len(test_peps) >= target:
                break
            test_peps |= comp
        if not test_peps:
            test_peps = set(comp_list[0])
        if len(test_peps) == len(by_pep) and len(comp_list) > 1:
            # everything landed in test: move the smallest component back
            test_peps -= min(comp_list, key=len)
        groups = {"__test__": [], "__train__": []}
        for i, p in enumerate(peptides):
            groups["__test__" if p in test_peps else "__train__"].append(i)
        test_tcrs = sorted({t for p in test_peps for t in by_pep[p]})
        train_tcrs = sorted({t for p in by_pep if p not in test_peps
                             for t in by_pep[p]})
        overlap = set(test_tcrs) & set(train_tcrs)
        if overlap:
            raise RuntimeError(
                f"unseen-both split: {len(overlap)} TCR clusters on both "
                f"sides — incidence closure failed")
        meta = {
            "n_groups": len(test_peps),
            "test_groups": len(test_peps),
            "test_tcr_clusters": test_tcrs,
            "train_tcr_clusters": train_tcrs,
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
