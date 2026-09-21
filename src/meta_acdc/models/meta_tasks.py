"""Few-shot task construction for meta-learning (proposal: per-TCR/epitope Tasks).

v0 task definition (PanPep-style, lit review 02 section 5.5):
  A task = "does this TCR bind epitope E?".
  - Support set: up to 50 positive TCRs (VDJdb-validated binders of E) +
    4x hard negatives (TCRs whose validated epitopes are sequence/structural
    neighbors of E — the proposal's 1:4 ratio; random negatives as fallback).
  - Query set: held-out TCRs (binders of E + non-binders), used for the
    few-shot AUROC evaluation.

Hard negatives come from models.hard_negatives (triple-criterion filter on
the EPITOPE level: TCRs bound to mimics of E are the difficult distractors).
"""

from __future__ import annotations

import random
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from meta_acdc.models.hard_negatives import (
    levenshtein,
    sequence_hard_negatives,
    structural_mimics,
)


def load_vdjdb(path: Path) -> tuple[defaultdict, dict]:
    """Return epitope -> list[(cdr3, label)] and allele info."""
    import csv

    epitope_tcrs: defaultdict[str, list[str]] = defaultdict(list)
    tcr_epitopes: dict[str, str] = {}  # CDR3 -> its (first) validated epitope
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cdr3, pep, label = row["cdr3_beta"], row["peptide"], row["label"]
            if not cdr3 or not pep:
                continue
            if label == "1":
                epitope_tcrs[pep].append(cdr3)
                tcr_epitopes.setdefault(cdr3, pep)
    return epitope_tcrs, tcr_epitopes


@dataclass
class Task:
    epitope: str
    support_pos: list[str]
    support_neg: list[str]
    query: list[str]
    query_labels: list[int]
    n_way: int = 2
    meta: dict = field(default_factory=dict)


class TaskSampler:
    """Episodic sampler of few-shot tasks from VDJdb."""

    def __init__(
        self,
        epitope_tcrs: defaultdict,
        tcr_epitopes: dict,
        min_pos: int = 10,
        n_support_pos: int = 50,
        neg_ratio: float = 4.0,
        seed: int = 42,
        use_hard_negatives: bool = True,
    ):
        self.rng = random.Random(seed)
        self.epitopes = sorted(
            e for e, tcrs in epitope_tcrs.items() if len(set(tcrs)) >= min_pos
        )
        self.epitope_tcrs = {e: list(set(tcrs)) for e, tcrs in epitope_tcrs.items()}
        self.tcr_epitopes = tcr_epitopes
        self.n_support_pos = n_support_pos
        self.neg_ratio = neg_ratio
        self.use_hard = use_hard_negatives
        # negative TCR pool: CDR3s validated against some OTHER epitope
        self.neg_pool = list(tcr_epitopes.keys())

    def _hard_negative_tcrs(self, epitope: str, k: int) -> list[str]:
        """TCRs validated against epitopes that mimic `epitope` (difficult)."""
        mimic_epitopes = set(sequence_hard_negatives(
            epitope, [e for e in self.epitope_tcrs if e != epitope], max_distance=3
        ))
        mimic_epitopes |= set(structural_mimics(
            epitope, [e for e in self.epitope_tcrs if e != epitope], threshold=0.8
        ))
        tcrs = [t for e in mimic_epitopes for t in self.epitope_tcrs[e]]
        self.rng.shuffle(tcrs)
        return tcrs[:k]

    def sample_task(self, epitope: str | None = None) -> Task:
        if epitope is None:
            epitope = self.rng.choice(self.epitopes)
        pos_all = self.epitope_tcrs[epitope]
        self.rng.shuffle(pos_all)
        n_sup = min(self.n_support_pos, max(2, int(0.6 * len(pos_all))))
        support_pos = pos_all[:n_sup]
        rest_pos = pos_all[n_sup:]

        n_neg_sup = int(round(n_sup * self.neg_ratio))
        neg_sup: list[str] = []
        if self.use_hard:
            hard = self._hard_negative_tcrs(epitope, n_neg_sup)
            neg_sup = [t for t in hard if t not in support_pos]
        if len(neg_sup) < n_neg_sup:
            pool = [t for t in self.neg_pool
                    if t not in support_pos and t not in neg_sup
                    and self.tcr_epitopes.get(t) != epitope]
            self.rng.shuffle(pool)
            neg_sup += pool[: n_neg_sup - len(neg_sup)]

        # query: held-out positives + random negatives (epitope-split regime)
        n_q = min(50, len(rest_pos))
        query = rest_pos[:n_q]
        labels = [1] * n_q
        pool_q = [t for t in self.neg_pool
                  if t not in support_pos and t not in neg_sup and t not in query]
        self.rng.shuffle(pool_q)
        query += pool_q[:n_q]
        labels += [0] * (len(query) - n_q)
        return Task(epitope=epitope, support_pos=support_pos, support_neg=neg_sup,
                    query=query, query_labels=labels,
                    meta={"n_support_pos": len(support_pos),
                          "n_support_neg": len(neg_sup)})


def split_epitopes(sampler: TaskSampler, val_frac: float = 0.2, seed: int = 7):
    rng = random.Random(seed)
    eps = sampler.epitopes[:]
    rng.shuffle(eps)
    n_val = max(1, int(val_frac * len(eps)))
    return set(eps[n_val:]), set(eps[:n_val])
