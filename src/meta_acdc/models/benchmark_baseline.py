"""Sequence-only baseline benchmark on VDJdb (KN-4 prep).

Purpose: quantify the generalization gap on OUR data before the EGNN arrives.
A logistic-regression baseline on k-mer + physicochemical features, evaluated
under the epitope-split (unseen epitopes) regime. Literature expectation
(02-ml-tcr-prediction.md): sequence-only models reach AUROC ~0.5-0.65 on
unseen epitopes (TITAN 0.62; ERGO TPP-III 0.669). EGNN must beat this.

Metrics: AUPRC primary (class imbalance), AUROC secondary.

Usage:
    .venv/bin/python src/meta_acdc/models/benchmark_baseline.py \
        --data data/processed/vdjdb.clean.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split

AA_PHYS = {
    # Kyte-Doolittle hydrophobicity
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2, "X": 0.0,
}


def featurize(cdr3: str, peptide: str) -> np.ndarray:
    """CDR3 + peptide -> fixed-size feature vector.

    Two blocks:
    * per-chain summaries (mean / std / extremes / terminal means / length)
    * TCR x peptide INTERACTION terms. Without these the model is exactly
      s(T,P) = a(T) + b(P): swapping the TCR shifts every peptide's logit by
      the same amount and cannot change the ranking of candidates for that
      receptor, so it cannot act as a per-receptor off-target surrogate
      (AUDIT.md #12).

    The interactions are products of chain-level physicochemical summaries
    (hydrophobicity, length-normalised bulk via extremes and means) plus a
    position-wise complementarity term.
    """
    f = []
    summaries = []
    for seq in (cdr3, peptide):
        vals = np.array([AA_PHYS.get(a, 0.0) for a in seq])
        summary = [
            vals.mean(), vals.std(), vals.min(), vals.max(),
            vals[:7].mean() if len(vals) >= 7 else vals.mean(),
            vals[-7:].mean() if len(vals) >= 7 else vals.mean(),
            len(vals),
        ]
        summaries.append((vals, summary))
        f += summary

    (tv, ts), (pv, ps) = summaries
    # chain-level products: 3 phys scalars (mean, min, max) x same on the other
    for i in (0, 2, 3):
        for j in (0, 2, 3):
            f.append(ts[i] * ps[j])
    # position-wise complementarity: matched hydrophobic core vs mismatch
    n = min(len(tv), len(pv))
    if n:
        f.append(float(np.mean(tv[:n] * pv[:n])))
        f.append(float(np.mean(np.abs(tv[:n] - pv[:n]))))
        f.append(float(np.sum(tv[:n] > 0) * np.sum(pv[:n] > 0)))
    else:
        f += [0.0, 0.0, 0.0]

    # amino-acid composition of the peptide (20-dim). NOTE: the previous
    # implementation counted 3-letter windows ("GIL") but then looked up
    # single letters ("G"), so all 20 values were always 0 and the baseline
    # was silently feature-starved. Keep the lookup and the count on the same
    # alphabet; this is asserted in tests/test_baseline_features.py.
    from collections import Counter
    counts = Counter(peptide)
    n_pep = max(len(peptide), 1)
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        f.append(counts.get(aa, 0) / n_pep)
    return np.asarray(f, dtype=float)


def load(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    """Read (cdr3_beta, peptide, label) from the unified TSV.

    Refuses to run on a file whose label column is empty. Since 2026-10-15 the
    VDJdb importer no longer writes a label at all: vdjdb.score is evidence
    confidence, not a binding outcome, so turning it into 1/0 manufactured
    185,117 "non-binders" and every AUROC computed on them (AUDIT.md #1).
    """
    cdr3s, peptides, labels = [], [], []
    n_missing = 0
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not (row["cdr3_beta"] and row["peptide"]):
                continue
            if not row["label"]:
                n_missing += 1
                continue
            cdr3s.append(row["cdr3_beta"])
            peptides.append(row["peptide"])
            labels.append(int(row["label"]))
    if n_missing and not cdr3s:
        raise ValueError(
            f"{path}: {n_missing} rows have no binding label. VDJdb's score is "
            f"an evidence-confidence field, so binding labels must come from "
            f"an assay-derived source. Nothing can be scored on this file "
            f"until the labels are rebuilt (see docs/AUDIT.md #1).")
    if n_missing:
        print(f"warning: skipped {n_missing} rows without a binding label",
              file=sys.stderr)
    return cdr3s, peptides, np.asarray(labels, dtype=float)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("data/processed/vdjdb.clean.tsv"))
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    cdr3s, peptides, y = load(args.data)
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    print(f"data: {len(y)} pairs (pos={n_pos}, neg={n_neg}, imbalance 1:{n_neg/max(n_pos,1):.0f})")

    X = np.vstack([featurize(c, p) for c, p in zip(cdr3s, peptides)])
    print(f"features: {X.shape}")

    # regime 1: random split (seen-TCR)
    rng = np.random.RandomState(args.seed)
    idx = rng.permutation(len(y))
    n_test = len(y) // 5
    for name, tr, te in [("random-split", idx[n_test:], idx[:n_test])]:
        clf = LogisticRegression(max_iter=1000, class_weight="balanced")
        clf.fit(X[tr], y[tr])
        prob = clf.predict_proba(X[te])[:, 1]
        print(f"{name}: AUROC={roc_auc_score(y[te], prob):.3f} "
              f"AUPRC={average_precision_score(y[te], prob):.3f}")

    # regime 2: epitope-split (unseen epitopes — the generalization gap)
    from meta_acdc.data.splits import group_split
    split = group_split(cdr3s, peptides, group_by="epitope", seed=args.seed)
    clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    clf.fit(X[split.train_idx], y[split.train_idx])
    prob = clf.predict_proba(X[split.test_idx])[:, 1]
    print(f"epitope-split: AUROC={roc_auc_score(y[split.test_idx], prob):.3f} "
          f"AUPRC={average_precision_score(y[split.test_idx], prob):.3f} "
          f"(test={len(split.test_idx)})")
    return 0


if __name__ == "__main__":
    main()
