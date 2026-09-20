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
    """CDR3 + peptide -> fixed-size feature vector (physicochemical profile)."""
    f = []
    for seq in (cdr3, peptide):
        vals = np.array([AA_PHYS.get(a, 0.0) for a in seq])
        # fixed-length profile: mean + moments + positional bins
        f += [
            vals.mean(), vals.std(), vals.min(), vals.max(),
            vals[:7].mean() if len(vals) >= 7 else vals.mean(),
            vals[-7:].mean() if len(vals) >= 7 else vals.mean(),
            len(vals),
        ]
    # 3-mer composition of the peptide (standard for epitopes)
    from collections import Counter
    kmers = Counter(peptide[i:i + 3] for i in range(len(peptide) - 2))
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        f.append(kmers.get(aa, 0) / max(len(peptide) - 2, 1))
    return np.asarray(f, dtype=float)


def load(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    cdr3s, peptides, labels = [], [], []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["cdr3_beta"] and row["peptide"] and row["label"]:
                cdr3s.append(row["cdr3_beta"])
                peptides.append(row["peptide"])
                labels.append(int(row["label"]))
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
