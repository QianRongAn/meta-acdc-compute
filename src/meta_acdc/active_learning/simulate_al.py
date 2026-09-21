"""Active-learning dry-run simulator (KN-7, proposal Module 2).

Simulates the dry-wet loop ON PUBLIC DATA before any wet-lab cycle:
- Pool: VDJdb pairs in the epitope-split TEST fold (unseen epitopes — the
  regime the model actually faces for proteome screening).
- Oracle: the true labels exist in the data; "wet validation" = revealing them.
- Model: logistic regression + bootstrap ensemble (MC-variance proxy for EIG).
- Rounds: select a batch with epsilon-greedy EIG, reveal labels, retrain.
- Metric: positive-recall vs sampling fraction; compare against random
  sampling; report log2 efficiency at 5% and 10% sampling.

Hypothesis under test (proposal): intelligent sampling recovers most of the
positive signal with 5-10% of the experiments.

Usage:
    .venv/bin/python src/meta_acdc/active_learning/simulate_al.py \
        --data data/processed/vdjdb.clean.tsv --rounds 4 --batch 3000
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from meta_acdc.active_learning.acquisition import (
    expected_information_gain,
    epsilon_greedy_batch,
)
from meta_acdc.data.splits import group_split
from meta_acdc.models.benchmark_baseline import featurize

N_ENSEMBLE = 20


def load(path: Path):
    cdr3s, peptides, labels = [], [], []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["cdr3_beta"] and row["peptide"] and row["label"]:
                cdr3s.append(row["cdr3_beta"])
                peptides.append(row["peptide"])
                labels.append(int(row["label"]))
    X = np.vstack([featurize(c, p) for c, p in zip(cdr3s, peptides)])
    return X, np.asarray(labels, dtype=float)


def fit_ensemble(X_tr, y_tr, X_pool, seed: int) -> np.ndarray:
    """Bootstrap ensemble -> (n_ensemble, n_pool) predicted probabilities."""
    rng = np.random.RandomState(seed)
    probs = []
    for k in range(N_ENSEMBLE):
        idx = rng.randint(0, len(y_tr), len(y_tr))
        clf = LogisticRegression(max_iter=1000, class_weight="balanced")
        clf.fit(X_tr[idx], y_tr[idx])
        probs.append(clf.predict_proba(X_pool)[:, 1])
    return np.vstack(probs)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("data/processed/vdjdb.clean.tsv"))
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--batch", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    X, y = load(args.data)
    cdr3s, peptides = None, None
    with open(args.data, newline="", encoding="utf-8", errors="replace") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["cdr3_beta"] and r["peptide"] and r["label"]]
        cdr3s = [r["cdr3_beta"] for r in rows]
        peptides = [r["peptide"] for r in rows]

    split = group_split(cdr3s, peptides, group_by="epitope", seed=42)
    X_tr, y_tr = X[split.train_idx], y[split.train_idx]
    pool_idx = split.test_idx
    X_pool, y_pool = X[pool_idx], y[pool_idx]
    n_pos_pool = int(y_pool.sum())
    print(f"pool (unseen epitopes): {len(pool_idx)} pairs, {n_pos_pool} positives "
          f"({n_pos_pool / len(pool_idx) * 100:.1f}%)", flush=True)

    # --- random sampling baseline ---
    rng = np.random.RandomState(args.seed)
    order = rng.permutation(len(pool_idx))
    batch_total = args.batch * args.rounds
    rand_recall = [y_pool[order[:b]].sum() / n_pos_pool
                   for b in range(args.batch, batch_total + 1, args.batch)]

    # --- active learning ---
    unlabeled = np.arange(len(pool_idx))
    labeled_extra: list[int] = []
    recall_history = []
    al_seed = 100 + args.seed
    for rnd in range(args.rounds):
        probs = fit_ensemble(X_tr, y_tr, X_pool[unlabeled], seed=al_seed + rnd)
        scores = expected_information_gain(probs)
        sel = epsilon_greedy_batch(scores, min(args.batch, len(unlabeled)),
                                   epsilon=0.15, rng=rng)
        chosen = unlabeled[sel.indices]
        labeled_extra.extend(chosen.tolist())
        # oracle: reveal labels, add to training set
        X_tr = np.vstack([X_tr, X_pool[chosen]])
        y_tr = np.concatenate([y_tr, y_pool[chosen]])
        unlabeled = np.setdiff1d(unlabeled, chosen)
        recall = y_pool[np.asarray(labeled_extra)].sum() / n_pos_pool
        frac = len(labeled_extra) / len(pool_idx)
        recall_history.append((frac, recall))
        print(f"round {rnd + 1}: sampled {len(labeled_extra)} ({frac * 100:.1f}%) "
              f"-> positive recall {recall * 100:.1f}% "
              f"[exploit={sel.n_exploit} explore={sel.n_explore}]", flush=True)

    # --- report ---
    print("\n=== recall vs sampling fraction ===")
    for i, r in enumerate(rand_recall):
        frac = (i + 1) * args.batch / len(pool_idx)
        al_recall = recall_history[i][1]
        print(f"  {frac * 100:5.1f}%:  random {r * 100:5.1f}%   "
              f"active {al_recall * 100:5.1f}%   "
              f"gain {al_recall / max(r, 1e-9):.2f}x", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
