"""Active-learning dry-run simulator (KN-7 / KN-7+, proposal Module 2).

Simulates the dry-wet loop ON PUBLIC DATA before any wet-lab cycle:
- Pool: VDJdb pairs in the epitope-split TEST fold (unseen epitopes — the
  regime the model actually faces for proteome screening).
- Oracle: the true labels exist in the data; "wet validation" = revealing them.
- Model: logistic regression + bootstrap ensemble (MC-variance proxy for EIG).
- Rounds: select a batch with epsilon-greedy EIG, reveal labels, retrain.
- Metric: positive-recall vs sampling fraction; compare against random
  sampling; report log2 efficiency at 5% and 10% sampling.

KN-7+ upgrades (2026-09-22):
- --pool-mode prefiltered: mimic the real ACDC library construction — keep
  only pool pairs whose peptide passes an MHCflurry presentation-predictor
  prefilter (per allele). Tests whether AL gains hold on the higher
  positive-density pool the platform will actually screen.
- --ablation: batch-size ablation at fixed total budget (literature: many
  small batches converge faster than few large ones).

Hypothesis under test (proposal): intelligent sampling recovers most of the
positive signal with 5-10% of the experiments.

Usage:
    .venv/bin/python src/meta_acdc/active_learning/simulate_al.py \
        --data data/processed/vdjdb.clean.tsv --rounds 4 --batch 3000
    .venv/bin/python src/meta_acdc/active_learning/simulate_al.py \
        --pool-mode prefiltered --rounds 6 --batch 1500 --seeds 3
    .venv/bin/python src/meta_acdc/active_learning/simulate_al.py --ablation
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
PREFILTER_ALLELE = "HLA-A*02:01"  # fallback allele when VDJdb allele unsupported


def load_rows(path: Path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["cdr3_beta"] and r["peptide"] and r["label"]]
    return rows


def load(path: Path):
    rows = load_rows(path)
    cdr3s = [r["cdr3_beta"] for r in rows]
    peptides = [r["peptide"] for r in rows]
    labels = [int(r["label"]) for r in rows]
    X = np.vstack([featurize(c, p) for c, p in zip(cdr3s, peptides)])
    return X, np.asarray(labels, dtype=float)


def mhcflurry_prefilter(rows, pool_positions: list[int], thresh: float,
                        cache_path: Path | None) -> np.ndarray:
    """MHCflurry presentation score per pool row (KN-7+ ACDC-pool mode).

    Scores unique (peptide, allele) pairs; caches to TSV so reruns are free.
    Alleles unsupported by MHCflurry fall back to PREFILTER_ALLELE (recorded
    in the cache as allele_used). Returns boolean keep-mask over
    pool_positions.
    """
    cache: dict[tuple[str, str], float] = {}
    if cache_path and cache_path.exists():
        with open(cache_path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cache[(r["peptide"], r["allele"])] = float(r["score"])
    todo = {}
    for i in pool_positions:
        pep, allele = rows[i]["peptide"], rows[i]["mhc_allele"]
        if (pep, allele) not in cache and (pep, allele) not in todo:
            todo[(pep, allele)] = True
    if todo:
        from mhcflurry import Class1PresentationPredictor
        predictor = Class1PresentationPredictor.load()
        by_allele: dict[str, list[str]] = {}
        for (pep, allele) in todo:
            by_allele.setdefault(allele, []).append(pep)
        n_fallback = 0
        for allele, peps in by_allele.items():
            use_allele = allele
            try:
                df = predictor.predict(peptides=peps, alleles=[use_allele],
                                       verbose=0)
            except Exception:
                use_allele = PREFILTER_ALLELE
                try:
                    df = predictor.predict(peptides=peps, alleles=[use_allele],
                                           verbose=0)
                    n_fallback += len(peps)
                except Exception:
                    for p in peps:  # unsupported even for fallback: drop
                        cache[(p, allele)] = -1.0
                    continue
            for p, s in zip(peps, df["presentation_score"].to_numpy()):
                cache[(p, allele)] = float(s)
        if n_fallback:
            print(f"mhcflurry: {n_fallback} peptides fell back to "
                  f"{PREFILTER_ALLELE} (original allele unsupported)",
                  flush=True)
        if cache_path:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["peptide", "allele", "score"],
                                  delimiter="\t")
                w.writeheader()
                for (p, a), s in cache.items():
                    w.writerow({"peptide": p, "allele": a, "score": s})
    scores = np.array([cache[(rows[i]["peptide"], rows[i]["mhc_allele"])]
                       for i in pool_positions])
    return scores >= thresh


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
    ap.add_argument("--mode", choices=["eig", "entropy", "variance", "random"],
                    default="eig")
    ap.add_argument("--seeds", type=int, default=1,
                    help="repeats with different seeds (summary over seeds)")
    ap.add_argument("--pool-mode", choices=["raw", "prefiltered"],
                    default="raw",
                    help="raw: full test-fold pool; prefiltered: ACDC-style "
                         "MHCflurry-presentation prefilter (KN-7+)")
    ap.add_argument("--prefilter-thresh", type=float, default=0.5,
                    help="MHCflurry presentation-score keep threshold")
    ap.add_argument("--mhcflurry-cache", type=Path,
                    default=Path("data/processed/mhcflurry_scores.tsv"))
    ap.add_argument("--ablation", action="store_true",
                    help="batch-size ablation at fixed total budget")
    ap.add_argument("--ablation-budget", type=float, default=0.12,
                    help="total sampling budget as pool fraction")
    ap.add_argument("--ablation-batches", type=str,
                    default="3600,1800,900,600,300",
                    help="comma list of batch sizes for the ablation")
    args = ap.parse_args()

    rows = load_rows(args.data)
    cdr3s = [r["cdr3_beta"] for r in rows]
    peptides = [r["peptide"] for r in rows]
    X, y = load(args.data)

    split = group_split(cdr3s, peptides, group_by="epitope", seed=42)
    X_tr, y_tr = X[split.train_idx], y[split.train_idx]
    pool_idx = np.asarray(split.test_idx)
    X_pool, y_pool = X[pool_idx], y[pool_idx]

    if args.pool_mode == "prefiltered":
        keep = mhcflurry_prefilter(rows, pool_idx.tolist(),
                                   args.prefilter_thresh, args.mhcflurry_cache)
        print(f"prefilter (mhcflurry presentation >= {args.prefilter_thresh}): "
              f"kept {keep.sum()}/{len(keep)} pool pairs "
              f"({keep.mean() * 100:.1f}%)", flush=True)
        pool_idx = pool_idx[keep]
        X_pool, y_pool = X[pool_idx], y[pool_idx]

    n_pos_pool = int(y_pool.sum())
    print(f"pool ({args.pool_mode}, unseen epitopes): {len(pool_idx)} pairs, "
          f"{n_pos_pool} positives ({n_pos_pool / len(pool_idx) * 100:.1f}%)",
          flush=True)

    def run_one(seed: int, rounds: int, batch: int, mode: str) -> list[float]:
        """One simulation; returns recall after each round."""
        rng = np.random.RandomState(seed)
        Xt, yt = X_tr.copy(), y_tr.copy()
        unlabeled = np.arange(len(pool_idx))
        labeled_extra: list[int] = []
        recalls = []
        for rnd in range(rounds):
            if mode == "random":
                n = min(batch, len(unlabeled))
                chosen = rng.choice(unlabeled, size=n, replace=False)
            else:
                probs = fit_ensemble(Xt, yt, X_pool[unlabeled], seed=seed * 100 + rnd)
                if mode == "eig":
                    scores = expected_information_gain(probs)
                elif mode == "entropy":
                    from meta_acdc.active_learning.acquisition import predictive_entropy
                    scores = predictive_entropy(probs)
                else:
                    from meta_acdc.active_learning.acquisition import predictive_uncertainty_mc
                    scores = predictive_uncertainty_mc(probs)
                sel = epsilon_greedy_batch(scores, min(batch, len(unlabeled)),
                                           epsilon=0.15, rng=rng)
                chosen = unlabeled[sel.indices]
            labeled_extra.extend(chosen.tolist())
            Xt = np.vstack([Xt, X_pool[chosen]])
            yt = np.concatenate([yt, y_pool[chosen]])
            unlabeled = np.setdiff1d(unlabeled, chosen)
            recalls.append(y_pool[np.asarray(labeled_extra)].sum() / n_pos_pool)
        return recalls

    if args.ablation:
        # batch-size ablation at fixed total budget (KN-7+, lit: many small
        # batches beat few large ones). Reports final recall per batch size.
        budget = int(round(args.ablation_budget * len(pool_idx)))
        batches = [int(b) for b in args.ablation_batches.split(",") if b]
        print(f"\nbatch ablation: budget {budget} samples "
              f"({args.ablation_budget * 100:.0f}% of pool), "
              f"seeds={max(args.seeds, 3)}", flush=True)
        n_seed = max(args.seeds, 3)
        for b in batches:
            rounds_b = max(1, round(budget / b))
            rec = [run_one(args.seed + k, rounds_b, b, args.mode)
                   for k in range(n_seed)]
            final = np.array([r[-1] for r in rec])
            print(f"  batch {b:5d} x {rounds_b:2d} rounds: recall "
                  f"{final.mean() * 100:5.1f}% (+/- {final.std() * 100:.1f})",
                  flush=True)
        return 0

    all_recalls = [run_one(args.seed + k, args.rounds, args.batch, args.mode)
                   for k in range(args.seeds)]
    rec = np.array(all_recalls)  # (seeds, rounds)
    print(f"\nmode={args.mode}, pool={args.pool_mode}, seeds={args.seeds}, "
          f"rounds={args.rounds}, batch={args.batch}", flush=True)
    for rnd in range(args.rounds):
        frac = (rnd + 1) * args.batch / len(pool_idx)
        mean, std = rec[:, rnd].mean(), rec[:, rnd].std()
        print(f"  {frac * 100:5.1f}% sampling: recall {mean * 100:5.1f}% "
              f"(+/- {std * 100:.1f})", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
