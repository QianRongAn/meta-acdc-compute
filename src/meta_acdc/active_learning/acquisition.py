"""Active-learning acquisition functions for the dry-wet loop (proposal Module 2).

Design decisions (lit reviews 01-03):
- Pooled ACDC experiments run in batches of ~50k peptides => we need BATCH
  acquisition, not pure sequential (AL literature: more batches with fewer
  samples each converges faster, but oligo-pool economics fix the batch size).
- Acquisition: Expected Information Gain (EIG) approximated via MC dropout or
  ensemble disagreement + epsilon-greedy exploration (proposal: 80-90%
  exploitation / 10-20% exploration).
- Hard-negative criterion from the proposal: prediction probability in
  [0.4, 0.6] with maximal mutual information.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass
class AcquisitionResult:
    indices: list[int]        # selected candidate indices (into the pool)
    scores: np.ndarray        # acquisition score per selected candidate
    n_exploit: int
    n_explore: int


def _entropy(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -(p * np.log(p) + (1 - p) * np.log(1 - p))


def predictive_uncertainty_mc(probs: np.ndarray, n_samples: int = 20) -> np.ndarray:
    """MC-dropout-style uncertainty: variance of sampled probabilities.

    probs: (n_samples, n_candidates) array of predicted P(binding).
    Returns per-candidate variance (proxy for EIG under symmetric noise).
    """
    return probs.var(axis=0)


def predictive_entropy(probs: np.ndarray) -> np.ndarray:
    """Entropy of the mean prediction (pure uncertainty sampling)."""
    return _entropy(probs.mean(axis=0))


def expected_information_gain(probs: np.ndarray) -> np.ndarray:
    """BALD-style EIG: H(E[p]) - E[H(p)]."""
    total = _entropy(probs.mean(axis=0))
    per_sample = _entropy(probs).mean(axis=0)
    return total - per_sample


def epsilon_greedy_batch(
    scores: np.ndarray,
    batch_size: int,
    epsilon: float = 0.15,
    rng: np.random.Generator | None = None,
) -> AcquisitionResult:
    """Select a batch: (1-eps) by score, eps uniformly at random (proposal Risk C)."""
    rng = rng or np.random.default_rng()
    n_pool = len(scores)
    batch_size = min(batch_size, n_pool)
    n_explore = int(round(epsilon * batch_size))
    n_exploit = batch_size - n_explore

    exploit_idx = np.argsort(-scores)[:n_exploit]
    remaining = np.setdiff1d(np.arange(n_pool), exploit_idx, assume_unique=False)
    explore_idx = rng.choice(remaining, size=min(n_explore, len(remaining)), replace=False)
    indices = np.concatenate([exploit_idx, explore_idx])
    return AcquisitionResult(
        indices=indices.tolist(),
        scores=scores[indices],
        n_exploit=len(exploit_idx),
        n_explore=len(explore_idx),
    )


def hard_negative_filter(
    mean_probs: np.ndarray,
    lo: float = 0.4,
    hi: float = 0.6,
) -> np.ndarray:
    """Proposal criterion: ambiguous predictions near the decision boundary.

    Returns a boolean mask over candidates in the [lo, hi] band.
    """
    return (mean_probs >= lo) & (mean_probs <= hi)


def select_round(
    probs: np.ndarray,
    batch_size: int,
    epsilon: float = 0.15,
    mode: str = "eig",
    rng: np.random.Generator | None = None,
) -> AcquisitionResult:
    """One dry-wet round selection from an MC-ensembled probability array.

    probs: (n_mc, n_pool). mode: 'eig' | 'entropy' | 'variance'.
    """
    if mode == "eig":
        scores = expected_information_gain(probs)
    elif mode == "entropy":
        scores = predictive_entropy(probs)
    elif mode == "variance":
        scores = predictive_uncertainty_mc(probs)
    else:
        raise ValueError(f"unknown mode: {mode}")
    return epsilon_greedy_batch(scores, batch_size, epsilon=epsilon, rng=rng)


def log2_fold_efficiency(covered: float, sampled: float) -> float:
    """Convenience: report coverage-per-sample in log2 scale."""
    return math.log2(max(covered, 1e-9) / max(sampled, 1e-9))
