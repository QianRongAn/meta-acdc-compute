"""Rashomon (cross-instance stability) analysis for FOLLOWUP-AUROC-RASHOMON-20261021.

Reads the per-(fold, seed) artifacts produced by
meta_acdc.models.train_fingerprint_v2 (embeddings, similarity matrices,
per-query rankings, metrics fragments) and computes the pre-registered
Rashomon metrics (docs/fingerprint_retrieval_20261021/
FOLLOWUP_AUROC_RASHOMON_20261021.md):

  1. performance band           AUROC / rank-1 spread across instances
  2. ranking stability          per-query Spearman rho of gallery rankings
                                across instance pairs (honest rebuild of the
                                withdrawn 0.644-style metric)
  3. retrieval agreement        rank-1 hit agreement, top-5/top-10 Jaccard
  4. pair flip rate             same/diff decisions at instance-median
                                thresholds; fraction flipping across pairs
  5. embedding agreement        same-graph cosine across instances vs the
                                cross-graph baseline
  6. pooled AUROC / rank-1      instance-averaged similarities (descriptive)

Reconstructs historical summary in memory only. Writes corrected retrospective
analysis to a fresh output directory, never to the source artifact directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import torch
from scipy.stats import spearmanr

from meta_acdc.models.train_fingerprint import auroc

SEEDS = [20261021, 1101, 4242]
FOLDS = [0, 1, 2, 3, 4]
SEED_PAIRS = list(combinations(SEEDS, 2))


def spearman_rho(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2 or np.all(a == a[0]) or np.all(b == b[0]):
        return float("nan")
    return float(spearmanr(a, b).statistic)


def jaccard(a: list[int], b: list[int]) -> float:
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return float("nan")
    return len(sa & sb) / len(sa | sb)


def flip_rate(sims: dict[int, np.ndarray], pair_population: list[tuple[int, int]],
              pairs: list[tuple[int, int]]) -> dict:
    """Fraction of (query, other) pairs whose same/diff decision (similarity
    above the instance-median over the pair population) flips between two
    instances."""
    n_pairs_total = len(pair_population)
    if n_pairs_total == 0:
        return {"mean": None, "n": 0}
    thresholds = {}
    for s, S in sims.items():
        vals = np.array([S[q, j] for q, j in pair_population])
        thresholds[s] = float(np.median(vals))
    flips = []
    for (a, b) in pairs:
        Sa, Sb = sims[a], sims[b]
        ta, tb = thresholds[a], thresholds[b]
        n_flip = 0
        for q, j in pair_population:
            if (Sa[q, j] > ta) != (Sb[q, j] > tb):
                n_flip += 1
        flips.append(n_flip / n_pairs_total)
    return {"mean": float(np.mean(flips)),
            "per_pair": {f"{a}_{b}": v for (a, b), v in zip(pairs, flips)},
            "n": n_pairs_total}


def bootstrap_ci(vals: np.ndarray, groups: list[str], n: int = 2000,
                 seed: int = 20261023) -> dict:
    """Average seed-pair values within query, then query means within family.

    Families, rather than repeated query/seed pairs, are the sampling units.
    Uncertainty is conditional on the cached models and fixed gallery.
    """
    values = np.asarray(vals, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if values.ndim != 2 or len(groups) != len(values):
        raise ValueError("one family ID required per query row")
    finite = np.isfinite(values)
    counts = finite.sum(axis=1)
    query_means = np.divide(np.where(finite, values, 0).sum(axis=1),
                            np.maximum(counts, 1)).astype(float)
    query_means[counts == 0] = np.nan
    labels = sorted(set(groups))
    family_means = []
    for label in labels:
        selected = query_means[np.asarray(groups) == label]
        selected = selected[np.isfinite(selected)]
        if len(selected):
            family_means.append(float(selected.mean()))
    out = {"mean": float(np.mean(family_means)) if family_means else None,
           "ci95": [None, None], "n": int(np.count_nonzero(counts)),
           "n_query_rows": len(values), "n_family_groups": len(labels),
           "n_defined_family_groups": len(family_means),
           "n_undefined_values": int(np.size(values) - finite.sum()),
           "seed_pair_columns": values.shape[1],
           "query_weighted_mean": float(np.mean(query_means[np.isfinite(query_means)]))
                                  if np.any(np.isfinite(query_means)) else None,
           "scope": "family-equal, conditional on fixed cached models and gallery",
           "status": "insufficient_family_groups"}
    if len(family_means) >= 2:
        rng = np.random.default_rng(seed)
        a = np.asarray(family_means)
        means = a[rng.integers(0, len(a), (n, len(a)))].mean(axis=1)
        out["ci95"] = [float(x) for x in np.quantile(means, [0.025, 0.975])]
        out["status"] = "conditional_cluster_bootstrap"
    return out


def fold_rashomon(fold: int, test_idx: list[int], records: list[dict],
                  sims: dict[int, np.ndarray],
                  embs: dict[int, np.ndarray],
                  rankings: dict[int, dict],
                  metrics: dict[int, dict]) -> dict:
    n = len(records)
    groups = [records[q]["family_id"] for q in test_idx]
    train_idx = [i for i in range(n) if i not in set(test_idx)]
    if not train_idx:
        raise ValueError("training-only alignment requires a nonempty training partition")
    gallery_of = {q: [j for j in range(n) if j != q] for q in test_idx}
    pair_population = [(q, j) for q in test_idx for j in gallery_of[q]]

    # --- 1. performance band ------------------------------------------------
    aurocs = [metrics[s]["retrieval"]["auroc"] for s in SEEDS]
    rank1s = [metrics[s]["retrieval"]["rank1"] for s in SEEDS]
    band = {"auroc": {"min": float(np.min(aurocs)),
                      "max": float(np.max(aurocs)),
                      "spread": float(np.max(aurocs) - np.min(aurocs)),
                      "per_seed": {str(s): a for s, a in zip(SEEDS, aurocs)}},
            "rank1": {"min": float(np.min(rank1s)),
                      "max": float(np.max(rank1s)),
                      "spread": float(np.max(rank1s) - np.min(rank1s)),
                      "per_seed": {str(s): r for s, r in zip(SEEDS, rank1s)}}}

    # --- 2. ranking stability (Spearman over gallery, per query) ------------
    rho_vals, rho_by_pair = [], {}
    for (a, b) in SEED_PAIRS:
        pair_rhos = []
        for q in test_idx:
            g = gallery_of[q]
            pair_rhos.append(spearman_rho(sims[a][q, g], sims[b][q, g]))
        rho_by_pair[f"{a}_{b}"] = bootstrap_ci(np.array(pair_rhos), groups)
        rho_vals.extend(pair_rhos)
    stability = {"pooled": bootstrap_ci(np.array(rho_vals).reshape(len(SEED_PAIRS), -1).T, groups),
                 "per_pair": rho_by_pair}

    # --- 3. retrieval agreement ---------------------------------------------
    hit = {s: {} for s in SEEDS}
    for s in SEEDS:
        for q in test_idx:
            top1 = rankings[s][str(q)]["gallery_order"][0]
            hit[s][q] = records[top1]["receptor_id"] == records[q]["receptor_id"]
    agree_rate, jac5, jac10 = [], [], []
    for (a, b) in SEED_PAIRS:
        agree_rate.extend([1.0 if hit[a][q] == hit[b][q] else 0.0
                           for q in test_idx])
        for q in test_idx:
            ga = rankings[a][str(q)]["gallery_order"]
            gb = rankings[b][str(q)]["gallery_order"]
            jac5.append(jaccard(ga[:5], gb[:5]))
            jac10.append(jaccard(ga[:10], gb[:10]))
    retrieval = {
        "rank1_agreement": bootstrap_ci(np.array(agree_rate).reshape(len(SEED_PAIRS), -1).T, groups),
        "rank1_agreement_definition": "binary hit-status agreement (includes both wrong); not identical retrieved identity",
        "rank1_hit_rate": {str(s): float(np.mean(list(hit[s].values())))
                           for s in SEEDS},
        "top5_jaccard": bootstrap_ci(np.array(jac5).reshape(len(SEED_PAIRS), -1).T, groups),
        "top10_jaccard": bootstrap_ci(np.array(jac10).reshape(len(SEED_PAIRS), -1).T, groups),
    }

    # --- 4. pair flip rate --------------------------------------------------
    flips_all = flip_rate(sims, pair_population, SEED_PAIRS)
    samefam_pop = [(q, j) for q, j in pair_population
                   if records[j]["family_id"] == records[q]["family_id"]]
    flips_samefam = flip_rate(sims, samefam_pop, SEED_PAIRS)

    # --- 5. embedding agreement ---------------------------------------------
    agree_graph, baseline = [], []
    for i in test_idx:
        for (a, b) in SEED_PAIRS:
            agree_graph.append(float(embs[a][i] @ embs[b][i]))
            others = embs[a][[j for j in range(n) if j != i]]
            baseline.append(float((embs[a][i][None, :] @ others.T).mean()))
    emb_agree = {
        "same_graph_across_instances": bootstrap_ci(np.array(agree_graph).reshape(len(test_idx), -1), groups),
        "cross_graph_within_instance": bootstrap_ci(np.array(baseline).reshape(len(test_idx), -1), groups),
        "caveat": ("raw cross-instance cosine; the contrastive objective is "
                   "invariant to global rotations of the embedding space, so "
                   "this metric is dominated by arbitrary per-instance "
                   "rotation, not semantic disagreement (see aligned "
                   "supplement below)"),
    }

    # --- 5b. supplement: Procrustes-aligned cross-instance agreement ---------
    # Align each instance's embedding matrix to the reference instance by the
    # orthogonal Procrustes rotation; then recompute same-graph agreement and
    # instance-pooled discrimination. Documented as a supplement to the
    # pre-registered raw metric (which remains reported as-is above).
    ref = embs[SEEDS[0]]
    aligned = {SEEDS[0]: ref}
    for s in SEEDS[1:]:
        u, _, vt = np.linalg.svd(embs[s][train_idx].T @ ref[train_idx], full_matrices=False)
        rot = u @ vt
        aligned[s] = embs[s] @ rot
    agree_aligned, baseline_aligned = [], []
    for i in test_idx:
        for (a, b) in SEED_PAIRS:
            agree_aligned.append(float(aligned[a][i] @ aligned[b][i]))
            others = aligned[a][[j for j in range(n) if j != i]]
            baseline_aligned.append(
                float((aligned[a][i][None, :] @ others.T).mean()))
    S_aligned = {s: aligned[s] @ aligned[s].T for s in SEEDS}
    pair_scores_al, pair_labels_al = [], []
    for q, j in pair_population:
        pair_scores_al.append(float(np.mean([S_aligned[s][q, j]
                                             for s in SEEDS])))
        pair_labels_al.append(1.0 if records[j]["receptor_id"]
                              == records[q]["receptor_id"] else 0.0)
    pooled_aligned_auroc = auroc(np.array(pair_scores_al),
                                 np.array(pair_labels_al))
    aligned_supplement = {
        "method": ("orthogonal Procrustes fitted on training graphs only; "
                   "agreement evaluated on held-out queries; "
                   "comparing; removes the arbitrary global rotation degree "
                   "of freedom of the contrastive objective"),
        "same_graph_across_instances": bootstrap_ci(
            np.array(agree_aligned).reshape(len(test_idx), -1), groups),
        "cross_graph_within_instance": bootstrap_ci(
            np.array(baseline_aligned).reshape(len(test_idx), -1), groups),
        "pooled_auroc": float(pooled_aligned_auroc),
        "maximum_cosine_change_float64": max(float(np.max(np.abs(
            S_aligned[s] - embs[s] @ embs[s].T))) for s in SEEDS),
        "precision_note": "cached float32 similarities and recomputed float64 dots can differ at ties; no gain from rotation",
        "fit_graph_indices": train_idx,
        "evaluation_graph_indices": test_idx,
        "orthogonal_alignment_changes_within_instance_cosines": False,
    }

    # --- 6. pooled similarities (instance-mean) -----------------------------
    S_pool = np.mean([sims[s] for s in SEEDS], axis=0)
    pair_scores, pair_labels = [], []
    for q, j in pair_population:
        pair_scores.append(float(S_pool[q, j]))
        pair_labels.append(1.0 if records[j]["receptor_id"]
                           == records[q]["receptor_id"] else 0.0)
    pooled_auroc = auroc(np.array(pair_scores), np.array(pair_labels))
    pooled_rank1s = []
    for q in test_idx:
        g = gallery_of[q]
        order = sorted(g, key=lambda j: -S_pool[q, j])
        rank = next((k + 1 for k, j in enumerate(order)
                     if records[j]["receptor_id"] == records[q]["receptor_id"]),
                    None)
        if rank is not None:
            pooled_rank1s.append(1.0 if rank == 1 else 0.0)
    pooled = {"auroc": float(pooled_auroc),
              "rank1": float(np.mean(pooled_rank1s)) if pooled_rank1s else None}

    return {
        "fold": fold,
        "n_test": len(test_idx),
        "query_statistics": {
            "query_indices": test_idx, "family_ids": groups,
            "ranking_rhos": np.array(rho_vals).reshape(len(SEED_PAIRS), -1).T.tolist(),
        },
        "n_pairs": len(pair_population),
        "performance_band": band,
        "ranking_stability": stability,
        "retrieval_agreement": retrieval,
        "pair_flip_rate": flips_all,
        "pair_flip_rate_same_family": flips_samefam,
        "embedding_agreement": emb_agree,
        "embedding_agreement_aligned": aligned_supplement,
        "pooled": pooled,
    }


def pool_folds(fold_results: list[dict]) -> dict:
    def gather(key_path):
        vals = []
        for fr in fold_results:
            d = fr
            for k in key_path:
                d = d[k]
            if isinstance(d, dict) and "mean" in d and d["mean"] is not None:
                vals.append(d["mean"])
            elif isinstance(d, (int, float)) and d is not None:
                vals.append(d)
        return float(np.mean(vals)) if vals else None

    return {
        "auroc_band_spread": gather(["performance_band", "auroc", "spread"]),
        "auroc_band_min": gather(["performance_band", "auroc", "min"]),
        "auroc_band_max": gather(["performance_band", "auroc", "max"]),
        "rank1_band_spread": gather(["performance_band", "rank1", "spread"]),
        "ranking_stability_rho": gather(["ranking_stability", "pooled"]),
        "rank1_agreement": gather(["retrieval_agreement", "rank1_agreement"]),
        "top5_jaccard": gather(["retrieval_agreement", "top5_jaccard"]),
        "top10_jaccard": gather(["retrieval_agreement", "top10_jaccard"]),
        "pair_flip_rate": gather(["pair_flip_rate"]),
        "pair_flip_rate_same_family": gather(["pair_flip_rate_same_family"]),
        "embedding_agreement_same_graph": gather(
            ["embedding_agreement", "same_graph_across_instances"]),
        "embedding_agreement_baseline": gather(
            ["embedding_agreement", "cross_graph_within_instance"]),
        "embedding_agreement_aligned_same_graph": gather(
            ["embedding_agreement_aligned", "same_graph_across_instances"]),
        "embedding_agreement_aligned_baseline": gather(
            ["embedding_agreement_aligned", "cross_graph_within_instance"]),
        "pooled_aligned_auroc": gather(
            ["embedding_agreement_aligned", "pooled_auroc"]),
        "pooled_auroc": gather(["pooled", "auroc"]),
        "pooled_rank1": gather(["pooled", "rank1"]),
    }


def rebuild_main_results(followup_dir: Path, config: dict) -> dict:
    """Rebuild followup_results.json (main result = seed 20261021) from
    per-(fold, seed) metric fragments. Idempotent."""
    fold_results = []
    for fold in FOLDS:
        fp = followup_dir / f"metrics_fold{fold}_seed{SEEDS[0]}.json"
        if fp.exists():
            fold_results.append(json.loads(fp.read_text()))
    fold_results.sort(key=lambda d: d["fold"])
    if not fold_results:
        raise SystemExit("no seed-20261021 metric fragments found")
    means = {
        "retrieval_rank1": float(np.nanmean(
            [r["retrieval"]["rank1"] for r in fold_results])),
        "retrieval_mrr": float(np.nanmean(
            [r["retrieval"]["mrr"] for r in fold_results])),
        "retrieval_chance_rank1": float(np.nanmean(
            [r["retrieval"]["chance_rank1"] for r in fold_results])),
        "retrieval_auroc": float(np.nanmean(
            [r["retrieval"]["auroc"] for r in fold_results])),
        "family_top1_k1": float(np.nanmean(
            [r["family"]["1"] for r in fold_results])),
        "domain_accuracy": float(np.nanmean(
            [r["domain"]["accuracy"] for r in fold_results])),
        "median_pairwise_cosine": float(np.nanmean(
            [r["median_pairwise_cosine"] for r in fold_results])),
    }
    agg = dict(config)
    agg["fold_results"] = fold_results
    agg["means"] = means
    return agg


def markdown(fold_results: list[dict], pooled: dict,
             main: dict, out_path: Path) -> None:
    def fmt(x):
        return "undefined" if x is None else f"{x:.4f}"
    lines = ["# Cached fingerprint artifacts: corrected retrospective analysis", "",
             "Historical v2 did not meet the preregistered AUROC >= 0.9 target. "
             "The six-change v1/v2 point difference is unverified; determinism replay failed. "
             "Training also used incorrect batched edge offsets. These scores are not validated gains.", "",
             "Family-cluster intervals average seed pairs within query, then queries within family. "
             "They are conditional on fixed cached models/gallery, not model-population uncertainty. "
             "Procrustes is fitted only on training graphs and evaluated only on test queries. "
             "Orthogonal alignment cannot improve within-instance cosine AUROC.", "",
             "| Fold | Families | Ranking rho | Conditional 95% CI | Held-out aligned cosine |",
             "|---|---:|---:|---|---:|"]
    for fr in fold_results:
        st = fr["ranking_stability"]["pooled"]
        al = fr["embedding_agreement_aligned"]["same_graph_across_instances"]
        lines.append(f"| {fr['fold']} | {st['n_defined_family_groups']} | {fmt(st['mean'])} | "
                     f"{[fmt(x) for x in st['ci95']]} | {fmt(al['mean'])} |")
    lines += ["", "Historical seed-20261021 mean AUROC: " + fmt(main["means"]["retrieval_auroc"]),
              "This is a retained historical observation, not a restored result or new independent test.", "",
              "| Descriptive fold mean | Value |", "|---|---:|"]
    lines += [f"| {k} | {fmt(v)} |" for k,v in pooled.items()]
    out_path.write_text("\n".join(lines) + "\n")


def finite_json(value):
    if isinstance(value, dict):
        return {k: finite_json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [finite_json(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--followup-dir", type=Path,
                    default=Path("docs/fingerprint_retrieval_20261021/followup"))
    ap.add_argument("--dataset", type=Path,
                    default=Path("data/processed/fingerprint_dataset.pt"))
    ap.add_argument("--fold-assignment", type=Path,
                    default=Path("docs/fingerprint_retrieval_20261021/"
                                 "fold_assignment.json"))
    ap.add_argument("--out-dir", type=Path, required=True,
                    help="new output directory; source artifacts are never overwritten")
    args = ap.parse_args(argv)
    if args.out_dir.exists():
        raise FileExistsError(args.out_dir)

    data = torch.load(args.dataset, map_location="cpu", weights_only=False)
    records = data["records"]
    fa = json.loads(args.fold_assignment.read_text())
    fold_of = fa["per_graph_fold"]
    if (len(fold_of) != len(records) or set(fold_of) != set(FOLDS)
            or fa["graph_records"] != [r["complex_id"] for r in records]):
        raise ValueError("fold assignment does not match dataset graph identities")
    for family in {r["family_id"] for r in records}:
        if len({f for f, r in zip(fold_of, records) if r["family_id"] == family}) != 1:
            raise ValueError("family crosses folds")

    config_path = args.followup_dir / "config.json"
    if not config_path.exists():
        raise SystemExit(f"missing {config_path}; run train_fingerprint_v2 first")
    config = json.loads(config_path.read_text())

    input_hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in [args.dataset, args.fold_assignment, config_path]}
    fold_results = []
    for fold in FOLDS:
        test_idx = [i for i, f in enumerate(fold_of) if f == fold]
        sims, embs, rankings, metrics = {}, {}, {}, {}
        for s in SEEDS:
            tag = f"fold{fold}_seed{s}"
            sim_path = args.followup_dir / "model_run_v2" / f"sim_{tag}.npy"
            emb_path = args.followup_dir / "model_run_v2" / f"emb_{tag}.npy"
            rank_path = args.followup_dir / f"rankings_{tag}.json"
            met_path = args.followup_dir / f"metrics_{tag}.json"
            if not (sim_path.exists() and emb_path.exists()
                    and rank_path.exists() and met_path.exists()):
                raise SystemExit(f"missing artifacts for {tag}")
            sims[s] = np.load(sim_path, allow_pickle=False)
            embs[s] = np.load(emb_path, allow_pickle=False).astype(np.float64)
            for path in [sim_path, emb_path, rank_path, met_path]:
                input_hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
            if (sims[s].shape != (len(records), len(records))
                    or embs[s].shape[0] != len(records)
                    or not np.isfinite(sims[s]).all() or not np.isfinite(embs[s]).all()
                    or not np.allclose(sims[s], embs[s] @ embs[s].T, atol=2e-6, rtol=0)):
                raise ValueError("cached similarities/embeddings fail input binding checks")
            rankings[s] = json.loads(rank_path.read_text())
            metrics[s] = json.loads(met_path.read_text())
            if metrics[s].get("fold") != fold or metrics[s].get("seed") != s:
                raise ValueError("metric fragment instance identity mismatch")
            for q in test_idx:
                order = rankings[s][str(q)]["gallery_order"]
                if sorted(order) != [j for j in range(len(records)) if j != q]:
                    raise ValueError("ranking gallery identity mismatch")
                if np.any(np.diff(sims[s][q, order]) > 2e-6):
                    raise ValueError("ranking order conflicts with saved similarities")
        fold_results.append(fold_rashomon(
            fold, test_idx, records, sims, embs, rankings, metrics))

    for source, expected in input_hashes.items():
        if hashlib.sha256(Path(source).read_bytes()).hexdigest() != expected:
            raise RuntimeError("source changed during analysis")
    pooled = pool_folds(fold_results)
    all_q = [row for fr in fold_results for row in fr["query_statistics"]["ranking_rhos"]]
    all_groups = [g for fr in fold_results for g in fr["query_statistics"]["family_ids"]]
    overall_ranking = bootstrap_ci(np.asarray(all_q), all_groups)
    main = rebuild_main_results(args.followup_dir, config)
    out = {
        "task_id": "FOLLOWUP-AUROC-RASHOMON-20261021",
        "seeds": SEEDS,
        "seed_pairs": [[a, b] for a, b in SEED_PAIRS],
        "fold_results": fold_results,
        "pooled": pooled,
        "overall_family_equal_ranking": overall_ranking,
        "main_results": main,
        "source_sha256": input_hashes,
        "input_hashes_unchanged": True,
        "analysis_scope": "retrospective diagnostics of defective cached training; no new acceptance",
    }
    args.out_dir.mkdir(parents=True)
    (args.out_dir / "rashomon_results.json").write_text(
        json.dumps(finite_json(out), indent=2, allow_nan=False) + "\n")
    markdown(fold_results, pooled, main, args.out_dir / "rashomon.md")
    print(json.dumps(pooled, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
