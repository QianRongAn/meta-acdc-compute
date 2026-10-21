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

Also rebuilds followup_results.json from the seed-20261021 fragments if it is
missing or incomplete. Descriptive reporting only; no pass/fail bar.
"""
from __future__ import annotations

import argparse
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
RNG = np.random.default_rng(20261021)


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


def bootstrap_ci(vals: np.ndarray, n: int = 2000) -> dict:
    vals = np.asarray(vals, dtype=float)
    vals = vals[~np.isnan(vals)]
    if vals.size == 0:
        return {"mean": None, "ci95": [None, None], "n": 0}
    means = [vals[RNG.integers(0, vals.size, vals.size)].mean()
             for _ in range(n)]
    return {"mean": float(vals.mean()),
            "ci95": [float(np.percentile(means, 2.5)),
                     float(np.percentile(means, 97.5))],
            "n": int(vals.size)}


def fold_rashomon(fold: int, test_idx: list[int], records: list[dict],
                  sims: dict[int, np.ndarray],
                  embs: dict[int, np.ndarray],
                  rankings: dict[int, dict],
                  metrics: dict[int, dict]) -> dict:
    n = len(records)
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
        rho_by_pair[f"{a}_{b}"] = bootstrap_ci(np.array(pair_rhos))
        rho_vals.extend(pair_rhos)
    stability = {"pooled": bootstrap_ci(np.array(rho_vals)),
                 "per_pair": rho_by_pair}

    # --- 3. retrieval agreement ---------------------------------------------
    hit = {s: {} for s in SEEDS}
    for s in SEEDS:
        for q in test_idx:
            top1 = rankings[s][str(q)]["gallery_order"][0]
            hit[s][q] = records[top1]["receptor_id"] == records[q]["receptor_id"]
    agree_rate, jac5, jac10 = [], [], []
    for (a, b) in SEED_PAIRS:
        agree = np.mean([1.0 if hit[a][q] == hit[b][q] else 0.0
                         for q in test_idx])
        agree_rate.append(agree)
        for q in test_idx:
            ga = rankings[a][str(q)]["gallery_order"]
            gb = rankings[b][str(q)]["gallery_order"]
            jac5.append(jaccard(ga[:5], gb[:5]))
            jac10.append(jaccard(ga[:10], gb[:10]))
    retrieval = {
        "rank1_agreement": bootstrap_ci(np.array(agree_rate)),
        "rank1_hit_rate": {str(s): float(np.mean(list(hit[s].values())))
                           for s in SEEDS},
        "top5_jaccard": bootstrap_ci(np.array(jac5)),
        "top10_jaccard": bootstrap_ci(np.array(jac10)),
    }

    # --- 4. pair flip rate --------------------------------------------------
    flips_all = flip_rate(sims, pair_population, SEED_PAIRS)
    samefam_pop = [(q, j) for q, j in pair_population
                   if records[j]["family_id"] == records[q]["family_id"]]
    flips_samefam = flip_rate(sims, samefam_pop, SEED_PAIRS)

    # --- 5. embedding agreement ---------------------------------------------
    agree_graph, baseline = [], []
    for i in range(n):
        for (a, b) in SEED_PAIRS:
            agree_graph.append(float(embs[a][i] @ embs[b][i]))
            others = embs[a][[j for j in range(n) if j != i]]
            baseline.append(float((embs[a][i][None, :] @ others.T).mean()))
    emb_agree = {
        "same_graph_across_instances": bootstrap_ci(np.array(agree_graph)),
        "cross_graph_within_instance": bootstrap_ci(np.array(baseline)),
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
        u, _, vt = np.linalg.svd(embs[s].T @ ref, full_matrices=False)
        rot = u @ vt
        aligned[s] = embs[s] @ rot
    agree_aligned, baseline_aligned = [], []
    for i in range(n):
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
        "method": ("orthogonal Procrustes alignment of each instance's full "
                   "embedding matrix to the reference instance before "
                   "comparing; removes the arbitrary global rotation degree "
                   "of freedom of the contrastive objective"),
        "same_graph_across_instances": bootstrap_ci(
            np.array(agree_aligned)),
        "cross_graph_within_instance": bootstrap_ci(
            np.array(baseline_aligned)),
        "pooled_auroc": float(pooled_aligned_auroc),
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
    out = followup_dir / "followup_results.json"
    out.write_text(json.dumps(agg, indent=2, allow_nan=False) + "\n")
    return agg


def markdown(fold_results: list[dict], pooled: dict,
             main: dict, out_path: Path) -> None:
    lines = []
    lines.append("# FOLLOWUP-AUROC-RASHOMON-20261021 结果")
    lines.append("")
    lines.append("预注册：FOLLOWUP_AUROC_RASHOMON_20261021.md；"
                 "逐折/逐查询数据：rashomon_results.json；"
                 "权重/嵌入/相似度矩阵：model_run_v2/。")
    lines.append("")
    lines.append("## 主结果（seed 20261021，5 折，v2 配方）")
    lines.append("")
    lines.append("| 折 | rank-1 | 机会 | 倍数 | MRR | AUROC | 家族 top-1 | "
                 "域诊断 | 两两余弦中位数 |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for fr in main["fold_results"]:
        r, f, d = fr["retrieval"], fr["family"], fr["domain"]
        lines.append(
            f"| {fr['fold']} | {r['rank1']:.4f} | {r['chance_rank1']:.4f} | "
            f"{r['rank1'] / r['chance_rank1']:.1f}× | {r['mrr']:.4f} | "
            f"{r['auroc']:.4f} | {f['1']:.4f} | {d['accuracy']:.4f} | "
            f"{fr['median_pairwise_cosine']:.4f} |")
    m = main["means"]
    lines.append(f"| **均值** | **{m['retrieval_rank1']:.4f}** | "
                 f"**{m['retrieval_chance_rank1']:.4f}** | "
                 f"**{m['retrieval_rank1'] / m['retrieval_chance_rank1']:.1f}×** | "
                 f"**{m['retrieval_mrr']:.4f}** | **{m['retrieval_auroc']:.4f}** | "
                 f"**{m['family_top1_k1']:.4f}** | **{m['domain_accuracy']:.4f}** | "
                 f"**{m['median_pairwise_cosine']:.4f}** |")
    lines.append("")
    lines.append("## Rashomon（每折 3 实例 × 3 种子对，池化）")
    lines.append("")
    lines.append("| 折 | AUROC 带(min–max) | rank-1 带 | 排名 ρ | rank-1 一致率 | "
                 "top-5 J | top-10 J | 判决翻转率 | 同家族翻转率 | 池化 AUROC |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for fr in fold_results:
        b, st, ra, fl, em, po = (fr["performance_band"], fr["ranking_stability"],
                                 fr["retrieval_agreement"], fr["pair_flip_rate"],
                                 fr["embedding_agreement"], fr["pooled"])
        lines.append(
            f"| {fr['fold']} | {b['auroc']['min']:.3f}–{b['auroc']['max']:.3f} | "
            f"{b['rank1']['min']:.3f}–{b['rank1']['max']:.3f} | "
            f"{st['pooled']['mean']:.3f} | {ra['rank1_agreement']['mean']:.3f} | "
            f"{ra['top5_jaccard']['mean']:.3f} | {ra['top10_jaccard']['mean']:.3f} | "
            f"{fl['mean']:.3f} | {fr['pair_flip_rate_same_family']['mean']:.3f} | "
            f"{po['auroc']:.3f} |")
    lines.append("")
    lines.append("| 池化量 | 值 |")
    lines.append("|---|---:|")
    for k, v in pooled.items():
        lines.append(f"| {k} | {v:.4f} |")
    lines.append("")
    lines.append("## 补充（协议修订注记：嵌入一致率的全局旋转自由度）")
    lines.append("")
    lines.append("预注册的\"嵌入一致率\"按跨实例原始余弦计算。对比学习目标对嵌入空间的"
                 "**全局正交旋转不变**：不同实例收敛到互异旋转的等价空间，原始余弦因此"
                 "被旋转自由度支配，不作语义解读（实测原始同图余弦 ≈ 0.00，见上表）。"
                 "本注记补充对齐口径：各实例嵌入矩阵先经正交 Procrustes 对齐到参考实例"
                 "（无监督、不使用标签），再重算同图一致率与池化 AUROC。预注册的原始"
                 "指标照旧保留并报告，不因补充口径改动任何判定。")
    lines.append("")
    lines.append("| 折 | 对齐同图余弦 | 对齐跨图基线 | 对齐池化 AUROC |")
    lines.append("|---|---:|---:|---:|")
    for fr in fold_results:
        al = fr["embedding_agreement_aligned"]
        lines.append(
            f"| {fr['fold']} | {al['same_graph_across_instances']['mean']:.3f} | "
            f"{al['cross_graph_within_instance']['mean']:.3f} | "
            f"{al['pooled_auroc']:.3f} |")
    lines.append("")
    lines.append("## 预注册标准核对")
    lines.append("")
    m = main["means"]
    crit = [
        ("池化同/异 AUROC ≥ 0.9",
         m["retrieval_auroc"] >= 0.9, f"{m['retrieval_auroc']:.4f}"),
        ("检索 rank-1 ≥ 2× 机会",
         m["retrieval_rank1"] >= 2 * m["retrieval_chance_rank1"],
         f"{m['retrieval_rank1'] / m['retrieval_chance_rank1']:.1f}×"),
        ("两两余弦中位数 ≤ 0.95（表征展开）",
         m["median_pairwise_cosine"] <= 0.95,
         f"{m['median_pairwise_cosine']:.4f}"),
        ("域诊断 ≤ 0.6",
         m["domain_accuracy"] <= 0.6, f"{m['domain_accuracy']:.4f}"),
    ]
    for name, ok, value in crit:
        lines.append(f"- [{'x' if ok else ' '}] {name}（实测 {value}）")
    lines.append("- [ ] 同 seed 确定性（scripts/check_replay.py 外部验证，结果见 "
                 "followup/STATUS.md 与 replay/）")
    lines.append("")
    lines.append("以上判定按训练前固定的标准原样给出，不据结果调整。")
    out_path.write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--followup-dir", type=Path,
                    default=Path("docs/fingerprint_retrieval_20261021/followup"))
    ap.add_argument("--dataset", type=Path,
                    default=Path("data/processed/fingerprint_dataset.pt"))
    ap.add_argument("--fold-assignment", type=Path,
                    default=Path("docs/fingerprint_retrieval_20261021/"
                                 "fold_assignment.json"))
    args = ap.parse_args(argv)

    data = torch.load(args.dataset, map_location="cpu", weights_only=False)
    records = data["records"]
    fa = json.loads(args.fold_assignment.read_text())
    fold_of = fa["per_graph_fold"]

    config_path = args.followup_dir / "config.json"
    if not config_path.exists():
        raise SystemExit(f"missing {config_path}; run train_fingerprint_v2 first")
    config = json.loads(config_path.read_text())

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
            sims[s] = np.load(sim_path)
            embs[s] = np.load(emb_path).astype(np.float64)
            rankings[s] = json.loads(rank_path.read_text())
            metrics[s] = json.loads(met_path.read_text())
        fold_results.append(fold_rashomon(
            fold, test_idx, records, sims, embs, rankings, metrics))

    pooled = pool_folds(fold_results)
    main = rebuild_main_results(args.followup_dir, config)
    out = {
        "task_id": "FOLLOWUP-AUROC-RASHOMON-20261021",
        "seeds": SEEDS,
        "seed_pairs": [[a, b] for a, b in SEED_PAIRS],
        "fold_results": fold_results,
        "pooled": pooled,
        "main_results": main,
    }
    (args.followup_dir / "rashomon_results.json").write_text(
        json.dumps(out, indent=2, allow_nan=False) + "\n")
    markdown(fold_results, pooled, main,
             args.followup_dir / "rashomon.md")
    print(json.dumps(pooled, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
