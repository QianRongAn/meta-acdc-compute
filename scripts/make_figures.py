"""Generate paper-1 figures (matplotlib, publication style).

Palette: reference dataviz palette (blue categorical slot, single-hue
sequential for the heatmap). Output: docs/figures/*.png (300 dpi) + .svg.

Usage:
    .venv/bin/python scripts/make_figures.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#2a78d6"
GRAY = "#b9b8b3"
INK = "#0b0b0b"
MUTED = "#52514e"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.edgecolor": "#d8d7d2",
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
})


def fig1_ablation() -> None:
    """Feature-ablation bar chart (graft task)."""
    configs = [
        ("Full (v9.1)", 0.933, True),
        ("− side-chain extent", 0.924, False),
        ("− residue identity", 0.899, False),
        ("− physicochemical", 0.606, False),
        ("− contact histogram", 0.501, False),
        ("Cα-only (v5-v8)", 0.500, False),
        ("Sequence baseline", 0.495, False),
    ]
    labels = [c for c, _, _ in configs]
    vals = [v for _, v, _ in configs]
    colors = [BLUE if full else GRAY for _, _, full in configs]

    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    bars = ax.barh(labels[::-1], vals[::-1], color=colors[::-1], height=0.62)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Graft-decoy AUROC (held-out complexes)")
    ax.axvline(0.5, color="#d8d7d2", linewidth=1, linestyle="--")
    ax.text(0.5, -0.42, "chance", color=MUTED, fontsize=7, ha="center")
    for bar, v in zip(bars, vals[::-1]):
        ax.text(bar.get_width() + 0.015, bar.get_y() + bar.get_height() / 2,
                f"{v:.3f}", va="center", fontsize=8, color=INK)
    ax.tick_params(axis="y", length=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_ablation.png", dpi=300)
    fig.savefig(OUT / "fig1_ablation.svg")
    plt.close(fig)
    print("fig1_ablation done")


def fig2_heatmap() -> None:
    """TCR x peptide ensemble-score heatmap (single-hue sequential)."""
    scores = {}
    with open(ROOT / "data/processed/prediction_scores_ensemble.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pdb, pep = r["job_id"].split("_", 1)
            scores.setdefault(pdb, {})[pep] = float(r["score"])

    tcrs = ["1ao7", "1qrn", "1g6r", "1mwa", "1oga", "5brz"]
    peps = sorted({p for d in scores.values() for p in d})
    import numpy as np
    mat = np.full((len(tcrs), len(peps)), np.nan)
    for i, t in enumerate(tcrs):
        for j, p in enumerate(peps):
            mat[i, j] = scores.get(t, {}).get(p, np.nan)

    fig, ax = plt.subplots(figsize=(8.5, 3.0))
    cmap = plt.cm.Blues.copy()
    cmap.set_bad("#ececea")
    im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(peps)), peps, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(tcrs)), tcrs)
    for i in range(len(tcrs)):
        for j in range(len(peps)):
            v = mat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=6.5,
                        color="white" if v > 0.55 else INK)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("compatibility score (5-model ensemble)")
    ax.set_xlabel("peptide")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_heatmap.png", dpi=300)
    fig.savefig(OUT / "fig2_heatmap.svg")
    plt.close(fig)
    print("fig2_heatmap done")


def fig3_resubmission() -> None:
    """First vs resubmitted AF3 scores (identity scatter)."""
    pairs = [
        ("LGYGFVNYI", 0.9864, 0.992), ("LLFGFPVYV", 0.9927, 0.999),
        ("LLFGKPVYV", 0.7821, 0.952), ("LLFGPVYV", 0.0018, 0.086),
        ("LLFGYAVYV", 0.9841, 0.991), ("LLFGYPRYV", 0.0486, 0.991),
        ("LLFGYPVAV", 0.9595, 0.887), ("MLWGYLQYV", 0.9884, 0.986),
        ("LLFGYPVYV(1qrn)", 0.9991, 0.996),
    ]
    fig, ax = plt.subplots(figsize=(3.8, 3.8))
    ax.plot([0, 1], [0, 1], color="#d8d7d2", linewidth=1)
    for name, v1, v2 in pairs:
        ax.scatter(v1, v2, s=22, color=BLUE, zorder=3, edgecolor="white",
                   linewidth=0.5)
        if abs(v2 - v1) > 0.15:
            ax.annotate(name.split("(")[0], (v1, v2),
                        textcoords="offset points", xytext=(4, 4),
                        fontsize=6.5, color=INK)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("first AF3 submission (model 0)")
    ax.set_ylabel("resubmission (model 0)")
    ax.set_title("Resubmission variance", fontsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_resubmission.png", dpi=300)
    fig.savefig(OUT / "fig3_resubmission.svg")
    plt.close(fig)
    print("fig3_resubmission done")


def fig5_instability() -> None:
    """Cross-instance instability (Rashomon): per-job score profiles for 4
    independently trained instances + pairwise Spearman matrix."""
    import numpy as np

    files = [
        ("v9.1 (trained w/ pLDDT)", "inst_scores_egnn_dataset_ensemble.tsv"),
        ("no-pLDDT seed A", "inst_scores_egnn_noplddt_s10_ensemble.tsv"),
        ("no-pLDDT seed B", "inst_scores_egnn_noplddt_s20_ensemble.tsv"),
        ("no-pLDDT seed C", "inst_scores_egnn_noplddt_s30_ensemble.tsv"),
    ]
    per = []
    for _, f in files:
        d = {}
        with open(ROOT / "data/processed" / f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["score"]:
                    d[r["job_id"]] = float(r["score"])
        per.append(d)
    jobs = sorted(per[0], key=lambda j: -per[0][j])

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.4),
                             gridspec_kw={"width_ratios": [2.1, 1]})
    ax = axes[0]
    x = np.arange(len(jobs))
    for k, (label, _) in enumerate(files):
        vals = np.array([per[k].get(j, np.nan) for j in jobs])
        ax.plot(x, vals, marker="o", markersize=3, linewidth=1.2,
                color=BLUE if k == 0 else GRAY, alpha=1.0 if k == 0 else 0.85,
                label=label, zorder=3 if k == 0 else 2)
    ax.set_xlabel("AF3 jobs (sorted by v9.1 score)")
    ax.set_ylabel("compatibility score")
    ax.set_ylim(-0.04, 1.04)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.set_title("A  Four instances, same 21 AF3 inputs", fontsize=9,
                 loc="left")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax = axes[1]
    from scipy.stats import spearmanr
    n = len(files)
    M = np.array([[per[k].get(j, np.nan) for j in jobs] for k in range(n)])
    R = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i != j:
                r, _ = spearmanr(M[i], M[j])
                R[i, j] = r
            else:
                R[i, j] = 1.0
    cmap = plt.cm.Blues.copy()
    cmap.set_bad("#ececea")
    im = ax.imshow(R, cmap=cmap, vmin=-0.5, vmax=1.0)
    ax.set_xticks(range(n), [f"S{k+1}" for k in range(n)], fontsize=8)
    ax.set_yticks(range(n), [f"S{k+1}" for k in range(n)], fontsize=8)
    for i in range(n):
        for j in range(n):
            v = R[i, j]
            txt = "n/d" if np.isnan(v) else f"{v:.2f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7,
                    color="white" if v > 0.55 else INK)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("Spearman r", fontsize=8)
    ax.set_title("B  Ranking agreement", fontsize=9, loc="left")
    fig.suptitle("Cross-instance instability (Rashomon effect)", fontsize=10,
                 x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "fig5_instability.png", dpi=300)
    fig.savefig(OUT / "fig5_instability.svg")
    plt.close(fig)
    print("fig5_instability done")


def figp2_recall() -> None:
    """Paper-2 Fig 1: recall vs sampling fraction, 4 acquisition functions."""
    data = {
        "random": [(6, 5.9), (24, 24.3)],
        "entropy": [(6, 5.8), (24, 24.2)],
        "MC variance": [(6, 9.3), (24, 28.8)],
        "EIG (BALD)": [(6, 10.4), (24, 33.0)],
    }
    colors = {"random": GRAY, "entropy": "#8ab6e8", "MC variance": "#5a91dd",
              "EIG (BALD)": BLUE}
    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    for name, pts in data.items():
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.plot(xs, ys, marker="o", markersize=4.5, linewidth=1.4,
                color=colors[name], label=name)
        ax.annotate(f"{ys[-1]:.1f}%", (xs[-1], ys[-1]),
                    textcoords="offset points", xytext=(5, 0), fontsize=7,
                    color=colors[name])
    ax.set_xlabel("sampling fraction of pool (%)")
    ax.set_ylabel("positive recall (%)")
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 38)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_title("Simulated dry-wet loop (held-out epitopes)", fontsize=9,
                 loc="left")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "figp2_recall.png", dpi=300)
    fig.savefig(OUT / "figp2_recall.svg")
    plt.close(fig)
    print("figp2_recall done")


def figp2_prefilter() -> None:
    """Paper-2: MHCflurry pre-filter recall/compression vs threshold."""
    thresholds = [0.5, 2.0, 5.0, 10.0]
    recall = [0.38, 0.56, 0.60, 0.72]
    compression = [200, 50, 20, 10]
    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    ax.plot(thresholds, recall, marker="o", markersize=4.5, linewidth=1.4,
            color=BLUE)
    for t, r, c in zip(thresholds, recall, compression):
        ax.annotate(f"{r:.0%} @ {c}x", (t, r), textcoords="offset points",
                    xytext=(4, 7), fontsize=7, color=MUTED)
    ax.set_xlabel("keep threshold (top x% of decoy-score distribution)")
    ax.set_ylabel("known-binder recall")
    ax.set_ylim(0, 0.85)
    ax.set_xticks(thresholds)
    ax.set_title("MHC presentation pre-filter (hard decoys)", fontsize=9,
                 loc="left")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "figp2_prefilter.png", dpi=300)
    fig.savefig(OUT / "figp2_prefilter.svg")
    plt.close(fig)
    print("figp2_prefilter done")


def main() -> int:
    fig1_ablation()
    fig2_heatmap()
    fig3_resubmission()
    fig5_instability()
    figp2_recall()
    figp2_prefilter()
    return 0


if __name__ == "__main__":
    sys.exit(main())
