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


def main() -> int:
    fig1_ablation()
    fig2_heatmap()
    fig3_resubmission()
    return 0


if __name__ == "__main__":
    sys.exit(main())
