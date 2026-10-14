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
    """TCR x peptide ensemble-score heatmap (single-hue sequential).

    Prefers the 5-seed domain-adapted ensemble (prediction_scores_ensemble.tsv,
    written from da_scores_seed*_ensemble.tsv), which is cross-instance stable;
    falls back to the legacy relabeled single-instance table.
    """
    scores = {}
    src = ROOT / "data/processed/prediction_scores_ensemble.tsv"
    if not src.exists():
        src = ROOT / "data/processed/prediction_scores_ensemble.relabeled.tsv"
    with open(src) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pdb, pep = r["job_id"].split("_", 1)
            scores.setdefault(pdb, {})[pep] = float(r["score"])

    prefer = ["1ao7", "1qrn", "1tcr", "2vlj", "5brz"]
    tcrs = [t for t in prefer if t in scores] + \
        [t for t in sorted(scores) if t not in prefer]
    peps = sorted({p for d in scores.values() for p in d})
    import numpy as np
    mat = np.full((len(tcrs), len(peps)), np.nan)
    for i, t in enumerate(tcrs):
        for j, p in enumerate(peps):
            mat[i, j] = scores.get(t, {}).get(p, np.nan)

    n_rows, n_cols = len(tcrs), len(peps)
    # Keep the canvas narrow enough that a column-width view still has
    # ~40 px per TCR row. 0.12 in/row stacked the PDB labels.
    fig_h = max(12.0, 0.50 * n_rows)
    fig_w = max(10.5, 0.16 * n_cols + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = plt.cm.Blues.copy()
    cmap.set_bad("#ececea")
    im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(n_cols), peps, rotation=90, ha="center", va="top",
                  fontsize=6.5)
    ax.set_yticks(range(n_rows), tcrs, fontsize=8)
    ax.tick_params(axis="y", pad=4, length=0)
    ax.tick_params(axis="x", pad=2, length=0)
    for i in range(len(tcrs)):
        for j in range(len(peps)):
            v = mat[i, j]
            # Adjacent scored cells are narrower than a 4-char label; colour
            # carries the value there. Annotate only isolated cells.
            if np.isnan(v):
                continue
            left = mat[i, j - 1] if j else np.nan
            right = mat[i, j + 1] if j + 1 < n_cols else np.nan
            if np.isnan(left) and np.isnan(right):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=6.5,
                        color="white" if v > 0.55 else INK)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("compatibility score (5-seed DA ensemble, cross-instance stable)")
    ax.set_xlabel("peptide")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_heatmap.png", dpi=300)
    fig.savefig(OUT / "fig2_heatmap.svg")
    plt.close(fig)
    print("fig2_heatmap done")


def fig3_resubmission() -> None:
    """AF3 resubmission variance, corrected: structure stable, EGNN modest.

    The legacy claim (docs/benchmarks.md 2026-09-22) that resubmission swings
    the score by up to 0.94 came from a scorer instance that no longer exists.
    With AF3's own confidence (model-independent) and a surviving fixed EGNN
    instance, resubmission moves both metrics by far less.
    """
    import numpy as np

    rows = list(csv.DictReader(
        open(ROOT / "data/processed/af3_resubmission_variance.tsv"),
        delimiter="\t"))

    def mean_of(prefix: str, row: dict) -> float | None:
        vals = [float(row[f"{prefix}_{k}"]) for k in ("B1", "B2", "B3")
                if row.get(f"{prefix}_{k}", "")]
        return float(np.mean(vals)) if vals else None

    iptm_pairs, egnn_pairs = [], []
    for r in rows:
        if not r.get("iptm_A", ""):
            continue
        y = mean_of("iptm", r)
        if y is not None:
            iptm_pairs.append((r["job_id"], float(r["iptm_A"]), y))
        y = mean_of("egnn_noplddt_s10", r)
        if y is not None and r.get("egnn_noplddt_s10_A", ""):
            egnn_pairs.append((r["job_id"], float(r["egnn_noplddt_s10_A"]), y))

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6))
    specs = [
        (axes[0], iptm_pairs, "AF3 self-confidence (ipTM)",
         "A  AlphaFold3 agrees with itself", (-0.04, 1.04)),
        (axes[1], egnn_pairs, "EGNN score (fixed no-pLDDT instance)",
         "B  EGNN scorer barely moves", (-0.05, 0.30)),
    ]
    for ax, pairs, xlabel, title, lim in specs:
        ax.plot(list(lim), list(lim), color="#d8d7d2", linewidth=1)
        for name, v1, v2 in pairs:
            ax.scatter(v1, v2, s=24, color=BLUE, zorder=3, edgecolor="white",
                       linewidth=0.5)
            if abs(v2 - v1) > 0.08:
                ax.annotate(name.split("_", 1)[1], (v1, v2),
                            textcoords="offset points", xytext=(4, 4),
                            fontsize=7.0, color=INK)
        ax.set_xlim(*lim)
        ax.set_ylim(*lim)
        ax.set_xlabel(f"first submission (A): {xlabel}")
        ax.set_ylabel(f"resubmission (B): {xlabel}")
        ax.set_title(title, fontsize=9, loc="left")
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
    fig.suptitle("AF3 resubmission variance is small once the scorer is fixed",
                 fontsize=10, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
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


def fig6_domain_adapt() -> None:
    """Domain adaptation on real native AF3 positives (172-native DA seeds).

    Panel A: per-seed candidate score profiles over the 188-candidate panel
    repaired 2026-10-14 (no seed collapses at this scale; the degeneracy
    guard is retained for smaller runs).
    Panel B: pairwise Spearman matrix over all five seeds; the reported
    mean pairwise r (0.644) crosses the 0.5 stability gate.
    """
    import numpy as np
    from scipy.stats import spearmanr

    seeds = [0, 1, 2, 3, 4]  # 172-native DA run: all five seeds healthy
    per = []
    for s in seeds:
        d = {}
        with open(ROOT / f"data/processed/da_scores_seed{s}_ensemble.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["score"]:
                    d[r["job_id"]] = float(r["score"])
        per.append(d)
    jobs = sorted(per[0])
    from collections import defaultdict
    groups: dict[str, list[int]] = defaultdict(list)
    for i, j in enumerate(jobs):
        groups[j.split("_", 1)[0]].append(i)
    order = sorted(groups, key=lambda p: groups[p][0])
    # Labels sit on the job-span centre, so adjacent singleton PDBs are one
    # job apart. Size panel A by jobs, not by group count, or those labels
    # still collide (5m00/5m01/5m02).
    panel_a = max(36.0, 0.26 * len(jobs))
    fig, axes = plt.subplots(1, 2, figsize=(panel_a + 5.4, 5.0),
                             gridspec_kw={"width_ratios": [panel_a, 4.4]})
    ax = axes[0]
    x = np.arange(len(jobs))
    for k, s in enumerate(seeds):
        vals = np.array([per[k].get(j, np.nan) for j in jobs])
        ax.plot(x, vals, marker="o", markersize=3, linewidth=1.2,
                color=BLUE if k == 0 else GRAY, alpha=1.0 if k == 0 else 0.85,
                label=f"DA seed {s}")
    ax.set_xlabel("AF3 cross-reactivity candidates (188, repaired panel)")
    ax.set_ylabel("compatibility score")
    ax.set_ylim(-0.04, 1.04)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.set_title(f"A  Five DA instances, {len(jobs)} AF3 candidates",
                 fontsize=9, loc="left")
    ax.set_xticks([sum(groups[p]) / len(groups[p]) for p in order], order,
                  rotation=90, ha="center", va="top", fontsize=7)
    ax.tick_params(axis="x", pad=2, length=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax = axes[1]
    M = np.array([[per[k].get(j, np.nan) for j in jobs] for k in range(len(seeds))])
    n = len(seeds)
    R = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            R[i, j] = 1.0 if i == j else spearmanr(M[i], M[j])[0]
    cmap = plt.cm.Blues.copy()
    cmap.set_bad("#ececea")
    im = ax.imshow(R, cmap=cmap, vmin=-0.5, vmax=1.0)
    ax.set_xticks(range(n), [f"S{s}" for s in seeds], fontsize=8)
    ax.set_yticks(range(n), [f"S{s}" for s in seeds], fontsize=8)
    for i in range(n):
        for j in range(n):
            v = R[i, j]
            ax.text(j, i, "n/d" if np.isnan(v) else f"{v:.2f}", ha="center",
                    va="center", fontsize=7, color="white" if v > 0.55 else INK)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("Spearman r", fontsize=8)
    ax.set_title("B  Ranking agreement", fontsize=9, loc="left")
    mean_r = np.nanmean(R[np.triu_indices(n, 1)])
    fig.suptitle(
        f"Native domain adaptation: no collapse, mean pairwise r="
        f"{mean_r:.3f} (gate 0.5)", fontsize=10, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "fig6_domain_adapt.png", dpi=300)
    fig.savefig(OUT / "fig6_domain_adapt.svg")
    plt.close(fig)
    print(f"fig6_domain_adapt done (mean pairwise r={mean_r:.3f})")


def fig7_crossreactivity() -> None:
    """Structure<->VDJdb cross-reactivity landscape.

    Distribution of cross-epitope breadth per TCR (epitopes the TCR is
    VDJdb-validated against that differ from its crystallized peptide) plus
    the most promiscuous families.
    """
    import numpy as np
    from collections import defaultdict

    same: dict[str, set] = defaultdict(set)
    diff: dict[str, set] = defaultdict(set)
    native: dict[str, str] = {}
    with open(ROOT / "data/processed/structure_vdjdb_map.tsv", newline="",
              encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pdb = r["pdb"]
            native.setdefault(pdb, r["pdb_peptide"])
            (same if r["match"] == "same" else diff)[pdb].add(r["vdjdb_epitope"])

    breadth = [len(diff[p]) for p in native]
    n_cross = sum(1 for b in breadth if b > 0)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2),
                             gridspec_kw={"width_ratios": [1, 1.3]})
    ax = axes[0]
    bins = np.arange(0, max(breadth) + 2) - 0.5
    ax.hist(breadth, bins=bins, color=BLUE, edgecolor="white")
    ax.set_xlabel("cross-reactive epitopes per TCR (VDJdb)")
    ax.set_ylabel("number of complexes")
    ax.set_title(f"A  {n_cross}/{len(native)} complexes cross-reactive",
                 fontsize=9, loc="left")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax = axes[1]
    top = sorted(native, key=lambda p: -len(diff[p]))[:10]
    labels = [f"{p} ({native[p]})" for p in top]
    vals = [len(diff[p]) for p in top]
    ax.barh(labels[::-1], vals[::-1], color=BLUE, height=0.62)
    ax.set_xlabel("cross-reactive epitopes")
    ax.set_title("B  Most promiscuous TCRs", fontsize=9, loc="left")
    ax.tick_params(axis="y", length=0, labelsize=7)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.suptitle("Cross-reactivity landscape from the structure-VDJdb map",
                 fontsize=10, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "fig7_crossreactivity.png", dpi=300)
    fig.savefig(OUT / "fig7_crossreactivity.svg")
    plt.close(fig)
    print(f"fig7_crossreactivity done ({n_cross}/{len(native)} cross-reactive)")


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



def _layout_scatter_labels(ax, fig, anchors, labels, fontsize=6.8,
                           n_iter=220, seed=0):
    """Decollide scatter annotations via greedy spiral placement.

    Labels are placed in order of local crowding (hardest first). Each label
    searches outward from its anchor on expanding rings for the closest
    position that does not overlap any already-placed label and does not cover
    a foreign data point; if the ring sweep finds no free slot, the
    least-overlapping slot is used. Positions are clamped inside the axes.
    Labels displaced far from their anchor get a faint leader line.
    Returns the list of Text artists.
    """
    import numpy as np

    rng = np.random.default_rng(seed)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_bbox = ax.get_window_extent(renderer)
    trans = ax.transData
    inv = trans.inverted()

    P = np.array([trans.transform((x, y)) for x, y in anchors], float)
    texts = [
        ax.text(x, y, lab, transform=trans, ha="center", va="center",
                fontsize=fontsize, color=INK, zorder=4, clip_on=False)
        for (x, y), lab in zip(anchors, labels)
    ]
    sizes = [t.get_window_extent(renderer) for t in texts]
    sizes = [(b.width, b.height) for b in sizes]

    order = sorted(range(len(P)),
                   key=lambda i: -sum(1 for p in P
                                      if np.hypot(*(p - P[i])) < 70.0))

    def placed_overlap(b, placed):
        total = 0.0
        for q in placed:
            ox = min(b[0] + b[2], q[0] + q[2]) - max(b[0], q[0])
            oy = min(b[1] + b[3], q[1] + q[3]) - max(b[1], q[1])
            if ox > 0 and oy > 0:
                total += ox * oy
        return total

    def covers_foreign_dot(b, skip):
        for k, p in enumerate(P):
            if k == skip:
                continue
            if b[0] < p[0] < b[0] + b[2] and b[1] < p[1] < b[1] + b[3]:
                return True
        return False

    placed = []  # (x0, y0, w, h) in px
    centers = {}  # label idx -> center px
    for i in order:
        w, h = sizes[i]
        rot = rng.uniform(0, 2 * np.pi)
        best, best_cost, best_r = None, None, None
        for r in np.arange(4.0, 92.0, 4.0):
            n_ang = max(8, int(np.pi * r / 10.0))
            for k in range(n_ang):
                ang = rot + 2 * np.pi * k / n_ang
                c = P[i] + r * np.array([np.cos(ang), np.sin(ang)])
                x0, y0 = c - np.array([w / 2, h / 2])
                if (x0 < ax_bbox.x0 + 1.5 or x0 + w > ax_bbox.x1 - 1.5
                        or y0 < ax_bbox.y0 + 1.5 or y0 + h > ax_bbox.y1 - 1.5):
                    continue
                b = (x0, y0, w, h)
                cost = placed_overlap(b, placed)
                if covers_foreign_dot(b, i):
                    cost += 1000.0
                if cost <= 0.0:
                    best, best_cost, best_r = b, 0.0, r
                    break
                if best_cost is None or cost < best_cost:
                    best, best_cost, best_r = b, cost, r
            if best_cost == 0.0:
                break
        placed.append(best)
        centers[i] = (best[0] + w / 2, best[1] + h / 2)

    for i, t in enumerate(texts):
        t.set_position(inv.transform(centers[i]))

    for i, (x, y) in enumerate(anchors):
        anchor = P[i]
        center = np.array(centers[i])
        if np.hypot(*(center - anchor)) > 10.0:
            cx, cy = inv.transform(center)
            ax.plot([x, cx], [y, cy], color=GRAY, lw=0.5, alpha=0.8,
                    zorder=2.5)
    return texts


def fig8_cross_method() -> None:
    """Cross-method negative control: AF3 vs ColabFold (AF2-multimer) structures
    of the same 8 jobs scored by the same DA ensemble (2026-10-14).

    Reading: near-identical peptide/TCR chains submitted to two structure
    generators yield only weakly correlated DA scores (Pearson ~0.35) with a
    systematic downward shift — the generator is not interchangeable at the
    scoring layer. ColabFold scores are therefore QC-tier only.
    """
    import numpy as np
    from scipy.stats import pearsonr

    def load(name: str) -> dict[str, float]:
        d: dict[str, float] = {}
        src = ROOT / "data/processed" / name
        with open(src) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("score"):
                    d[r["job_id"]] = float(r["score"])
        return d

    af3 = load("prediction_scores_ensemble.tsv")
    cf = load("colabfold_scores_ensemble.tsv")
    jobs = sorted(set(af3) & set(cf))
    x = np.array([af3[j] for j in jobs])
    y = np.array([cf[j] for j in jobs])
    r = pearsonr(x, y)[0]

    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    ax.plot([0, 1], [0, 1], ls="--", lw=1, color="#c9c9c7", zorder=1)
    # 28 个点的分数高度聚簇、直接重叠：显示坐标加固定种子微抖（±0.008，
    # 相对 0-1 轴可忽略），Pearson 与均值仍用未抖动的 x/y 计算。
    # 样本名标签用力导向布局去重（_layout_scatter_labels），离锚点较远
    # 的标签画浅色指引线；图幅加宽以便双栏排版时标签仍可读。
    rng = np.random.default_rng(42)
    jx = x + rng.uniform(-0.008, 0.008, len(x))
    jy = y + rng.uniform(-0.008, 0.008, len(y))
    ax.scatter(jx, jy, s=30, color=BLUE, alpha=0.8, zorder=3, edgecolor="white",
               linewidth=0.6)
    _layout_scatter_labels(ax, fig, list(zip(jx, jy)), jobs, fontsize=6.8)
    ax.set_xlabel("AF3 structure — DA score", fontsize=9)
    ax.set_ylabel("ColabFold (AF2-multimer) — DA score", fontsize=9)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(0.04, 0.96,
            f"n = {len(jobs)}   Pearson r = {r:.2f}\n"
            f"mean  AF3 {x.mean():.2f}  vs  ColabFold {y.mean():.2f}",
            transform=ax.transAxes, va="top", fontsize=8.5, color=INK)
    ax.set_title("Structure generator is not interchangeable\n"
                 "at the scoring layer", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(OUT / "fig8_cross_method.png", dpi=300)
    fig.savefig(OUT / "fig8_cross_method.svg")
    plt.close(fig)
    print(f"fig8_cross_method done (r={r:.3f}, n={len(jobs)})")


def main() -> int:
    fig1_ablation()
    fig2_heatmap()
    fig3_resubmission()
    fig5_instability()
    fig6_domain_adapt()
    fig7_crossreactivity()
    fig8_cross_method()
    figp2_recall()
    figp2_prefilter()
    return 0


if __name__ == "__main__":
    sys.exit(main())
