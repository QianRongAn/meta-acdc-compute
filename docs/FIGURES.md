# Figures — manifest & submission checklist

> Repo-generated figures for paper1 (Meta-TCR-GNN) and paper2 (AL safety loop).
> Regenerate everything with `python scripts/make_figures.py` (venv, CPU).
> All figures are 300-dpi PNG + vector SVG (text as paths, submission-safe).
> Last full rebuild: **2026-10-15** (post-attribution-repair data).

## Manifest

| # | File(s) | Used in | Panel content | Data source (regenerable) | Fit |
|---|---|---|---|---|---|
| Fig 1 | `fig1_ablation.{png,svg}` | paper1 §2.2 | Cα inseparability: graft AUROC vs data size + hypothesis-elimination series | `data/processed` KN-4 ablation logs | single |
| Fig 2 | `fig2_heatmap.{png,svg}` | paper1 §2.1/§2.3 | per-TCR × peptide compatibility heatmap (26 families) | `prediction_scores_ensemble.tsv` (188-panel ensemble, post-repair) | double |
| Fig 3 | `fig3_resubmission.{png,svg}` | paper1 §2.5 | AF3 resubmission variance: ipTM stable vs fixed-instance EGNN range | `af3_resubmission_variance.tsv` | single |
| Fig 4 | `fig5_instability.{png,svg}` (legacy number) | paper1 §2.6 | Rashomon crisis: 4 independent instances' score profiles + pairwise Spearman matrix | `inst_scores_*_ensemble.tsv` (Sept snapshot) | double |
| Fig 5 | `fig6_domain_adapt.{png,svg}` (legacy number) | paper1 §2.6 | DA resolution: 5 seeds × 188 candidates, no collapse; Spearman matrix, mean r=0.644 | `da_scores_seed{0..4}_ensemble.tsv` | double |
| Fig 6 | `fig7_crossreactivity.{png,svg}` (legacy number) | paper1 §2.4/intro | VDJdb cross-reactivity landscape (cross-epitope breadth per TCR) | `structure_vdjdb_map.tsv` | single |
| Fig 7 | `fig8_cross_method.{png,svg}` | paper1 §2.6 (negative control) | AF3 vs ColabFold (AF2-multimer) DA scores, same 28 jobs: Pearson 0.184 | `prediction_scores_ensemble.tsv` × `colabfold_scores_ensemble.tsv` | single |
| Fig P2-1 | `figp2_recall.{png,svg}` | paper2 §2.1 | recall vs sampling fraction, 4 acquisition functions (AL simulator) | benchmarks AL-ablation numbers | single |
| Fig P2-2 | `figp2_prefilter.{png,svg}` | paper2 §2.7 | HLA prefilter compression / per-TCR density (dataset-internal) | benchmarks prefilter numbers | single |

## Notes for camera-ready

- **Legacy numbering**: `fig5`/`fig6` filenames predate the current paper1 figure list
  (paper1 lists them as Fig 5/Fig 6 with the content above); the numbering in the
  manuscript captions is authoritative, filenames are stable identifiers.
- **Column fit / print widths (Nature-style: 180 mm full width = 7.09 in;
  88 mm single column = 3.46 in)**: import every figure at exactly one of the
  widths below. The SVG is vector and scales losslessly (Word 2016+ imports SVG
  directly; the 300-dpi PNG alternative is already print-ready at 7-9 in
  native size, i.e. 2100+ px, so scaling to 88-180 mm is a document-side
  setting, not a re-render).

| Figure | import width | LaTeX |
|---|---|---|
| fig1_ablation | 88 mm | `\includegraphics[width=0.48\textwidth]` |
| fig2_heatmap | 180 mm | `\includegraphics[width=\textwidth]` |
| fig3_resubmission | 88 mm | `\includegraphics[width=0.48\textwidth]` |
| fig5_instability | 180 mm | `\includegraphics[width=\textwidth]` |
| fig6_domain_adapt | 180 mm | `\includegraphics[width=\textwidth]` |
| fig7_crossreactivity | 88 mm | `\includegraphics[width=0.48\textwidth]` |
| fig8_cross_method | 88 mm | `\includegraphics[width=0.48\textwidth]` |
| figp2_recall | 88 mm | `\includegraphics[width=0.48\textwidth]` |
| figp2_prefilter | 88 mm | `\includegraphics[width=0.48\textwidth]` |
- **Do not re-crop**: the 2026-10-14 repair changed panel contents (154→188) and
  all scores; figures regenerated from the repaired data supersede the Sept-23
  versions (kept only `figp2_*` unchanged: paper-2 AL numbers are unaffected).
- Reproduce: `.venv/bin/python scripts/make_figures.py` (~1 min, CPU-only).
