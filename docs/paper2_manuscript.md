# Manuscript draft v0.1 — Active-learning safety assessment

> Working draft, 2026-09-22. Dry-lab scope. Numbers from docs/benchmarks.md.
> Companion to paper1 (Meta-TCR-GNN); cites its model as [Meta-TCR-GNN].

## Title (candidate)

**How few experiments suffice? An active-learning digital twin for
proteome-wide safety assessment of therapeutic TCRs**

## Abstract

Proteome-wide screening of therapeutic TCR cross-reactivity is physically
intractable (>10^7 peptides x >10^15 TCRs). We build a digital twin of the
dry-wet loop: a structural compatibility model [Meta-TCR-GNN] scores
AlphaFold 3-predicted TCR-peptide interfaces, and an active-learning layer
(EIG acquisition, epsilon-greedy) selects which pairs to "validate." In a
simulated closed loop on held-out epitopes, EIG sampling achieves 2.0x
positive-recall gain over random sampling at 6% of the pool, outperforming
MC-variance (1.6x) and entropy sampling (no gain — entropy of the mean
prediction carries no model-uncertainty signal on weak models). Applied to
the clinically lethal MAG-IC3 TCR, the model ranks the native MAGE-A3 target
first and scores the titin mimic — responsible for two deaths — at the
target's own level under 5-model ensembling. We quantify the dominant
uncertainty source: AF3 prediction variance (resubmission can swing a
single-model score by 0.94 AUROC; 5-model ensemble stds reach 0.43). The
result is a ranking-based safety-assessment pipeline (TCR-Safety-Radar) and
a reproducible protocol for in silico off-target screening.

## 1. Introduction (skeleton)

- Safety crisis & screening economics [paper1 intro + proposal Section III].
- Active learning in protein engineering: ALDE, ALSEBO, EVOLVEpro [CIT].
- Our framing: ranking, not absolute scoring (from paper1's 5BRZ/5BS0
  negative result).

## 2. Results

### 2.1 Acquisition-function comparison (simulated closed loop)

Table 1: recall at 6/12/18/24% sampling for random, entropy, MC-variance,
EIG (2 seeds). EIG 10.4% vs random 5.9% at 6% (2.0x); variance 9.3%;
entropy 5.8% (~random). Diminishing gains with sampling fraction
(expected depletion of easy positives).

### 2.2 Clinical gold-standard ranking (MAG-IC3)

Table 2: ensemble ranking of MAGE-A3 target / titin mimic / MAGE-A6 /
MAGE-B18 / random control. Native #1; titin at target level (0.706 vs
0.683); control unexpectedly high (0.792) — caveat discussed.

### 2.3 Uncertainty quantification

- Resubmission variance (Fig): 1/9 jobs swing 0.94; ensemble std 0.01-0.43.
- ipTM vs model score r = 0.333: predictor confidence cannot replace the
  model; borderline cases carry the largest ensemble std (an acquisition
  signal complementary to EIG).

### 2.4 (pending) Four-structure native-control test

With native controls for 1AO7/1QSE/1QSF (A6/B7 family), pooled test of
native-rank significance. [data pending user submissions]

## 3. Methods

- Simulator: pool = 49,696 held-out epitope pairs (7.0% positive); rounds
  of 3,000 with epsilon=0.15; oracle = VDJdb labels; model = bootstrap
  logistic ensemble (20); acquisition = EIG (BALD), predictive entropy,
  MC variance.
- Structure pipeline: AF3 five-chain submissions (web, terms-compliant),
  mmCIF parsing, interface graphs, EGNN scoring [Meta-TCR-GNN], 5-model
  ensemble mean +/- std.
- Dashboard: TCR-Safety-Radar (zero-dependency web server + heatmap/ranking
  UI).

## 4. Discussion (skeleton)

- Digital-twin validity: simulator uses weak sequence oracle — the gains
  bound what the structure model can do with better scores.
- Entropy-sampling failure as a methodological warning.
- AF3 variance as the field's reproducibility issue; ensemble reporting
  as the minimum standard.
- Toward real ACDC loops: the pipeline is drop-in for wet-lab data
  (technical route preserved in KEY-NODES.md).

## Figures & Tables

- Fig 1: recall-vs-sampling curves (4 acquisition functions).
- Fig 2: clinical ranking (bar with std whiskers).
- Fig 3: resubmission scatter + std-vs-score.
- Table 1: acquisition comparison.
- Table 2: MAG-IC3 ensemble ranking.
