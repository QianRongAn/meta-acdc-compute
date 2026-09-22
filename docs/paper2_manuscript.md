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
TCR-peptide interfaces, and an active-learning layer (EIG acquisition,
epsilon-greedy) selects which pairs to "validate." In a simulated closed
loop on held-out epitopes, EIG sampling achieves 2.0x positive-recall gain
over random sampling at 6% of the pool, outperforming MC-variance (1.6x) and
entropy sampling (no gain — entropy of the mean prediction carries no
model-uncertainty signal on weak models). An MHCflurry presentation
pre-filter — mirroring how the ACDC display library is constructed —
compresses the candidate pool ~50x at a top-2% threshold while retaining
56% of known binders under a composition-shuffled (hard) decoy protocol.
We further quantify the two dominant uncertainty sources: MHC-presentation
prediction and structure-prediction variance (AF3 resubmission can swing a
single-model score by 0.94 AUROC; 5-model ensemble stds reach 0.43), and we
report that cross-reactivity rankings on AF3-predicted structures are not
yet reproducible across model instances at the current data scale — a
stability requirement any OOD-scoring pipeline must satisfy before
deployment. The result is a ranking-based safety-assessment pipeline
(TCR-Safety-Radar) and a reproducible protocol for in silico off-target
screening.

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

### 2.2 Clinical gold-standard ranking (MAG-IC3) — single-instance caveat

Table 2: ensemble ranking of MAGE-A3 target / titin mimic / MAGE-A6 /
MAGE-B18 / random control. Native #1; titin at target level (0.706 vs
0.683); control unexpectedly high (0.792) — caveat discussed. **Downgraded
to single-instance observation**: cross-instance Spearman 0.17/-0.24/NaN
(Rashomon effect, [Meta-TCR-GNN] §2.6); verdict pending domain adaptation.

### 2.3 Uncertainty quantification

- Resubmission variance (Fig): 1/9 jobs swing 0.94; ensemble std 0.01-0.43.
- ipTM vs model score r = 0.333: predictor confidence cannot replace the
  model; borderline cases carry the largest ensemble std (an acquisition
  signal complementary to EIG).

### 2.4 MHC-presentation pre-filter benchmark (ACDC library construction)

- MHCflurry 2.2.1 presentation predictor; 50 IEDB-validated HLA-A*02:01
  binders vs 50 composition-shuffled decoys (hard protocol).
- Top-2% threshold: 56% binder recall at ~50x compression; top-10%: 72% at
  ~10x; mean rank binders 30 vs decoys 71/100.
- VDJdb-pool emulation: presentation filtering at 0.5 keeps 81.7% of pairs
  with positive density unchanged (7.0% -> 6.9%) — VDJdb negatives are
  themselves presented peptides, so the density boost expected on a real
  immunopeptidome pool cannot be emulated on VDJdb; strict-threshold and
  batch-ablation simulations quantify the accessible regimes.

### 2.5 (pending) Batch-size ablation & strict-prefilter pool

- Literature predicts many-small-batches beat few-large at fixed budget;
  simulation grid running (fixed 12% budget; batch 3600/1800/600/300).
- [results pending]

### 2.6 (pending) Four-structure native-control test

With native controls for 1AO7/1QSE/1QSF (A6/B7 family), pooled test of
native-rank significance. [data pending user submissions]

## 3. Methods

- Simulator: pool = 49,696 held-out epitope pairs (7.0% positive); rounds
  of 3,000 with epsilon=0.15; oracle = VDJdb labels; model = bootstrap
  logistic ensemble (20); acquisition = EIG (BALD), predictive entropy,
  MC variance. Pool modes: raw and MHCflurry-presentation-prefiltered
  (KN-7+); batch-size ablation at fixed total budget.
- Structure pipeline: AF3 five-chain submissions (web, terms-compliant),
  mmCIF parsing, interface graphs, EGNN scoring [Meta-TCR-GNN], 5-model
  ensemble mean +/- std.
- MHC pre-filter: MHCflurry 2.2.1 Class1PresentationPredictor; per-allele
  scoring with HLA-A*02:01 fallback; thresholds calibrated on decoy
  quantiles (composition-shuffled binders).
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
