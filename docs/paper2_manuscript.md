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
model-uncertainty signal on weak models). A boundary result qualifies this
headline: the EIG gain is a *function of the pre-filter threshold*, not a
constant. On a strictly pre-filtered pool (MHC presentation threshold 0.9)
the weak surrogate's gain vanishes (1.0x), while on a leniently pre-filtered
pool (threshold 0.5) it survives (1.76x at 7.4% sampling) — showing part of
the raw-pool gain reflected presentation correlation rather than
TCR-specific signal and making the structure-model EIG a necessary
condition. A batch-size ablation further inverts the field's "many small
batches" heuristic: at fixed budget, large batches (3600x2, 24.5%) beat
small ones (300x20, 21.6%). An MHCflurry presentation
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
(Rashomon effect, [Meta-TCR-GNN] §2.6); verdict pending domain adaptation
(mixing predicted-style structures into training — route now AF3 results,
as the TCRmodel2 queue stalled; see §4).

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
- VDJdb-pool emulation: presentation filtering at 0.5 / 0.9 keeps 81.7% /
  67.1% of pairs with positive density unchanged (7.0% -> 6.9% / 7.0%) —
  VDJdb negatives are themselves presented peptides, so the density boost
  expected on a real immunopeptidome pool cannot be emulated on VDJdb.

### 2.5 Pre-filtered pools modulate weak-model EIG gains (boundary result)

The EIG gain is a **function of the pre-filter threshold**, not a constant
(Table 3). On the 0.9-presentation pool, EIG sampling converges to random
(recall 27.0% vs 27.2% at 27% sampling; gain 1.0x across the curve), while
the raw pool showed 2.0x at 6% sampling. Critically, on the *lenient*
0.5-presentation pool the gain *survives* (1.76x at 7.4% sampling, falling
monotonically to 1.28x at 22.2%). Interpretation: part of the raw-pool gain
was presentation-correlation (physicochemical features recognize
low-presentation negatives), which vanishes once the pool is homogenized by
a strict filter — but a lenient filter leaves enough presentation gradient
for the weak surrogate to exploit. Since the real ACDC platform pre-filters
at an intermediate stringency, the gain's survival is an empirical function
of that threshold and cannot be stated as a single headline number.
**Consequence: intelligent sampling with the weak sequence surrogate is
unreliable on the real platform; the structure-model EIG (post domain
adaptation) is a necessary condition, not an optimization.** This boundary
must be reported alongside the 2.0x raw-pool headline.

### 2.6 Proteome-scale candidate list (KN-8 dry-run)

Pipeline complete: 20,431 Swiss-Prot proteins -> 22.5M 9/10-mers ->
MHCflurry -> 4.65M-peptide pool (threshold 0.0122) -> EIG -> top-50k
list with UniProt mapping. Negative result: the weak-model EIG does not
transfer to proteome scale — known A2 binders show zero enrichment in the
top-50k (0/24 vs 0.3 expected, hypergeometric P=1.0) and the Kimmtrak
target peptide ranks in the bottom third. The list is a pipeline
deliverable; the scientifically meaningful ranking requires the
domain-adapted EGNN re-rank.

### 2.7 Batch-size ablation (inverts the field heuristic)

At fixed budget (12% of pool, EIG, 3 seeds, raw pool), large batches beat
small ones — **opposite to the literature's "many small batches" rule**:

| batch size | rounds | final recall |
|---|---|---|
| 3600 | 2 | 24.5% ± 2.3 |
| 1800 | 3 | 19.7% ± 1.4 |
| 900 | 7 | 20.8% ± 0.7 |
| 600 | 10 | 21.9% ± 0.5 |
| 300 | 20 | 21.6% ± 1.0 |

Interpretation: with a weak model refit from a bootstrap ensemble each
round, each batch must be large enough to convey a distinguishable signal;
over-fragmented batches add resampling noise rather than adaptation.
Combined with the low-sampling gain decay ("the earlier the more"), this
shows batch size has a *lower bound* — it is not "the smaller the better".
Large batches carry higher round-to-round variance (±2.3, only 2 rounds),
so both recall and variance must be reported.

### 2.8 (pending) Four-structure native-control test

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
- Domain adaptation (fixing the Rashomon ranking instability, [Meta-TCR-GNN]
  §2.6): mix predicted-style structures into training so the model sees OOD
  inputs during fit. Source structures = AF3 predictions (the TCRmodel2
  queue stalled — 100+ jobs pending on a 2021-era template server — so the
  DA route uses AF3 results directly).
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
- The threshold-dependent gain and the inverted batch-size heuristic are
  both *negative-capability* findings: the field's convenient rules
  ("smaller batches converge faster", "intelligent sampling always helps")
  do not survive contact with weak surrogates and pre-filtered pools.
- Toward real ACDC loops: the pipeline is drop-in for wet-lab data
  (technical route preserved in KEY-NODES.md).

## Figures & Tables

- Fig 1: recall-vs-sampling curves (4 acquisition functions).
- Fig 2: clinical ranking (bar with std whiskers).
- Fig 3: resubmission scatter + std-vs-score.
- Table 1: acquisition comparison.
- Table 2: MAG-IC3 ensemble ranking.
- Table 3: EIG gain vs pre-filter threshold (raw / 0.5 / 0.9 pools).
- Table 4: batch-size ablation (3600/1800/900/600/300 at fixed budget).
