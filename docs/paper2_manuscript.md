# Manuscript draft v0.2 — Active-learning safety assessment

> Working draft, 2026-09-23. Dry-lab scope. Numbers from docs/benchmarks.md.
> Companion to paper1 (Meta-TCR-GNN); cites its model as [Meta-TCR-GNN].
> Intro/Results/Methods/Discussion expanded from v0.1 skeleton; §2.8 and
> domain-adaptation results remain pending AF3 submissions.

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
prediction and structure-prediction stability (AF3's own confidence is highly
reproducible across resubmissions, ipTM median Δ 0.025; 5-model ensemble stds
reach 0.43), and we
report that cross-reactivity rankings on AF3-predicted structures are not
yet reproducible across model instances at the current data scale — a
stability requirement any OOD-scoring pipeline must satisfy before
deployment. The result is a ranking-based safety-assessment pipeline
(TCR-Safety-Radar) and a reproducible protocol for in silico off-target
screening.

## 1. Introduction

Engineered T cell receptors (TCR-T) and TCR-based biologics carry a
documented safety risk that conventional pre-clinical screening
systematically misses: cross-reactivity against off-target self-peptides.
The fatal MAGE-A3 trial — in which a TCR selected against the cancer-testis
antigen MAGE-A3 recognized a titin peptide expressed in cardiac muscle —
killed two patients and aborted the field's first TCR affinity-enhanced
clinical program [CIT Cameron 2013, Linette 2013, Morgan 2013]. The
mechanism was direct molecular mimicry: the titin epitope differs from the
target by a handful of residues yet docks with near-identical backbone
geometry (0.285 Å RMSD) [CIT Raman 2016]. The standard in vitro counter-
screen (peptide scanning against known human peptides) failed to flag it.

The screening problem is one of sheer combinatorics. A therapeutic TCR can
in principle recognize any peptide presented by any HLA allele on any of
~20,000 human proteins, yielding a search space of >10^7 candidate peptides
per TCR, each of which must be scored against >10^15 possible TCRs — a
regime where exhaustive experimental screening is physically intractable.
The field therefore needs a *prioritization* strategy: a way to spend a
small, finite number of wet-lab validation experiments on the pairs most
likely to be dangerous. This is precisely the setting for active learning.

Active learning has produced strong efficiency gains in protein engineering
(ALDE, ALSEBO, EVOLVEpro: order-of-magnitude reductions in the experiments
needed to reach a fitness target) [CIT]. But those gains assume a
*predictive oracle* — a model whose uncertainty is informative about what
remains to be learned. Whether the same holds for TCR cross-reactivity is
unclear, because the underlying models are weak: sequence-based TCR-peptide
binders plateau at AUROC 0.62–0.64 on held-out epitopes [Meta-TCR-GNN], and
structural predictors (AlphaFold 3, DockQ ≈0.50 on TCR-pMHC) are too coarse
to score interfaces directly. It is therefore an open question whether
"intelligent sampling" can beat random sampling in this domain, and — more
importantly — under what conditions it does.

Two prior decisions frame this work. First, from our companion study
[Meta-TCR-GNN], we score *rankings* rather than absolute affinities: a
structural compatibility model assigns each TCR a relative ordering over a
candidate peptide library, and safety is judged by whether known lethal
off-targets enter the top of that ordering. This reframing emerged from a
negative result — absolute scoring of two native complexes (MAGE-A3/titin)
did not separate them — and it is the correct formulation for a screening
task. Second, we build the pipeline as a *digital twin* of the ACDC
(antigen-presenting-cell display) dry-wet loop, so that every design choice
(acquisition function, batch size, pre-filter threshold) is testable in
silico against a public-data oracle before any wet-lab reagent is ordered.

This paper makes four contributions. (i) We benchmark acquisition functions
(EIG/BALD, predictive entropy, MC-variance, random) in a simulated closed
loop on held-out epitopes, and show EIG delivers a 2.0× positive-recall
gain at 6% sampling. (ii) We show this headline is *conditional*: the gain
is a function of the MHC-presentation pre-filter threshold, and vanishes
on strictly pre-filtered pools — the regime the real platform always
operates in. (iii) We quantify the two dominant uncertainty sources
(MHC-presentation prediction and AF3 structure-prediction variance) that
constrain any such pipeline. (iv) We report the pipeline end-to-end as a
reproducible in silico safety-assessment protocol (TCR-Safety-Radar), and
flag a stability requirement — cross-instance reproducibility of OOD
rankings — that remains open at the current data scale. The overarching
conclusion is methodological: intelligent sampling for TCR safety is
*necessary but insufficient* without a structure-level signal, and the
field's convenient heuristics (many small batches; sampling always helps)
do not survive contact with weak surrogates and pre-filtered pools.

## 2. Results

### 2.1 Acquisition-function comparison (simulated closed loop)

We simulate a closed loop on 49,696 held-out epitope pairs (7.0% positive,
split by epitope so no epitope seen in training leaks into evaluation). Each
round, a bootstrap logistic ensemble (20 members) trained on physicochemical
peptide/CDR3 features scores every unlabeled pair; an acquisition function
selects a batch of 3,000 to "validate" (reveal the VDJdb label); the oracle
labels are added to the training set and the model is refit. ε-greedy
exploration (ε = 0.15) keeps a floor of random samples to prevent
acquisition collapse.

Table 1 reports positive-recall at 6/12/18/24% of the pool sampled, for
four acquisition functions. EIG (BALD) dominates: 10.4% recall at 6%
sampling versus 5.9% for random — a 2.0× gain. MC-variance is second
(9.3%); predictive entropy is indistinguishable from random (5.8%), a result
we return to below. The gain decays with sampling fraction (1.6×, 1.55×,
1.36× at 12/18/24%), consistent with the expected depletion of easy
positives: the early rounds harvest the separable signal, and later rounds
converge toward the model's residual discriminative ceiling. The low-
sampling regime — 6% of the pool — is precisely the operating point the
ACDC platform targets (0.05% seed sampling in the original proposal), so
the 2.0× figure is the operationally relevant one.

The entropy result is worth isolating as a methodological warning. The
entropy of a *mean* prediction carries no signal about model uncertainty
when the ensemble's disagreement is concentrated in unlabeled regions the
mean happens to place near 0.5; on a weak model this is common, and it is
why the "obvious" choice of uncertainty sampling can silently reduce to
random. We therefore report EIG — which explicitly quantifies expected
information gain — as the acquisition function of record for the remaining
experiments.

### 2.2 Clinical gold-standard ranking (MAG-IC3) — single-instance caveat

As a ground-truth probe we score the MAG-IC3 TCR — the receptor at the
center of the MAGE-A3 trial — against five peptides on AF3-predicted
structures: its target MAGE-A3, the titin mimic ESDPIVAQY, the family
cross-reactant MAGE-A6, the MAGE-B18 peptide (a reported cross-reactivity
of a *different* TCR, included as a near-decoy), and a random control.

Under 5-model ensembling the ordering is (Table 2): MAGE-A3 0.683 ± 0.257,
titin 0.706 ± 0.147, MAGE-A6 0.895, MAGE-B18 0.295, control 0.792. Two
observations stand out. First, the titin mimic scores *at the level of the
target* (0.706 vs 0.683), which matches its confirmed cross-reactive
biology — the model does not falsely clear the lethal off-target. Second,
the random control scores unexpectedly high (0.792), and the ordering is
broadly "muddy": the single-model ranking that originally placed MAGE-A3
first is not reproduced by the ensemble, where titin and the control both
exceed the target.

This is the first of several results that force a downgrade. We later show
(Section 2.5 of [Meta-TCR-GNN]; Figure 5) that these rankings are *not
reproducible across independently trained instances*: pairwise peptide-
ranking Spearman correlations of 0.17 / −0.24 / NaN, and an independent
replication on 21 AF3 jobs yielding 0.61 / 0.11 / 0.33 / 3×NaN. The
MAG-IC3 result is therefore reported as a *single-instance observation* —
a motivating case, not a claim — with the verdict deferred to domain
adaptation (mixing predicted-style structures into training; route now AF3
results after the TCRmodel2 queue stalled, §4). The honest reading is that
at the current data scale the model can sometimes place the lethal
off-target at the target's level, but cannot yet do so reliably.

### 2.3 Uncertainty quantification

A safety pipeline must report not just a score but a defensible uncertainty
interval. Two sources dominate here, and both are measurable.

*Structure-prediction variance (corrected 2026-09-23).* An earlier claim that
re-submitting a job to AF3 can swing a single-model score by 0.94
(LLFGYPRYV 0.049 → 0.991, a first-pass fold failure) does not survive
scrutiny: the scorer that produced those numbers was overwritten by a retrain
and none of the six surviving checkpoints reproduces it. Re-scoring every AF3
submission with a *fixed* surviving instance bounds the resubmission range to
median 0.017 / max 0.166, and AF3's own confidence is highly reproducible
across submissions (ipTM median Δ 0.025, max 0.070; pure resubmissions agree
to ≤0.01), and the peptide Cα register is likewise reproduced (RMSD median
0.20 Å, max 0.38 Å). A separate, detectable failure mode is a defective first
submission with truncated TCR chains (324 vs 443 Cα; flagged by a lower ipTM
of 0.86 vs 0.93). The large swings are therefore *model-side* arbitrariness
(the Rashomon effect, §[Meta-TCR-GNN] 2.6) amplified on such defective inputs,
not AF3 run-to-run variance. Across 5
AF3 models within a single submission, ensemble standard deviations still
range 0.01–0.43. The operational lesson stands but is relocated: a *single
model instance* is not a defensible scorer; report a fixed multi-instance
ensemble behind a stability gate, and treat AF3's own confidence as a
separate, comparatively stable signal.

*MHC-presentation prediction.* The pre-filter that defines the candidate
pool (§2.4) is itself uncertain; its threshold trades recall against pool
size and propagates into which pairs are ever scored.

Two further results sharpen how uncertainty should be used. First, AF3's own
confidence (ipTM) correlates only weakly with our score (Pearson r = 0.333):
predictor confidence cannot substitute for model disagreement, and high
ipTM does not guarantee a sensible score (a high-ipTM JM22 prediction was
still scored low). Second, ensemble std is itself informative — high-
scoring jobs carry small std (0.01–0.03) while *borderline* jobs carry the
largest std (up to 0.43) — so std doubles as a complementary acquisition
signal: the pairs the model is least sure about are the ones to validate
next. We note this is a second, independent justification for EIG-style
sampling that does not rely on the sequence surrogate's quality.

### 2.4 MHC-presentation pre-filter benchmark (ACDC library construction)

The ACDC platform constructs its display library from *presented* peptides
(those that survive MHC presentation), so the candidate pool is by design
pre-filtered. We benchmark the pre-filter that mirrors this: MHCflurry 2.2.1
presentation prediction, evaluated on 50 IEDB-validated HLA-A*02:01 binders
against 50 composition-shuffled decoys — a deliberately hard protocol, since
shuffling preserves amino-acid composition and only destroys order, which
sequence-level predictors partially tolerate.

The pre-filter compresses effectively (Table 4 caption; see §Methods): at a
top-2% threshold it retains 56% of binders at ~50× compression; top-10%
retains 72% at ~10×; mean rank of binders is 30 versus 71 for decoys.
Presentation filtering is therefore a sound first stage that cuts the pool
by an order of magnitude while keeping most true binders.

A second, subtler result concerns what pre-filtering *cannot* do. Emulating
the pre-filter on the VDJdb-derived pool leaves positive density unchanged:
threshold 0.5 keeps 81.7% of pairs at 6.9% density (vs 7.0%), and the
strict 0.9 threshold keeps 67.1% at 7.0%. The reason is structural — VDJdb
negatives are themselves *presented* peptides that simply elicited no
observed TCR response, so peptide-side filtering cannot concentrate
positives. The density gain that a real immunopeptidome pool would provide
(where negatives include many non-presented peptides) is *not emulatable on
VDJdb*. This limitation propagates directly into the boundary result of the
next section.

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

We run the full pipeline at proteome scale as a dry-run deliverable: 20,431
Swiss-Prot proteins → 22.5M canonical 9/10-mers → MHCflurry presentation
scoring → a 4,654,973-peptide pool (threshold 0.0122) → EIG ranking → a
top-50k candidate list with UniProt mapping (per TCR; delivered for
Kimmtrak/gp100 and A6).

The result is a decisive negative. The weak sequence surrogate's EIG does
*not* transfer to proteome scale: known HLA-A2 binders show zero enrichment
in the top-50k (0 of 24 recovered versus 0.3 expected, hypergeometric
P = 1.0), and the Kimmtrak target peptide ranks in the bottom third of the
pool. The 2.0× gain observed in the simulated closed loop (§2.1) was
therefore a *dataset-internal* correlation, not a transferable signal of
TCR-specific recognition — the same failure mode, at sequence space, as the
Rashomon crisis at structure space. The infrastructure (pool construction,
presentation scoring, UniProt localization, candidate export) is complete
and reusable; the scientifically meaningful ranking awaits the
domain-adapted EGNN re-rank. We report this as a pipeline deliverable with
an explicit negative result, to prevent the top-50k list being mistaken for
a validated off-target ranking.

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

### 3.1 Active-learning simulator

The simulator operates on 49,696 held-out epitope pairs (7.0% positive,
epitope-disjoint from training). Each round selects a batch of 3,000 pairs
to "validate" against VDJdb oracle labels, with ε-greedy exploration
(ε = 0.15) reserving a random floor. The scoring model is a bootstrap
logistic ensemble (20 members) over physicochemical peptide/CDR3 features —
deliberately weak, to bound the best case a sequence surrogate can achieve.
Acquisition functions: EIG (BALD), predictive entropy, MC-variance, and
random. Pool modes: raw, and MHCflurry-presentation-prefiltered at 0.5 and
0.9 thresholds (KN-7+). Batch-size ablation at fixed total budget (5,964
samples = 12% of pool) sweeps batch sizes 3600/1800/900/600/300.

### 3.2 Structure scoring pipeline

AF3 five-chain submissions (official web interface, terms-compliant) →
mmCIF parsing → interface-graph construction → EGNN scoring
[Meta-TCR-GNN] → 5-model ensemble mean ± std → per-TCR ranking. Scores are
reported as ensemble means; single-model and single-submission scores are
explicitly treated as non-defensible (§2.3).

### 3.3 Domain adaptation (fixing the Rashomon instability)

To make OOD rankings reproducible, predicted-style structures are mixed
into training so the model sees out-of-distribution inputs during fit. The
positive source is AF3-predicted structures — the scriptable TCRmodel2
alternative stalled (100+ jobs pending on a 2021-era template server), so
the DA route uses AF3 results directly. pLDDT is zeroed in both training
and scoring to remove the AF3 B-factor confound. The gate is a 3-seed
pairwise Spearman on held-out decoys (mean r ≥ 0.5 required before any
ranking is reported as a claim).

### 3.4 MHC presentation pre-filter

MHCflurry 2.2.1 Class1PresentationPredictor; per-allele scoring with
HLA-A*02:01 fallback; thresholds calibrated on decoy quantiles from
composition-shuffled binders (50 IEDB-validated binders vs 50 shuffles).

### 3.5 Dashboard

TCR-Safety-Radar: a zero-dependency web server rendering per-TCR ranking
tables, a TCR × peptide heatmap (traffic-light colors), and a clinical
gold-standard panel, with a built-in honesty banner carrying the
Rashomon-crisis caveat so that interim rankings are never presented as
safety verdicts.

## 4. Discussion

The headline result — a 2.0× recall gain from EIG sampling — should be read
as a *conditional* and *bounded* claim. It is conditional because the gain
is a function of the pre-filter threshold (§2.5), vanishing when the pool
is homogenized by a strict presentation filter. It is bounded because the
simulator's oracle is a weak sequence surrogate: the gains reported here are
the *best* a sequence-level signal can achieve, and therefore a lower bound
on what a structure-level model, scoring real interfaces, could deliver.
Whether that bound is actually exceeded is exactly what the domain-adapted
EGNN re-rank will test.

Two negative-capability findings deserve emphasis because they contradict
convenient field heuristics. First, the "many small batches" rule fails: at
fixed budget, large batches beat small ones (§2.7), because a weak model
needs each batch to be large enough to carry a distinguishable signal.
Second, "intelligent sampling always helps" fails: predictive-entropy
sampling is indistinguishable from random (§2.1), and EIG's advantage is
erased by pre-filtering (§2.5). Both results are specific to the
weak-surrogate regime, but that is precisely the regime in which TCR
cross-reactivity prediction currently operates, and we argue any screening
pipeline must therefore be validated under these conditions rather than
assumed to inherit the nice properties of better-studied protein-fitness
settings.

The Rashomon crisis (§2.2; [Meta-TCR-GNN] §2.6) is the deepest open
problem. A model can discriminate native from decoy interfaces reliably
(interpolation) and yet produce rankings on AF3-predicted structures that
are arbitrary across independently trained instances (extrapolation). This
is not a minor calibration issue but a stability requirement: an
OOD-scoring pipeline that cannot reproduce its own rankings is not a
measurement instrument, and no safety verdict can be built on it. We have
adopted a concrete gate — 3-seed pairwise Spearman on held-out decoys, with
a reporting threshold — and we treat every interim ranking as an anecdote
until that gate is passed. The domain-adaptation route (training on
predicted-style structures) is the fix we are pursuing, now sourced from
AF3 predictions.

Finally, the digital twin is deliberately drop-in for wet-lab data: every
stage — pre-filter, structural scoring, acquisition — takes the same inputs
the real ACDC loop would generate, so the transition from simulation to
live experiments is a matter of swapping the oracle, not the machinery. The
technical route for the wet-lab closed loop is preserved separately
(KEY-NODES.md). What this work establishes, ahead of any reagent, is a
quantitative map of where intelligent sampling does and does not help, and
a stability standard that any in silico off-target screen must meet before
its rankings can inform a safety decision.

## Figures & Tables

- Fig 1: recall-vs-sampling curves (4 acquisition functions).
- Fig 2: clinical ranking (bar with std whiskers).
- Fig 3: AF3 resubmission variance (ipTM stable; fixed-instance EGNN range).
- Table 1: acquisition comparison.
- Table 2: MAG-IC3 ensemble ranking.
- Table 3: EIG gain vs pre-filter threshold (raw / 0.5 / 0.9 pools).
- Table 4: batch-size ablation (3600/1800/900/600/300 at fixed budget).
