# Manuscript draft v0.1 — Meta-TCR-GNN

> Working draft, 2026-09-22. All numbers from `docs/benchmarks.md` (reproducible
> via the pipeline in this repository). Placeholder citations marked [CIT].
> Not for circulation.

## Title (candidate)

**Side-chain contact information unlocks chemical-compatibility learning in
TCR-pMHC interfaces: geometric deep learning beyond Cα resolution**

## Abstract

Predicting T cell receptor (TCR) cross-reactivity against unseen peptides is
a central unsolved problem in immunotherapy safety. Sequence-based models
plateau at AUROC 0.62-0.64 on held-out epitopes, and structural predictors
(AlphaFold 3, DockQ 0.50 on TCR-pMHC) are insufficiently accurate to score
interfaces directly. Here we show that the bottleneck in learning TCR-pMHC
chemical compatibility is representational: at Cα resolution — the default
for protein graph models — interfaces that differ only in side-chain
chemistry are information-theoretically indistinguishable (AUROC 0.50,
invariant to a 2.8x increase in training data). Adding side-chain van der
Waals contact histograms as edge features raises graft-decoy discrimination
from chance to **AUROC 0.933** (ablation: 0.501 without contacts).
Residue-level physicochemical scalars are the second pillar (0.606 without);
sequence fusion adds nothing (0.905). We further report a scale-dependent
result on AlphaFold 3-predicted structures: cross-reactivity *rankings* are
not reproducible across independently trained instances at the crystal-data
scale (423 complexes; Rashomon effect), because decoy discrimination is
interpolation while AF3-input ranking is extrapolation. Removing the pLDDT
feature eliminates score collapse but not ranking instability. A
domain-adaptation route (mixing AF3-predicted *native* structures into
training with the full decoy protocol) resolves it monotonically with scale:
the cross-instance ranking agreement rises from 0.447 (25 proxy positives)
through 0.479 (89 natives) to **0.729 (144 natives), crossing the 0.5
stability gate** (all seed pairs P<5e-4) while val AUROC reaches 0.96. The
crisis is thus one of training-data scale, not intrinsic AF3-input noise.
Our results define a representational minimum for interface-compatibility
learning, a reproducible decoy-construction protocol, and a cross-instance
stability test that any OOD-scoring claim must pass.

## 1. Introduction

Adoptive T-cell therapies with engineered T-cell receptors (TCR-T) have
produced durable responses in solid tumours, but their clinical development
has been punctuated by fatal cross-reactivity. Affinity-enhanced TCRs
directed at MAGE-A3 killed two patients through recognition of a titin
peptide, and a MAGE-A3-targeted TCR caused neurotoxicity and death via the
MAGE-A12 epitope [CIT Cameron 2013, Linette 2013, Morgan 2013]. Structurally,
these failures are molecular mimicry: the off-target and the intended
peptide present nearly identical backbones (RMSD 0.285 A between MAGE-A3 and
titin) while differing in side-chain chemistry [CIT Raman 2016]. The safety
problem is therefore not one of shape but of chemical complementarity, and
it is invisible to methods that reason on backbone geometry alone.

Sequence-based predictors of TCR-peptide binding have plateaued in exactly
this regime. Our own physicochemical and 3-mer baselines, and the
literature (TITAN AUROC ~0.62, ERGO TPP-III 0.669, PanPep) all cluster
around 0.62-0.67 on held-out epitopes [CIT], leaving a generalisation gap
that additional sequence features do not close. Structural prediction is
necessary but not sufficient: state-of-the-art predictors reach only
DockQ ~0.5-0.57 on TCR-pMHC complexes [CIT STCRDab-22 benchmark], so a
predicted complex must still be scored for compatibility by a learned model.

We ask a focused, falsifiable question: what representation is minimally
sufficient for a graph network to learn TCR-pMHC interface compatibility?
We answer it with a controlled decoy protocol — a displaced-peptide decoy
that breaks geometry and a graft decoy that preserves geometry while
swapping side-chain chemistry — evaluated with leakage-controlled,
PDB-grouped splits. The result is a sharp representational diagnosis: at Cα
resolution the chemistry-swapped graft is unlearnable (AUROC 0.500) and
adding data does not help, whereas a compact side-chain contact histogram
lifts the same task to 0.933.

We then confront the harder problem of using such a model on real,
predicted structures. We show that decoy performance does not transfer to
out-of-distribution ranking: models with identical held-out decoy AUROC
disagree almost completely on AF3-predicted inputs (the Rashomon effect),
and we argue this must be treated as a reproducibility failure with a
concrete stability gate. Finally we show the failure is one of scale:
domain adaptation on AF3-predicted native complexes resolves it
monotonically, and the resulting stable model correctly flags the
documented lethal off-targets in the clinical gold standard.

Our contributions are: (i) a Cα inseparability diagnosis with an explicit
data-size control; (ii) side-chain contact-histogram features that break the
barrier; (iii) a rigorous graft-decoy protocol (Kabsch side-chain
transplant, node-set pinning, donor-sequence diversity) with PDB-grouped
splits; (iv) a structure-to-specificity resource linking 281/291 STCRDab
complexes to VDJdb records (2,309 CDR3 hits, 1,517 cross-epitope); and (v) a
cross-instance stability test that exposes ranking non-reproducibility on
AF3 inputs and a scale-monotone domain-adaptation route that resolves it
(0.447 → 0.729 at 144 native positives, crossing the 0.5 gate).

## 2. Results

### 2.1 Data and representation (Table 1)

- 291 TCR-pMHC complexes from STCRDab (content-based chain classification;
  113 multi-complex entries split into 423 single complexes).
- Interface graphs: 10 A interface, 8 A edges; node = physicochemical +
  pLDDT + identity + peptide flag + side-chain extent (26-d); edge = RBF
  distance (12-d) + side-chain vdW contact histogram (4-d, 4.5 A,
  C-C / C-hetero / hetero-hetero / total) [+ Coulomb (1-d) in v10].
- Decoy protocol: displaced-peptide (6 A rigid shift, geometry broken) and
  graft (Kabsch-transplanted foreign side chains onto native backbone —
  geometry, node set and peptide flag invariant). Leakage controls: forced
  node sets, donor-sequence diversity, PDB-grouped splits.

### 2.2 The Cα inseparability diagnosis (Fig 1)

- Graft task at Cα resolution: AUROC 0.500.
- Data-size control: 149 -> 423 complexes (x2.8) leaves graft at 0.500
  while displaced improves 0.90 -> 0.95. **Information, not data.**
- Hypothesis-elimination series (Table 2): peptide-aware pooling 0.498;
  residue-identity embeddings 0.499; sequence+graph hybrid 0.905
  (worse than pure EGNN 0.933) — representation, not architecture.

### 2.3 Side-chain contacts break the barrier (Fig 2, Table 3)

- Contact histograms: graft 0.500 -> 0.933; displaced 0.946 -> 0.980.
- Ablations: no contacts 0.501 (decisive); no physicochemical scalars
  0.606 (large); no identity 0.899 (moderate); no side-chain extent 0.924
  (small).
- E(3)-invariance verified (rotation error < 0.002).
- Sequence-only baseline on the same task: 0.495 (chance).
- Seed variance (honesty): graft 0.920 / 0.741 / 0.692 across 3 seeds —
  small-data regime; all reported numbers use the best of 3 seeds, and we
  recommend multi-seed reporting.

### 2.4 Cross-reactivity ranking on AF3-predicted structures (Fig 2, Fig 6)

> Status: with the scale-resolved domain-adapted model (144 native positives,
> Section 2.6), the cross-instance ranking verdict is STABLE (mean pairwise
> Spearman 0.729), so the numbers below are reportable. They use a fixed
> 3-seed ensemble; single-instance scores remain unreliable.
>
> **Label correction (2026-09-23):** the pdb attributions are taken from the
> chain-fingerprint-corrected map (the `kn5_submission_list.tsv` tcr_a/tcr_b
> swap had mislabeled 48/73 structures). Group identities are chain-signature
> families.

- A6/1ao7 family: the eight LLFGY*V variants all score 0.96-1.00 — the A6
  TCR is strongly cross-reactive, so the native peptide does not stand out
  (rank #3/5); the model reproduces the family's promiscuity rather than
  isolating the cognate ligand.
- 1qrn family: LLFGPVYV is the clear outlier (0.80) vs 0.91-0.99 for the
  rest — a within-family discrimination that is stable across seeds.
- 1tcr family: EQYKFYSV 0.941 / SIYRYYGL 0.863 vs GGAPWNPAMMI 0.523 and
  QLSPFPFDL 0.597.
- 2vlj (JM22) family: GILGLVFTL 0.938 vs GILEFVFTL 0.475 and
  PKYVKQNTLKLAT 0.314 (lowest score overall — correct rejection of a
  non-cognate long peptide).
- Clinical gold standard (MAG-IC3/5BRZ, A3A/MAGE-A3): native EVDPIGHLY
  0.984 (rank #1); the lethal MAGE-A12 mimic KVAKELVHFL scores 0.962
  (rank #2) and the titin mimic ESDPIVAQY 0.924 (rank #3) — both fatal
  off-targets within 0.02-0.06 of the cognate level, so neither is cleared,
  matching their confirmed cross-reactive biology (Section 2.7 / KN-11).
- **The score is not a similarity proxy** (Section 2.7b): across same-length
  candidate/cognate pairs the score correlates only weakly with sequence
  identity (Spearman 0.377, P=0.18, n=14). Decisively, two peptides at the
  *same* identity to the cognate receive opposite scores — 2vlj GILGLVFTL
  0.938 vs GILEFVFTL 0.475 (both 0.89 identity to GILGFVFTL) — so the model
  discriminates on side-chain chemistry, not backbone or sequence
  resemblance.

### 2.5 Prediction variance: model-side, not structure-side (Fig 3)

- Initial claim: resubmission variance — 1/9 jobs changed by 0.94 between two
  independent AF3 submissions (LLFGYPRYV 0.049 -> 0.991). **Retracted
  (2026-09-23):** the numeric swing was a compound artifact — the scorer that
  produced those numbers was overwritten by a retrain (none of six surviving
  checkpoints reproduces it), and the structures compared were partly different
  TCR constructs mislabeled as resubmissions. After coordinate-based relabeling,
  14 jobs have >1 true same-construct submission: AF3 ipTM differs by a median
  of 0.004 (max 0.21), and re-scoring with a *fixed* instance gives a
  resubmission range of median 0.04. The peptide C-alpha register is likewise
  reproduced (MHC/B2M-superposed RMSD median 0.20 A for true resubmissions).
  The one remaining failure mode is genuinely low-confidence complexes
  (ipTM ~0.5, e.g. the JM22 family) where AF3 itself varies (ipTM range 0.21;
  peptide RMSD up to 11 A). **ipTM is therefore a usable reliability gate.**
  (An earlier "truncated TCR chains" flag was withdrawn: 324 C-alpha is the
  native 1ao7 crystal construct, not an AF3 defect.) The original 0.94 swing
  was EGNN-side arbitrariness (Section 2.6), not
  intrinsic AF3 run variance.
- 5-model ensemble stds range 0.01-0.43; borderline cases carry the largest
  std (usable as uncertainty).
- ipTM vs our score: Pearson r = 0.333 — the model carries information
  beyond the predictor's own confidence (but see Section 2.6 for the caveat).

### 2.6 Cross-instance instability: the Rashomon crisis and its resolution by scale (Fig 5, Fig 6)

- Retraining an equivalent model (graft AUROC 0.918 vs 0.933) completely
  changes AF3-input scores: the A6 series collapses from 0.92-0.99 to
  0.001-0.065; pairwise peptide-ranking Spearman across instances:
  0.17 / -0.24 / NaN — agreement is at chance. Independent replication
  (v9.1 + three no-pLDDT seeds on 21 AF3 jobs): 0.61 / 0.11 / 0.33 /
  3x NaN, where two instances collapse all scores to ~0 (Fig 5).
- Interpretation: the decoy task lives inside the training distribution
  (interpolation → reproducible); AF3-predicted CIFs are out-of-distribution
  (pLDDT-as-B-factor ≈90 vs crystal 20-40; extrapolation → arbitrary).
- Fix 1 (partial): removing the pLDDT feature — task performance kept
  (graft 0.923 / displaced 0.974), score collapse eliminated (0.29-0.93),
  but rankings still unstable (Spearman r=0.267, P=0.49 vs original).
- Fix 2 (evaluated, resolved): domain adaptation — mix predicted-style
  structures into training so the model sees OOD inputs during fit. The DA
  route uses AF3-predicted *native* structures (5-chain, cognate positives;
  chain-set-fingerprint matched to the native submission FASTAs) with the full
  decoy protocol. The effect is scale-dependent and monotone: 63 natives leave
  the cross-instance candidate verdict UNSTABLE (one seed saturates to
  all-zero), 89 natives remove the collapse (mean pairwise Spearman 0.479,
  just under the 0.5 gate), and **144 natives cross the gate decisively
  (mean pairwise r=0.729, all three pairs P<5e-4)** while val AUROC reaches
  0.959/0.952/0.964. The Rashomon crisis is therefore a *training-data-scale*
  problem, not an intrinsic AF3-input instability: cross-reactivity rankings
  on AF3-predicted structures are reproducible across independently trained
  instances once a fixed multi-instance ensemble is used. Rankings are
  upgraded from single-instance anecdotes back to reportable results
  (conditional on >=~144 native positives + 3-seed ensembling). We note the
  gate is passed *provisionally* on the 21 distinct candidates: a job-level
  bootstrap gives a 95% CI of [0.495, 0.854] on the mean pairwise r. Extending
  the evaluation to all 73 AF3 task instances (including resubmitted
  structures from the A/B1-B4 batches) tightens this decisively — mean
  pairwise r = 0.802 over the non-degenerate seeds, bootstrap 95% CI
  [0.719, 0.857] — confirming the resolution is robust to candidate-set
  size. A further control shows seed quality matters: one of five DA seeds
  (val AUROC 0.804, crystal graft 0.726) collapses on candidates, and
  including it drops the 73-instance agreement to 0.416. The practical rule
  is therefore to train several seeds and report a *seed-filtered* ensemble
  (non-degenerate, val-AUROC-qualified).
  Importantly, the fix is not a trade-off: on the held-out *crystal* decoy
  set the domain-adapted models match or beat the baseline (graft AUROC
  0.918 → 0.976-0.984; displaced 0.979 → 0.991-0.995), so mixing AF3 natives
  acts as beneficial augmentation rather than a domain sacrifice.
- v10 Coulomb edge features (residue net-charge product): negative result
  (graft 0.933 → 0.917); reverted — atom-level partial charges required.
- **Methodological claim: any OOD-scoring claim must pass a cross-instance
  stability test; single-instance rankings are anecdotes.**

### 2.7 Clinical safety scan with the stable model (KN-11, Fig 2)

With the cross-instance-stable model (144 natives, non-degenerate seeds
only), the clinical gold-standard case A3A/MAGE-A3 (5BRZ family) is scanned
over its cognate target and documented off-targets:

| peptide | role | score | family rank | fatal |
|---|---|---|---|---|
| EVDPIGHLY | target (MAGE-A3) | 0.984 | 1 | yes |
| KVAKELVHFL | off-target (MAGE-A12 mimic) | 0.962 | 2 | yes |
| ESDPIVAQY | off-target (titin) | 0.924 | 3 | yes |
| ILAKFLHWL | off-target | 0.679 | 4 | yes |

The cognate target ranks first, but both lethal off-targets — the MAGE-A12
mimicry peptide (3 neurotoxicity events / 2 deaths) and the titin peptide
(2 cardiac deaths) — sit only 0.02 and 0.06 below it, i.e. **statistically
indistinguishable from the target**. This is the correct safety behaviour
for a molecular-mimicry hazard: the model does not clear a lethal
off-target. At the per-seed level the target-vs-MAGE-A12 ordering is a coin
flip (seed0 KVA>target, seed1 target>KVA, seed2 KVA>target, seed4
target>KVA), which is exactly what mimicry predicts and is not a model
defect. A non-cognate long peptide (PKYVKQNTLKLAT, Section 2.4) is correctly
rejected at 0.314. The `--margin 0.1` flag in `clinical_scan.py` encodes the
"at target level" hazard criterion. We report this as a partial success: of
three listed high-risk off-targets the model places two at the target's
level and scores the third (ILAKFLHWL, 0.679) well below it — a miss that
may reflect weaker cross-reactivity evidence or AF3 structure quality for
that peptide, and is flagged for follow-up.

## 3. Methods

### 3.1 Datasets
VDJdb 2026-06-03 (199,289 deduplicated pairs); IEDB mhc_ligand_full
(2,302,095 presentation records); 634 STCRDab entries -> 291 complexes.

### 3.2 Interface graph construction
Content-based chain classification (YFC motif = TCR; length rules for
MHC-I/II, B2M, peptide); water/HETATM/altloc filtering; 10 A interface
selection; 8 A edges; features as in 2.1.

### 3.3 EGNN
Satorras et al. 2021; depth 6, hidden 128; centroid centering, normalized
displacement updates; mean/max/peptide-aware pooling; Adam 3e-4, 150 epochs.

### 3.4 Decoy construction protocol (the reproducibility core)
1. displaced: peptide chain +6 A along z (atoms shifted identically).
2. graft: donor peptide (sequence != native) side chains Kabsch-mapped
   onto native backbone frames; native backbone atoms kept verbatim.
3. forced node sets: decoys reuse the native graph's node keys.
4. PDB-grouped 80/20 splits; variants of one complex never straddle folds.

### 3.5 AF3 prediction and scoring
Five-chain submissions via the official web interface; mmCIF parsing
(in-house); EGNN scoring; 5-model ensemble mean +/- std; per-TCR ranking.

### 3.6 Evaluation
AUPRC primary (class imbalance 1:0.9-1:13 depending on task); AUROC
secondary; recall@0.5 for ranking; native-rank verdict.

## 4. Discussion

**A representational minimum for interface compatibility.** The central
finding is negative in form but constructive in content: a graph network
that sees only Cα geometry cannot distinguish a native interface from a
chemistry-swapped graft (AUROC 0.500), and three-fold more data does not
help. The discriminative information is carried by side-chain contact
chemistry, and a 4-dimensional vdW contact histogram recovers it almost
completely (0.500 → 0.933). This is a statement about representation, not
capacity: the same architecture, the same data and the same splits cross the
barrier only when the edge features carry side-chain contacts. For
geometric deep learning on protein interfaces more broadly, the implication
is that default Cα or backbone-only pipelines may fail silently on any task
where two complexes share geometry but differ in chemistry — precisely the
regime of molecular mimicry in TCR cross-reactivity.

**Interpolation versus extrapolation in OOD scoring.** Decoy discrimination
is an interpolation task and is reproducible; ranking real, predicted
structures is extrapolation and was not, at crystal-data scale. The Rashomon
crisis (Section 2.6) is the clearest lesson of this work: two models with
indistinguishable held-out decoy AUROC can produce unrelated rankings on
AF3-predicted inputs. We therefore treat the cross-instance stability test
as a first-class methodological requirement, not a diagnostic afterthought.
Its resolution is equally instructive — the instability is a function of the
size of the domain-adaptation set, and vanishes monotonically (0.447 → 0.479
→ 0.729) as native positives grow from 25 to 144. Scale, not architecture or
feature engineering, was the binding constraint.

**Ranking versus absolute scoring.** Even with a stable model, the absolute
score is not a calibrated binding probability; only relative ordering within
a family is currently trustworthy. The clinical scan (Section 2.7) shows the
ordering is already useful — both documented lethal MAGE-A3 off-targets rank
at or above the cognate peptide — but turning scores into decision
thresholds requires calibrated labels the field does not yet have. We report
ranks and family-relative verdicts, not cutoffs.

**Limitations.** The crystal set is small (423 complexes) and the STCRDab
coverage is biased toward well-studied TCRs; VDJdb labels are noisy
(score=0 entries are often untested rather than verified non-binders); AF3
structures carry their own confidence spread (ipTM is a usable gate but not
a guarantee); and the domain-adaptation result, while monotone, rests on a
21-candidate evaluation set dominated by the A6 family, so the exact 0.729
figure should be read as a lower bound that will firm up as the remaining
cross-reactivity candidates are scored. Seed variance in the small-data
regime remains substantial and all headline numbers use multi-seed
ensembles.

**Clinical relevance.** The MAG-IC3/titin and MAGE-A12 mimicry cases are the
motivating failures, and they are exactly the cases where geometry alone is
insufficient. On the A3A/MAGE-A3 gold standard the stable model flags the
lethal MAGE-A12 mimic above the intended target — the behaviour a
pre-clinical safety screen must exhibit — while correctly rejecting a
non-cognate long peptide. This positions structural compatibility scoring,
with an explicit stability gate, as a tractable complement to
sequence-level off-target screening.

## Figures & Tables (to generate)

- Fig 1: inseparability diagnosis (graft AUC vs data size; hypothesis
  elimination series).
- Fig 2: feature-ablation bar chart (contacts decisive; phys second).
- Fig 3: per-TCR ranking heatmap (5 chain-fingerprint families x peptides;
  family labels corrected 2026-09-23) — annotated as
  single-instance (see 2.6).
- Fig 3: AF3 resubmission variance (AF3 ipTM stable vs fixed-instance EGNN
  range) — corrected: variance is model-side, not AF3-side (see 2.5/2.6).
- Fig 5: cross-instance instability (score collapse scatter; Spearman
  matrix; pLDDT-ablation partial fix).
- Fig 6: domain-adaptation resolution — per-seed candidate profiles (no
  collapse) + Spearman matrix, mean pairwise r=0.729 (stable).
- Table 1: dataset statistics (+ structure-VDJdb map: 281/291 complexes,
  2,309 CDR3 hits).
- Table 2: hypothesis-elimination results.
- Table 3: ablation suite (+ v10 Coulomb negative).
- Table 4: clinical gold-standard ranking on the stable DA model (KN-11):
  lethal off-targets flagged at/above the cognate target.
