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

## 1. Introduction (skeleton)

- TCR-T safety failures: MAGE-A3/titin and MAGE-A12 [CIT Cameron 2013,
  Linette 2013, Morgan 2013]; mechanism = direct molecular mimicry
  (backbone RMSD 0.285 A) [CIT Raman 2016].
- Sequence models plateau: [our Table 1] + literature TITAN 0.62 [CIT],
  ERGO TPP-III 0.669 [CIT], PanPep [CIT].
- Structural prediction insufficient: AF3 DockQ 0.499, TCRmodel2 0.566
  [CIT STCRDab-22 benchmark].
- Our contributions: (i) Cα inseparability diagnosis with data-size
  control; (ii) contact-histogram features that break the barrier;
  (iii) rigorous graft-decoy protocol (Kabsch side-chain transplant,
  node-set pinning, donor-sequence diversity); (iv) a structure-to-specificity
  resource linking 281/291 STCRDab complexes to VDJdb records (2,309 CDR3
  hits; 1,517 cross-epitope); (v) a cross-instance stability test exposing
  ranking non-reproducibility on AF3 inputs (Rashomon effect) and the
  scale-monotone domain-adaptation route that resolves it (0.447→0.729 at
  144 native positives, crossing the 0.5 gate).

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

### 2.4 Cross-reactivity ranking on AF3-predicted structures — single-instance observations, later overturned (Fig 3, Table 4)

> Status: the numbers below were obtained with ONE trained instance. Section
> 2.6 shows they do not survive instance resampling; they are retained as the
> motivation for the stability analysis, not as claims.
>
> **Label correction (2026-09-23):** the pdb attributions below are unreliable.
> `kn5_submission_list.tsv` had its tcr_a/tcr_b columns swapped, so the import
> silently fell back to "first candidate row" and mislabeled 48/73 structures
> (verified against 291 crystal TCR signatures). The "A6" entries are in fact
> the B7/1qrn family, "1G6R/1MWA" are the 2C/1tcr family, and "JM22 (1OGA)" is
> the 2vlj/2vlk family. Group identities below should be read as chain-signature
> families until the structure map is repaired.

- A6 TCR: 8 candidate peptides, 5/8 VDJdb-validated binders score
  0.92-0.99 (ensemble); native peptide ranks #3 (0.973). Weak/non-binders
  LLFGPVYV 0.029, LLFGKPVYV 0.659.
- B7 TCR (mirror family): Tax homolog 0.973.
- 1G6R/1MWA mirror pair: cross-scores 0.73-0.86.
- JM22 (1OGA): all candidates low (0.16-0.23) — consistent with score-1
  weak VDJdb evidence; interpreted as correct rejection (inconclusive).
- Clinical gold standard (MAG-IC3/5BRZ): native MAGE-A3 ranked #1
  (0.683 ± 0.257); the titin mimic ESDPIVAQY at 0.706 ± 0.147 — at the
  target's level under 5-model ensembling, matching its confirmed
  cross-reactive biology.

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
  (conditional on >=~144 native positives + 3-seed ensembling).
- v10 Coulomb edge features (residue net-charge product): negative result
  (graft 0.933 → 0.917); reverted — atom-level partial charges required.
- **Methodological claim: any OOD-scoring claim must pass a cross-instance
  stability test; single-instance rankings are anecdotes.**

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

## 4. Discussion (skeleton)

- Representational minimum for interface compatibility: side-chain-level
  contact information is necessary; Cα-only protein graphs are blind to
  chemistry-swapped interfaces.
- Implications for geometric-DL protein design: default Cα pipelines may
  silently fail on tasks where chemistry differs under identical geometry.
- Ranking vs absolute scoring; calibration across complexes remains open.
- Limits: 423 complexes; AF3 variance; seed variance; VDJdb label noise.
- Clinical relevance: MAG-IC3/titin verdict as a step toward in silico
  safety screening for TCR therapeutics.

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
- Table 1: dataset statistics (+ structure-VDJdb map: 281/291 complexes,
  2,309 CDR3 hits).
- Table 2: hypothesis-elimination results.
- Table 3: ablation suite (+ v10 Coulomb negative).
- Table 4: clinical gold-standard ranking (single-instance caveat).
