# What an audit teaches: pre-registered, family-disjoint benchmarking of structure-based TCR representations — a reconstruction case study

**Status**: preprint draft v0.2 (2026-10-21). All numbers derive from committed,
hash-pinned artifacts. This manuscript claims **no binding prediction and no
clinical utility**. Sections marked [TCRen] are pending the external-baseline
run documented in Sec. 6.2.

---

## Abstract

Structure-based models of TCR–pMHC interaction are increasingly common, yet the
field's validation practice is under scrutiny: recent structure-based tools
have themselves reported that structural inputs have "not demonstrated
significant gains" over sequence methods, and reviews of TCR binding models
repeatedly identify validation-design biases. Here we present a case study in
what rigorous validation actually requires, and the price of not doing it. Our
pipeline's previous quantitative results were fully withdrawn after an
external audit of 28 items and a second delivery review of 39 items (5 fatal)
found that input-encoding, splitting, QC-gating and receptor-identity defects
survived into the code paths that produced the reported numbers — including a
clinical case where an engineered ImmTAC reagent (MAG-IC3) had been scored and
reported as the lethal clinical construct (a3a). We report (i) the
reconstruction: byte-level, hash-pinned repair of the pipeline, with
identity-verified inputs, atomic screening, provenance receipts and
family-aware splitting; (ii) a pre-registered representation-quality benchmark
whose labels are verified structural metadata — receptor identity — on 331
interface graphs (160 receptors, 110 families) from a hash-manifested corpus.
An EGNN embedding trained with a supervised-contrastive objective retrieves the
correct receptor at rank-1 accuracy 0.166 versus a per-query chance level of
0.0123 (13.4×; minimum 7.2× over five family-grouped cross-validation folds;
226 queries), with family classification at 3.2× chance, no crystal/predicted
domain leakage (0.496 vs 0.5 chance), and a same-vs-different-receptor pair
AUROC of 0.767 — below our pre-registered 0.9 bar, reported as a limitation.
We also report a numerical stability finding: the equivariant
coordinate-update branch of the EGNN diverges under the contrastive objective
on this corpus (bifurcation-edge dynamics), documented in a protocol
amendment. Two external controls bound the result: official-weight ERGO-I
performs at chance on 10,000 measured functional observations
(per-peptide macro-AUROC 0.511 vs 0.500 baseline), and [TCRen — pending]. We
argue that the primary transferable contributions are not the numbers but the
protocol elements — pre-registration, family-disjoint splits, per-query chance
calibration, identity verification against role-annotated references, and
compulsory negative reporting — which the field should adopt before
interpreting any structure-derived score.

## 1. Introduction

### 1.1 The reproducibility problem in TCR–pMHC modeling

Binding prediction for T cell receptors is dominated by sequence-based
supervised models (TITAN, NetTCR, ERGO, pMTnet, MixTCRpred and many others).
Structure-based approaches promise access to interface geometry and chemistry
(TCRen, STAG, NetTCR-struc, SageTCR, t2pmhc), yet the field's own assessments
are candid: structure inputs have not yet shown clear gains over sequence
methods, and published reviews repeatedly identify validation-design biases —
dataset leakage, optimistic splits, missing negatives, and performance reported
on protocols that cannot distinguish learned binding rules from learned
artifacts of data construction.

### 1.2 This paper

This paper is a case study of that problem, conducted on ourselves. We do not
present a new architecture (ours is a standard EGNN), and we do not claim
binding prediction. We report:

1. the **audit**: how 28 + 39 externally identified defects invalidated every
   quantitative claim of a TCR–pMHC structural pipeline (Sec. 2);
2. the **reconstruction**: what byte-level repair and identity verification
   actually require (Sec. 3);
3. the **benchmark**: a pre-registered, family-disjoint, chance-calibrated
   representation-quality protocol, run exactly as registered, including its
   partial failure against our own criteria (Sec. 4–5);
4. the **controls**: an official-weight published tool at chance on measured
   functional data, and a second structure-based tool on the same corpus
   [TCRen — pending] (Sec. 6);
5. the **recommendations**: which protocol elements the field should adopt
   (Sec. 7).

## 2. The audit: 67 findings that invalidated everything

### 2.1 Round one (28 items)

Label semantics (VDJdb evidence confidence treated as binding labels,
fabricating 185,117 "negatives"; IEDB endpoints conflated), model-input
defects (three-letter residue names mapped by first character, collapsing
ARG/ASN/ASP to A; peptide-node pooling reading an amino-acid identity column
instead of the peptide flag; batch-wide coordinate normalization making a
graph's features depend on its batchmates), splitting without complex
grouping, decoy construction with stale centroids and fixed-axis displacement,
and clinical reference errors (wrong PMIDs; wrong structure attribution).

### 2.2 Round two (39 items; 5 fatal)

An independent delivery review found, among others: two delivered result tables
contradicting each other for the same jobs with reversed qualitative
conclusions; interpreted score differences (0.02–0.06) an order of magnitude
below the scores' own dispersion (median per-job SD 0.2512); a clinical index
case where the structure of the ImmTAC reagent MAG-IC3 was scored and labelled
as the clinical a3a construct that caused patient deaths; QC flags bypassed in
4 of 7 delivered clinical rows; and at least eight "fixed" defects still
present in the exact code paths that had produced the reported numbers.

### 2.3 Why this happens

The audit's consistent finding was a mismatch between methodological intent
(correct, well documented) and engineering delivery (fixes not rerun, reruns
not delivered, deliveries not reconciled). The same review praised the honesty
of the self-audit while concluding that no quantitative result survived. This
is the failure mode our protocol elements (Sec. 4) are designed to prevent.

## 3. Reconstruction

- **Identity verification**: chain roles are never inferred from chain letters,
  lengths, or motifs. Predicted structures are annotated only by exact sequence
  matching against role-annotated reference records (alpha/beta/mhc/b2m).
- **Input hygiene**: every input file is hash-pinned into manifests; releases
  carry fingerprint receipts; scoring requires per-model QC, release and
  checkpoint hashes; manifests reject duplicates and announce legacy
  swapped-column files.
- **Construction screening**: synthetic controls pass the same heavy-atom and
  van-der-Waals overlap screens as native structures; decoy builders raise on
  length mismatch, incomplete backbones, or degenerate geometry instead of
  silently falling back.
- **Splitting and statistics**: family-aware splits; family-stratified Fisher
  pooling; bootstrap clustering by receptor, not PDB prefix.
- **Forensic verification of the reconstruction**: five-seed replay of the
  withdrawn run reproduces 938/940 job×seed means to four decimals; the
  347-item pending panel is byte-identical to its first-commit state (98 exact
  αβ receptor pairs, 45 MHC, 11 β2m sequences). These confirm the repaired
  pipeline's inputs are what they claim to be.

## 4. The pre-registered benchmark

### 4.1 Task design rationale

Receptor identity is the one label that is simultaneously (a) consequential —
any binding model must at minimum distinguish receptors; (b) verifiable — it is
structural metadata, not experimental binding data; and (c) abundant — 55 of
160 receptors in our corpus have multiple structures. The benchmark therefore
asks the minimal question whose answer any structure-based TCR model must get
right before its scores can be interpreted: **does the representation encode
who the receptor is?**

### 4.2 Corpus

331 interface graphs: 168 atom-screened native structures (145 crystal, 23
AF3-predicted) from the verified inventory, plus 163 predicted candidates
admitted only on exact αβ sequence match to role-annotated references (16
excluded and recorded). 160 receptors; 110 families; peptide lengths 8–13;
confidence features masked everywhere. No binding labels exist anywhere in the
corpus.

### 4.3 Protocol (registered before training)

Model: EGNN (6 layers, hidden 128, feature dims 26/16 — the scoring-path
schema) + 64-dim projection head; supervised-contrastive objective (positives =
same receptor), temperature 0.1; Adam (3e-4, wd 1e-4); 60 epochs; batch 16.
Splits: 5-fold CV grouped by receptor family. Tasks: (A) receptor retrieval
(rank-1, MRR, pair AUROC, per-query chance); (B) family k-NN classification;
(C) domain diagnostic. Success criteria, fixed in advance: rank-1 ≥ 2× chance;
AUROC ≥ 0.9; family top-1 ≥ 2× chance; same-seed reproducibility. Failures are
reported as-is.

### 4.4 Protocol amendment: a stability finding

The equivariant variant (coordinate-update branch on, matching the scoring
path) completed fold 0 (rank-1 0.238, AUROC 0.710) then diverged
deterministically in fold 1 (non-finite loss at epoch ~13–15). Isolated
replays with identical seeds diverged or survived depending on CUDA reduction
noise — the dynamics sit at a bifurcation edge; the mechanism is contrastive
gradients flowing back through the pooled coordinate-norm feature into the
coordinate-update branch, amplifying distances and messages. We did not tune
around the failure; we adopted the pipeline's sanctioned invariant variant
(coordinate updates off) for the benchmark and record the amendment and both
variants' outputs verbatim. The instability of the equivariant branch under
contrastive training is itself a finding.

## 5. Results

### 5.1 Receptor retrieval (Task A)

| Fold | rank-1 | chance | ratio | MRR | AUROC | queries |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0.2857 | 0.0105 | 27.1× | 0.3814 | 0.7580 | 42 |
| 1 | 0.1167 | 0.0162 | 7.2× | 0.2454 | 0.7518 | 60 |
| 2 | 0.2532 | 0.0194 | 13.0× | 0.3672 | 0.7589 | 79 |
| 3 | 0.1250 | 0.0101 | 12.4× | 0.2389 | 0.7076 | 24 |
| 4 | 0.0476 | 0.0055 | 8.7× | 0.1895 | 0.8569 | 21 |
| mean | **0.1656** | **0.0123** | **13.4×** | **0.2845** | **0.7666** | 226 |

### 5.2 Family classification (Task B)

Top-1 accuracy 0.224 (chance 0.071, 3.2×); k=3: 0.246; k=5: 0.224.

### 5.3 Domain diagnostic (Task C)

Crystal-vs-predicted nearest-centroid accuracy 0.496 (chance 0.5): no
detectable domain signature under confidence masking.

### 5.4 Verdict against the pre-registered criteria

Two of three quantitative criteria met (13.4× retrieval; 3.2× family
classification); the AUROC criterion (≥0.9) not met (0.767) and reported as a
limitation rather than redefined; reproducibility met for the invariant
variant with the equivariant divergence recorded as an amendment. The
interpretation we defend is narrow: the reconstructed pipeline demonstrably
encodes receptor identity in structure-derived representations, under a
family-disjoint protocol with no domain leakage — the necessary, but not
sufficient, precondition for structure-based binding claims.

## 6. External controls

### 6.1 ERGO-I on measured functional data (committed, 2026-10-21)

Official ERGO-I weights (fixed McPAS LSTM checkpoint, git-blob-verified) on
10,000 measured MIRA functional observations, under a pre-registered protocol
(input-order controls, per-peptide metrics, no threshold search):
per-peptide macro-AUROC 0.5112 vs 0.5001 for the SPR-derived baseline;
macro-AP 0.1166 vs 0.1197. Both are at chance for this endpoint. We do not
claim a cross-task method comparison; we report the bound on what
sequence-only inputs carry here.

### 6.2 TCRen on the verified corpus [pending]

[TCRen — see Sec. 8 reproducibility for the protocol. To be filled when the
run completes: cognate-vs-candidate ranking per receptor on the verified
corpus, reported either way.]

## 7. What the field should adopt

The transferable output of this case study is not the 13.4× number — a working
structure model should be able to identify receptors — but the protocol
elements whose absence caused our own failures and which we now demonstrate are
operationally feasible:

1. **Pre-registration with fixed success criteria**, evaluated as-is;
2. **family-disjoint splitting** (receptor families never span train/test);
3. **per-query chance calibration** (performance reported against the
   composition-correct chance level, not against 0.5 or 1/N);
4. **identity verification against role-annotated references**, never chain
   letters, lengths or motifs;
5. **compulsory negative reporting**, including numerical-instability findings
   and at-chance published baselines;
6. **hash-pinned delivery** so audits can verify which code produced which
   number.

## 8. Reproducibility

All artifacts committed with hashes: corpus manifest
(`data/processed/fingerprint_dataset.manifest.json`), pre-registration
(`docs/fingerprint_retrieval_20261021/plan.json`), protocol amendment, fold
assignment, per-query results (`results.json`), weights and logs (`model_run/`),
ERGO protocol and outputs (`docs/six_blockers_20261021/`).

```bash
# fingerprint benchmark
PYTHONPATH=src .venv/bin/python scripts/build_fingerprint_dataset.py \
    --inventory data/releases/20261017-structural-v7-atomic-repair/native_inventory_graphs.pt
PYTHONPATH=src .venv/bin/python -m meta_acdc.models.train_fingerprint \
    --epochs 60 --batch-size 16 --fixed-coordinates
# [TCRen comparison: command and outputs to be added]
```

## References (to complete)

- Karnaukhov et al., TCRen, Nat Comput Sci 2024 (DOI 10.1038/s43588-024-00653-0)
- STAG (IEEE/ACM TCBB 2025); NetTCR-struc (Front Immunol 2025); SageTCR (Brief
  Bioinform 2025); t2pmhc (bioRxiv 2026)
- "T-cell receptor binding prediction: a machine learning revolution",
  ImmunoInformatics 2024 (validation-design critique)
- ERGO-I (Springer et al.); ZeroBind (Nat Commun 2023, DOI
  10.1038/s41467-023-43597-1)
- Cameron et al. (PMID 23926201); Linette et al. (Blood 2013); Morgan et al.
  (PMID 23377668); Raman et al. 2016 (5BRZ/MAG-IC3, PMID 26758806); Ma et al.
  Nat Commun 2025 (TCR4 SPR, DOI 10.1038/s41467-025-56004-8)
- External delivery review (docs/delivery_audit_20261018/external_review.md);
  six-blocker audit (docs/six_blockers_20261021/README.md)
