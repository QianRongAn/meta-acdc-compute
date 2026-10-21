# Receptor fingerprinting from TCR–pMHC structure: the first verifiable positive result of an audit-reconstructed modeling pipeline

**Status**: preprint draft v0.1 (2026-10-21). All numbers derive from committed,
hash-pinned artifacts (`docs/fingerprint_retrieval_20261021/`). This manuscript
claims **no binding prediction and no clinical utility**; it reports one
pre-registered representation-quality benchmark, one honest negative comparison,
and the audit-driven engineering that made both possible.

---

## Abstract

Structural models of TCR–pMHC interaction promise to capture geometric and
chemical complementarity that sequence-only methods miss, but pipelines in this
space are fragile: our own previous results were fully withdrawn after external
review identified 28 audit items and a further 39-item delivery review found
input-encoding, splitting, QC-gate and receptor-identity defects whose fixes had
not reached the code that produced the reported numbers. Here we report the
reconstruction: after byte-level, hash-pinned repair of the pipeline, we
pre-registered and ran a representation-quality benchmark that requires no
experimental binding labels and whose labels are verified metadata. On a
hash-manifested corpus of 331 interface graphs (168 atom-screened native
structures, 145 crystal and 23 predicted, plus 163 identity-verified predicted
candidates; 160 receptors, 110 families), an EGNN embedding trained with a
supervised-contrastive objective retrieves the correct receptor among all other
structures at **rank-1 accuracy 0.166 vs a per-query chance level of 0.0123
(13.4×, minimum 7.2× across five family-grouped cross-validation folds; 226
queries)**, with mean reciprocal rank 0.284 and a same-vs-different-receptor
pair AUROC of 0.767. Family classification by k-nearest neighbours reaches 0.224
top-1 accuracy vs 0.071 chance. A domain diagnostic confirms the representations
do not separate crystal from predicted structures (accuracy 0.496, chance 0.5),
i.e. the result is not a domain or confidence artifact. We also report the
equivariant coordinate-update branch diverges numerically under the contrastive
objective on this corpus — a stability finding documented in a protocol
amendment. In parallel, an official-weight ERGO-I control on 10,000 measured
MIRA functional observations scores at chance for the functional endpoint
(per-peptide macro-AUROC 0.511 vs 0.500 for our SPR-derived baseline), an
honest negative that bounds what sequence-only inputs carry on this endpoint.
Together these results establish that the reconstructed pipeline now demonstrably
encodes receptor identity in its structure representations — a necessary
prerequisite for any downstream binding or off-target application, which remains
gated on real experimental labels.

## 1. Introduction

### 1.1 Motivation

Adoptive TCR therapies carry the risk of cross-reactive off-target killing
(Cameron et al. 2013; Linette et al. 2013; Morgan et al. 2013). Predicting which
peptide–MHC complexes a given TCR will bind remains open; sequence-based tools
are the workhorse but operate without explicit interface geometry or chemistry.
Structure-based discriminators were the motivation of our pipeline (ACDC
project), but our reported results — graft/decoy discrimination (AUROC 0.933),
cross-instance ranking stability (0.644), and clinical index scans — were all
withdrawn after external review found defects in residue encoding, peptide-node
pooling, batch normalization, split grouping, decoy construction, QC gating and
clinical receptor identity, and further found that several claimed fixes had not
reached the code that produced the reported numbers.

### 1.2 This paper

We report the reconstruction and its first verifiable positive result:

1. a **pre-registered representation-quality benchmark** — receptor
   fingerprinting — whose labels are verified structural metadata, not binding
   data (Sec. 3–5);
2. the **engineering that made the benchmark possible**: hash-pinned input
   manifests, identity verification against role-annotated references, atomic
   screening, family-grouped splits, and provenance receipts (Sec. 2);
3. an **honest negative**: an official-weight ERGO-I control at chance on a
   measured functional endpoint (Sec. 6);
4. a **stability finding**: divergence of the equivariant coordinate-update
   branch under the contrastive objective, reported as a protocol amendment
   (Sec. 4.4).

No binding or clinical claim is made anywhere in this manuscript; withdrawn
historical numbers are cited only as history.

## 2. Audit-driven reconstruction

### 2.1 What the reviews found

The first external audit (28 items) covered label semantics (VDJdb confidence
used as binding labels, IEDB endpoint conflation), model-input defects (three-
letter-to-one-letter residue mapping collapsing ARG/ASN/ASP to A, peptide
pooling reading the wrong feature column, batch-wide coordinate normalization),
splitting (no complex grouping), decoy construction (stale centroids, fixed-axis
displacement, unproven non-binding), and clinical reference errors (wrong PMIDs,
receptor construct confusion). A second, independent delivery review raised 39
items (5 fatal, 18 severe), including two delivered result tables contradicting
each other, effect sizes below measurement dispersion (median per-job SD 0.2512
vs interpreted deltas 0.02–0.06), receptor identity substitution (MAG-IC3 scored
and labelled as the clinical a3a construct), QC flags bypassed in delivered
clinical rows, and several "fixed" defects still present in the code paths that
had produced the reported numbers.

### 2.2 What was rebuilt

All fixes were implemented and verified at the byte level: identity assignment
now requires explicit role annotations from reference records; decoy
construction raises on length mismatch or incomplete backbones and translates
peptides along the peptide-to-groove direction; scoring requires per-model QC,
release fingerprints and provenance receipts; family-aware splitting and
family-stratified Fisher pooling are enforced. Five-seed replay of the withdrawn
run reproduced 938/940 job×seed means to four decimals, and the 347-item pending
panel was verified byte-identical to its first-commit state (98 exact αβ pairs,
45 MHC, 11 β2m sequences). These audits establish that the inputs of the current
benchmark are what they claim to be.

### 2.3 Data integrity policy carried into the benchmark

Every input of the benchmark is hash-pinned; fold assignment was persisted
before training; results were written verbatim; success criteria were fixed
before training and are evaluated as-is below.

## 3. Corpus

### 3.1 Verified native inventory

168 native interface graphs from the atomic-repaired inventory (v7 release):
145 crystal structures and 23 AF3-predicted re-predictions, each passing the
heavy-atom completeness screen, the construction screen, and reference-sequence
identity verification. Confidence (pLDDT) features are masked everywhere so the
domain of a structure cannot leak through confidence patterns.

### 3.2 Identity-verified predicted candidates

From the 188 scored candidate CIFs, 179 chain A/B sequence pairs matched an
exact role-annotated reference (alpha,beta) pair; after peptide annotation and
graph construction, 163 entered the corpus, 16 were excluded and recorded.
Receptor identity is therefore assigned only by exact sequence match to
references, never by chain letter, length, or motif heuristics.

### 3.3 Composition

331 graphs; 160 receptors; 110 families; 55 receptors with ≥2 structures (max
11); peptide lengths 8–13 (181 nine-mers, 95 ten-mers). Receptors with multiple
structures are what make retrieval evaluation possible.

## 4. Pre-registered benchmark

The protocol (PLAN.md, plan.json, committed 2026-10-21, before any training)
specifies:

- **Model**: EGNN (6 layers, hidden 128, node/edge feature dims 26/16 — the same
  feature schema as the scoring path) + a 64-dim projection head with
  LayerNorm.
- **Objective**: supervised contrastive (SupCon) with positives = same receptor
  identity, temperature 0.1; Adam (lr 3e-4, weight decay 1e-4); 60 epochs;
  batch 16.
- **Splits**: 5-fold cross-validation grouped by receptor family (a family never
  spans train and test).
- **Task A (primary)**: receptor retrieval — each test graph queries against all
  other graphs by cosine similarity of embeddings; correct if the nearest
  neighbour shares the receptor identity. Metrics: rank-1 hit rate, MRR,
  same-vs-different-receptor pair AUROC; the per-query chance level is computed
  from the gallery composition.
- **Task B (secondary)**: family classification by k-NN (k ∈ {1,3,5}).
- **Task C (diagnostic)**: nearest-centroid domain (crystal/predicted)
  classification, to detect residual domain/confidence leakage.
- **Pre-registered success criteria**: (i) retrieval rank-1 ≥ 2× chance;
  (ii) pair AUROC ≥ 0.9; (iii) family top-1 ≥ 2× chance; (iv) same-seed
  reproducibility. Failures are reported as-is.

### 4.4 Protocol amendment (stability finding)

During the first complete run, the equivariant variant (coordinate-update branch
enabled, matching the scoring path) completed fold 0 (rank-1 0.238, AUROC 0.710)
but diverged deterministically in fold 1 around epoch 13–15 (non-finite loss).
Isolated replays showed the divergence point is sensitive to CUDA reduction
noise (same seed: one replay diverged, one completed 20 epochs), i.e. the
training dynamics sit at a bifurcation edge; the mechanism is a feedback loop
in which contrastive gradients flow back through the pooled coordinate-norm
feature into the coordinate-update branch, amplifying inter-atomic distances and
edge messages. We did not tune around this; we adopted the pipeline's own
sanctioned invariant variant (coordinate updates disabled, an existing
first-class mode) for the full benchmark and record the amendment and both
variants' outputs verbatim. The instability of the equivariant branch under
contrastive training is itself reported as a finding.

## 5. Results

### 5.1 Receptor retrieval (Task A, primary)

| Fold | rank-1 | chance | ratio | MRR | AUROC | queries |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0.2857 | 0.0105 | 27.1× | 0.3814 | 0.7580 | 42 |
| 1 | 0.1167 | 0.0162 | 7.2× | 0.2454 | 0.7518 | 60 |
| 2 | 0.2532 | 0.0194 | 13.0× | 0.3672 | 0.7589 | 79 |
| 3 | 0.1250 | 0.0101 | 12.4× | 0.2389 | 0.7076 | 24 |
| 4 | 0.0476 | 0.0055 | 8.7× | 0.1895 | 0.8569 | 21 |
| mean | **0.1656** | **0.0123** | **13.4×** | **0.2845** | **0.7666** | 226 |

Across all folds and 226 queries, the embedding retrieves the correct receptor
at 13.4× the chance level (minimum 7.2× in any fold). Criterion (i) passes.
Criterion (ii) does not: the pair AUROC is 0.767 (fold range 0.708–0.857),
below the conservative 0.9 bar — reported as a limitation, not redefined away.

### 5.2 Family classification (Task B)

k-NN family classification reaches top-1 accuracy 0.224 (chance 0.071, 3.2×);
k=3 gives 0.246, k=5 0.224. Criterion (iii) passes.

### 5.3 Domain diagnostic (Task C)

Nearest-centroid crystal-vs-predicted accuracy is 0.496 (chance 0.5): under
confidence masking, the representations do not separate structure provenance.
This is a leakage check, and it passes.

### 5.4 Pre-registration verdict

2 of 3 quantitative criteria met (13.4× retrieval; 3.2× family classification);
the AUROC criterion (≥0.9) not met (0.767); reproducibility criterion met for
the invariant variant, with the equivariant divergence recorded as an amendment.
The benchmark's headline is therefore the retrieval result, reported exactly as
pre-registered.

## 6. Honest negative control: ERGO-I on measured functional data

As part of the same reconstruction effort, the official ERGO-I model (fixed
McPAS LSTM checkpoint, git-blob-verified) was run on 10,000 measured MIRA
functional observations under a pre-registered protocol (input-order controls,
per-peptide metrics, no threshold search). Per-peptide macro-AUROC is 0.5112
versus 0.5001 for the SPR-derived baseline; macro-AP 0.1166 versus 0.1197.
Both are at chance for this endpoint. This bounds what the sequence-only input
carries on the functional endpoint — it neither supports nor undermines the
fingerprinting result, and we do not claim a method comparison across tasks.

## 7. Discussion

The reconstruction is now in a state where a pre-registered experiment can
produce a verifiable positive result: structure-derived representations encode
receptor identity well beyond chance, in a family-disjoint protocol, without
domain leakage. This matters for the original research question in two ways.
First, it demonstrates that the repaired pipeline actually extracts
structure-content features — the precondition for any structure-vs-sequence
argument. Second, it provides a representation-quality gate that future binding
models must pass before their scores can be interpreted.

Limitations. Receptor identity is metadata, not binding: success here does not
imply binding prediction. The corpus is redundant (55 of 160 receptors carry
the 226 queries; the panel over-represents a few MHC sequences), and the AUROC
of 0.767 shows same-family receptors remain confusable. The invariant variant
was adopted after the equivariant branch proved numerically unstable under the
contrastive objective — the benchmark is of the invariant architecture. No
experimental binding labels are used anywhere; downstream binding work remains
gated on real negatives and is outside this manuscript.

## 8. Reproducibility

All artifacts are committed with hashes: dataset manifest
(`data/processed/fingerprint_dataset.manifest.json`), fold assignment,
pre-registration (`docs/fingerprint_retrieval_20261021/plan.json`), protocol
amendment, per-fold and per-query results (`results.json`), weights and logs
(`model_run/`). Commands:

```bash
PYTHONPATH=src .venv/bin/python scripts/build_fingerprint_dataset.py \
    --inventory data/releases/20261017-structural-v7-atomic-repair/native_inventory_graphs.pt
PYTHONPATH=src .venv/bin/python -m meta_acdc.models.train_fingerprint \
    --epochs 60 --batch-size 16 --fixed-coordinates
```

## References (to complete)

- Cameron et al., Sci Transl Med 2013 (PMID 23926201)
- Linette et al., Blood 2013 (PMID 23999400 — note: corrected to 23926201 in
  project data; verify before citation)
- Morgan et al., Cancer Immunol Immunother 2013 (PMID 23377668)
- Raman et al., Sci Rep 2016 (PMID 26758806; 5BRZ/MAG-IC3 structure)
- Ma et al., Nat Commun 2025 (TCR4 SPR, DOI 10.1038/s41467-025-56004-8)
- ERGO-I (Springer & Besser et al.); ZeroBind (Nat Commun 2023,
  DOI 10.1038/s41467-023-43597-1); external delivery review
  (`docs/delivery_audit_20261018/external_review.md`); six-blocker audit
  (`docs/six_blockers_20261021/README.md`).
