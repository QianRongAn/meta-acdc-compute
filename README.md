# Meta-ACDC

**Deciphering the Proteome-wide TCR Cross-reactivity Landscape via Active Meta-learning and Synthetic Immune Cells**

PhD research project — applicant Nan Wang, supervisor Prof. Sai T. Reddy,
Department of Biosystems Science and Engineering (D-BSSE), ETH Zurich.

Project started: 2026-09-20 (Year 1, Month 1 of the 48-month timeline).

## Repository layout

```
docs/
  proposal/       Decoded proposal text + section-by-section extracts
  literature/     Literature review (Year 1 Q1 deliverable)
  milestones/     Timeline tracker, key decision nodes, experiment log
src/meta_acdc/
  data/           Data curation pipelines (VDJdb, IEDB, ACDC in-house)
  models/         EGNN encoder, meta-learning (MAML / Prototypical Networks)
  structure/      AlphaFold 3 / ColabFold interface, interface-graph construction
  active_learning/  EIG acquisition, epsilon-greedy, dry-wet loop orchestration
  dashboard/      TCR-Safety-Radar web dashboard
notebooks/        Exploration and benchmarking notebooks
data/             raw/ (gitignored) and processed/ datasets
```

## Three objectives

1. **Meta-TCR-GNN** — geometric deep learning (EGNN) on AlphaFold 3-predicted
   TCR-pMHC interface graphs, trained with hard-negative mining and meta-learning
   for few-shot generalization.
2. **AI-driven "Dry-Wet" closed loop** — active learning (Expected Information
   Gain, epsilon-greedy) to select the most informative 0.05% of candidates for
   ACDC validation; 3-5 iterations.
3. **Proteome-wide Cross-reactivity Risk Dashboard** — clinical-grade safety
   assessment for therapeutic TCRs (e.g., Kimmtrak analogs).

## Key resources

- Proposal text: `docs/proposal/proposal_decoded.txt`
- Timeline & key nodes: `docs/milestones/KEY-NODES.md`
- Literature review: `docs/literature/`

## Environment

Development target: local machine (GTX 1050 Ti 4GB, 8GB RAM). Prototyping and
benchmarks run locally; large training runs and AlphaFold 3 structure
prediction use cloud GPUs / the AlphaFold Server API (budgeted in the proposal).
