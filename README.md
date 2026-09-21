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
- Benchmark records: `docs/benchmarks.md`
- Paper outlines: `docs/paper1_outline.md` (Meta-TCR-GNN), `docs/paper2_outline.md` (AL safety assessment)

## Current status (2026-09-22)

**Scope**: dry-lab only; wet-lab route preserved as documentation.

**Key results so far** (details in `docs/benchmarks.md`):
- Sequence-model ceiling quantified: ProtoNet/MAML++ cluster at AUROC 0.62-0.64
  on held-out epitopes (matches TITAN/ERGO literature).
- **Side-chain contact breakthrough (v9.1)**: C-alpha-level graphs cannot
  distinguish chemistry-swapped peptides at fixed geometry (AUROC 0.50); adding
  side-chain vdW contact histograms raises graft-decoy discrimination to
  **0.933** (ablation-confirmed: 0.501 without contacts). Geometry task 0.980.
- **AF3 cross-reactivity ranking experiment (first real data)**: for the A6 TCR,
  5/8 VDJdb-validated cross-reactive peptides scored 0.96-0.99 on AF3-predicted
  structures; B7 TCR + Tax variant scored 0.999. Native-peptide controls and
  more TCR series in progress (submissions handled by the user via the AF3 web
  UI; `data/processed/af3_starter/` and `af3_json/` contain the jobs).
- Active learning simulator: EIG acquisition gives 1.8-2x positive-recall gain
  over random sampling at 6% sampling fraction; EIG > variance > entropy.

**Pipeline**: VDJdb/IEDB ingestion → interface graphs (10A interface, 8A edges,
atom-contact features) → EGNN training → AF3/TCRmodel2 structure scoring
(`src/meta_acdc/structure/score_predictions.py`) → per-TCR ranking
(`rank_analysis.py`).

## Environment

Development target: local machine (GTX 1050 Ti 4GB, 8GB RAM), Python 3.14 venv
at `.venv/` (torch 2.9.1+cu126). See `scripts/setup_env.sh`. AlphaFold 3
structure prediction runs through the official web UI (user-submitted);
TCRmodel2 is the license-safe alternative.
