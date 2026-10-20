#!/usr/bin/env bash
# One-shot DA retrain after the collate fix (train/inference z-score alignment).
#
# The pre-fix models (da_seed0-4, trained 2026-10-11) scored their own
# training positives at ~0 under per-graph normalization — batch-level
# z-score overfit. train_domain_adapt.py now normalizes per graph in
# training (collate_per_graph). This script re-runs the full post-train
# flow (5 seeds, 150 epochs) on the LOCAL GPU with existing data:
# no new AF3 download needed.
#
# Run:  scripts/da_retrain_collatefix.sh   (background via nohup)
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python

echo "== 1/5 retrain DA, collate fix (5 seeds x 150ep) $(date +%F\ %T)"
$PY src/meta_acdc/models/train_domain_adapt.py \
    --base data/processed/egnn_dataset.pt \
    --af3 data/raw/af3_native_predictions \
    --af3-manifest data/processed/af3_native/manifest.tsv \
    --seeds 0 1 2 3 4 --epochs 150
echo "train exit: $?"

echo "== 2/5 score candidate AF3 structures per seed $(date +%F\ %T)"
for s in 0 1 2 3 4; do
  $PY src/meta_acdc/structure/score_predictions.py \
      --cifs data/raw/af3_predictions \
      --model "data/processed/da_seed${s}.model.pt" --mask-plddt \
      --out "data/processed/da_scores_seed${s}.tsv" || echo "score seed${s} FAILED"
done

echo "== 3/5 cross-instance stability verdict $(date +%F\ %T)"
$PY src/meta_acdc/structure/da_stability.py \
    --scores data/processed/da_scores_seed0_ensemble.tsv \
              data/processed/da_scores_seed1_ensemble.tsv \
              data/processed/da_scores_seed2_ensemble.tsv \
              data/processed/da_scores_seed3_ensemble.tsv \
              data/processed/da_scores_seed4_ensemble.tsv || echo "stability FAILED"

echo "== 4/5 fixed seed-filtered ensemble + clinical safety scan (KN-11) $(date +%F\ %T)"
$PY src/meta_acdc/structure/ensemble_mean.py \
    --scores data/processed/da_scores_seed0_ensemble.tsv \
              data/processed/da_scores_seed1_ensemble.tsv \
              data/processed/da_scores_seed2_ensemble.tsv \
              data/processed/da_scores_seed3_ensemble.tsv \
              data/processed/da_scores_seed4_ensemble.tsv \
    --out data/processed/prediction_scores_ensemble.tsv || echo "ensemble FAILED"
$PY src/meta_acdc/structure/clinical_scan.py \
    --out data/processed/kn11_clinical_scan.tsv || echo "clinical_scan FAILED"

echo "== 5/5 DONE $(date +%F\ %T)"
