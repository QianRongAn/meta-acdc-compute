#!/usr/bin/env bash
# GPU (AutoDL 3080 Ti) production DA flow, 2026-10-12.
# 1. import the 150 downloaded AF3 results (5 folds batches, model_0 only)
# 2. retrain domain-adapted EGNN with the collate_per_graph fix (5 seeds, cuda)
# 3. score candidate structures per seed
# 4. cross-seed stability verdict (Rashomon gate)
# 5. fixed seed-filtered ensemble + KN-11 clinical safety scan
set -uo pipefail
cd /root/autodl-tmp/meta-acdc-compute || exit 1
export PATH=/root/miniconda3/bin:$PATH
PY=python
FOLDS=/root/autodl-tmp/meta-acdc-compute/folds_import

echo "== 1/6 import folds batches $(date +%F\ %T)"
for d in "$FOLDS"/folds_*; do
  [ -d "$d" ] || continue
  echo "--- import: $(basename "$d")"
  $PY src/meta_acdc/structure/native_manifest.py \
      --src "$d" \
      --out data/raw/af3_native_predictions \
      --manifest data/processed/af3_native/manifest.tsv \
      --report data/processed/af3_native/import_report.tsv \
      --candidate-out data/raw/af3_predictions \
      --list data/processed/kn5_submission_list.tsv \
    || echo "import FAILED: $d"
done
echo "natives after import: $(ls data/raw/af3_native_predictions/*_model_0.cif 2>/dev/null | wc -l)"

echo "== 2/6 retrain DA cuda 5 seeds x 150ep $(date +%F\ %T)"
$PY src/meta_acdc/models/train_domain_adapt.py \
    --base data/processed/egnn_dataset.pt \
    --af3 data/raw/af3_native_predictions \
    --af3-manifest data/processed/af3_native/manifest.tsv \
    --seeds 0 1 2 3 4 --epochs 150 --device cuda
echo "train exit: $?"

echo "== 3/6 score candidates per seed $(date +%F\ %T)"
for s in 0 1 2 3 4; do
  $PY src/meta_acdc/structure/score_predictions.py \
      --cifs data/raw/af3_predictions \
      --model "data/processed/da_seed${s}.model.pt" --mask-plddt \
      --out "data/processed/da_scores_seed${s}.tsv" || echo "score seed${s} FAILED"
done

echo "== 4/6 stability verdict $(date +%F\ %T)"
$PY src/meta_acdc/structure/da_stability.py \
    --scores data/processed/da_scores_seed0_ensemble.tsv \
              data/processed/da_scores_seed1_ensemble.tsv \
              data/processed/da_scores_seed2_ensemble.tsv \
              data/processed/da_scores_seed3_ensemble.tsv \
              data/processed/da_scores_seed4_ensemble.tsv || echo "stability FAILED"

echo "== 5/6 ensemble + clinical scan $(date +%F\ %T)"
$PY src/meta_acdc/structure/ensemble_mean.py \
    --scores data/processed/da_scores_seed0_ensemble.tsv \
              data/processed/da_scores_seed1_ensemble.tsv \
              data/processed/da_scores_seed2_ensemble.tsv \
              data/processed/da_scores_seed3_ensemble.tsv \
              data/processed/da_scores_seed4_ensemble.tsv \
    --out data/processed/prediction_scores_ensemble.tsv || echo "ensemble FAILED"
$PY src/meta_acdc/structure/clinical_scan.py \
    --out data/processed/kn11_clinical_scan.tsv || echo "clinical_scan FAILED"

echo "== 6/6 ALL DONE $(date +%F\ %T)"
