#!/usr/bin/env bash
# One-shot domain-adaptation driver for the AF3-positive route.
#
# Waits for nothing: run this after the user drops the native AF3 results
# (a folds_* folder downloaded from alphafoldserver.com). It:
#   1. imports the native AF3 results -> data/raw/af3_native_predictions
#   2. trains the domain-adapted EGNN (3 seeds, pLDDT masked, full decoy protocol)
#   3. scores the candidate AF3 structures with each DA seed (--mask-plddt)
#   4. reports the cross-instance Spearman verdict (the Rashomon gate)
#
# Usage:
#   scripts/run_da_af3.sh "/mnt/c/Users/Administrator/Downloads/folds_YYYY_MM_DD_HH_MM"
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1

SRC="${1:?usage: run_da_af3.sh <folds_dir_with_native_AF3_results>}"
PY=.venv/bin/python
NATIVE_DIR=data/raw/af3_native_predictions
SCORE_CIFS="${SCORE_CIFS:-data/raw/af3_predictions}"

echo "== 1/4 build native AF3 positives (chain-fingerprint matched)"
$PY src/meta_acdc/structure/native_manifest.py \
    --src "$SRC" \
    --out "$NATIVE_DIR" \
    --manifest data/processed/af3_native/manifest.tsv \
    --report data/processed/af3_native/import_report.tsv

echo "== 2/4 train domain-adapted EGNN (3 seeds)"
$PY src/meta_acdc/models/train_domain_adapt.py \
    --af3 "$NATIVE_DIR" \
    --af3-manifest data/processed/af3_native/manifest.tsv \
    --epochs "${EPOCHS:-150}"

echo "== 3/4 score candidate AF3 structures with each DA seed"
for s in 0 1 2; do
  $PY src/meta_acdc/structure/score_predictions.py \
      --cifs "$SCORE_CIFS" \
      --model "data/processed/da_seed${s}.model.pt" --mask-plddt \
      --out "data/processed/da_scores_seed${s}.tsv"
done

echo "== 4/6 cross-instance stability verdict"
$PY src/meta_acdc/structure/da_stability.py \
    --scores data/processed/da_scores_seed0_ensemble.tsv \
              data/processed/da_scores_seed1_ensemble.tsv \
              data/processed/da_scores_seed2_ensemble.tsv

echo "== 5/6 fixed 3-seed ensemble + clinical safety scan (KN-11)"
$PY src/meta_acdc/structure/ensemble_mean.py \
    --scores data/processed/da_scores_seed0_ensemble.tsv \
              data/processed/da_scores_seed1_ensemble.tsv \
              data/processed/da_scores_seed2_ensemble.tsv \
    --out data/processed/prediction_scores_ensemble.tsv
$PY src/meta_acdc/structure/clinical_scan.py \
    --out data/processed/kn11_clinical_scan.tsv

echo "== 6/6 refresh figures"
$PY scripts/make_figures.py >/dev/null 2>&1 || echo "  (figures skipped)"

echo "== done (verdict above; STABLE iff mean pairwise Spearman r >= 0.5)"
