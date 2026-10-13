#!/usr/bin/env bash
# Score ColabFold (AF2-multimer) TCR-pMHC predictions with the 5-seed DA model.
# Mirrors da_rescore_gpu_20261012.sh: same seeds, same --mask-plddt, T=1.0.
# Input:  data/raw/colabfold_predictions/*_model_1.cif (from colabfold_import.py)
# Output: data/processed/da_scores_colabfold_seed{0-4}.tsv (+ _ensemble per seed)
#         data/processed/colabfold_scores_ensemble.tsv (cross-seed mean)
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python

echo "== 1/2 score per seed $(date +%F\ %T)"
for s in 0 1 2 3 4; do
  $PY src/meta_acdc/structure/score_predictions.py \
      --cifs data/raw/colabfold_predictions \
      --model "data/processed/da_seed${s}.model.pt" --mask-plddt \
      --out "data/processed/da_scores_colabfold_seed${s}.tsv" \
    || { echo "score seed${s} FAILED"; exit 1; }
done

echo "== 2/2 cross-seed ensemble $(date +%F\ %T)"
$PY src/meta_acdc/structure/ensemble_mean.py \
    --scores data/processed/da_scores_colabfold_seed0_ensemble.tsv \
              data/processed/da_scores_colabfold_seed1_ensemble.tsv \
              data/processed/da_scores_colabfold_seed2_ensemble.tsv \
              data/processed/da_scores_colabfold_seed3_ensemble.tsv \
              data/processed/da_scores_colabfold_seed4_ensemble.tsv \
    --out data/processed/colabfold_scores_ensemble.tsv \
  || { echo "ensemble FAILED"; exit 1; }

echo "== ALLDONE $(date +%F\ %T)"
