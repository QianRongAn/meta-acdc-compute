#!/usr/bin/env bash
# Domain-adaptation orchestration: poll TCRmodel2 queue -> when all resolved,
# trigger DA retrain -> score AF3 -> stability verdict.
#
# Idempotent: safe to run every cron tick. Only advances stages when their
# inputs are ready. Designed to be called by crontab every 3h.
set -u
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python
STAMP=data/processed/.da_orchestrated
LOCK=/tmp/da_orchestrate.lock

log() { echo "[$(date '+%F %T')] $*"; }

# -- single-instance guard --
if [ -e "$LOCK" ]; then
  pid=$(cat "$LOCK" 2>/dev/null)
  if kill -0 "$pid" 2>/dev/null; then
    log "already running (pid $pid), exit"
    exit 0
  fi
fi
echo $$ > "$LOCK"
trap 'rm -f "$LOCK"' EXIT

# -- stage 0: poll TCRmodel2 (download any finished jobs) --
log "polling TCRmodel2 queue"
$PY src/meta_acdc/structure/tcrmodel_poll.py \
    --state data/processed/tcrmodel_jobs.tsv \
    --out data/raw/tcrmodel \
    --rounds 1 --interval 1 --workers 8

n_done=$(awk -F'\t' 'NR>1 && $3=="1"{c++} END{print c+0}' data/processed/tcrmodel_jobs.tsv)
n_pend=$(awk -F'\t' 'NR>1 && $3!="1" && $3!="failed"{c++} END{print c+0}' data/processed/tcrmodel_jobs.tsv)
n_fail=$(awk -F'\t' 'NR>1 && $3=="failed"{c++} END{print c+0}' data/processed/tcrmodel_jobs.tsv)
log "downloaded=$n_done pending=$n_pend failed=$n_fail"

if [ "$n_pend" -gt 0 ]; then
  log "queue still has $n_pend pending; nothing else to do this tick"
  exit 0
fi

# -- stage 1: domain-adaptation retrain (once all resolved) --
if [ -e "$STAMP" ]; then
  log "DA already orchestrated ($STAMP exists); skipping retrain"
  exit 0
fi
log "all TCRmodel2 jobs resolved -> starting domain-adaptation retrain"
$PY src/meta_acdc/models/train_domain_adapt.py \
    --base data/processed/egnn_dataset.pt \
    --tcrmodel data/raw/tcrmodel \
    --epochs 150

# -- stage 2: score AF3 predictions with each DA seed (pLDDT masked) --
# (AF3 CIFs are expected in data/raw/af3_predictions/; skip if absent)
if [ -d data/raw/af3_predictions ]; then
  for s in 0 1 2; do
    $PY src/meta_acdc/structure/score_predictions.py \
        --model data/processed/da_seed${s}.model.pt --mask-plddt \
        --out data/processed/da_scores_seed${s}_ensemble.tsv
  done
  # -- stage 3: cross-instance stability verdict --
  $PY src/meta_acdc/structure/da_stability.py \
      --scores data/processed/da_scores_seed0_ensemble.tsv \
                data/processed/da_scores_seed1_ensemble.tsv \
                data/processed/da_scores_seed2_ensemble.tsv
fi

touch "$STAMP"
log "DA orchestration complete"
