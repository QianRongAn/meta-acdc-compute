#!/usr/bin/env bash
# Domain-adaptation orchestration (AF3-native route).
#
# Polls the AF3 web-downloads folder for *new* folds_* directories; when one
# appears it runs the full native DA flow via scripts/run_da_af3.sh
# (chain-fingerprint import -> 3-seed train -> score candidates -> verdict).
#
# Idempotent: processed source dirs are recorded in a stamp file, so a cron
# tick never re-runs the same download.
#
# NOTE: the legacy TCRmodel2 queue route is retired (server stalled >100
# jobs); the positive source is now AF3-predicted native structures that the
# user submits via the AlphaFold web UI (terms-compliant, manual).
#
# Usage: call every few hours from cron, or run manually with a folds dir:
#   scripts/orchestrate_da.sh [folds_dir ...]
set -u
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python
STAMP=data/processed/.da_af3_processed
LOCK=/tmp/da_orchestrate.lock
DOWNLOADS="${AF3_DOWNLOADS:-/mnt/c/Users/Administrator/Downloads}"

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
touch "$STAMP"

# -- collect candidate folds dirs (explicit args or scan downloads) --
candidates=("$@")
if [ "${#candidates[@]}" -eq 0 ]; then
  while IFS= read -r d; do candidates+=("$d"); done \
    < <(ls -d "$DOWNLOADS"/folds_* 2>/dev/null | grep -v '\.zip$')
fi

ran=0
for d in "${candidates[@]}"; do
  [ -d "$d" ] || continue
  if grep -qxF "$d" "$STAMP"; then
    continue
  fi
  log "new AF3 download detected: $d -> running native DA flow"
  if scripts/run_da_af3.sh "$d"; then
    echo "$d" >> "$STAMP"
    log "completed DA flow for $d"
    ran=$((ran + 1))
  else
    log "DA flow failed for $d (will retry next tick)"
  fi
done

if [ "$ran" -eq 0 ]; then
  log "no new AF3 downloads; nothing to do this tick"
fi