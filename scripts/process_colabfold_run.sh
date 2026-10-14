#!/usr/bin/env bash
# Generic ColabFold dump -> repo import + DA scoring + ledger update.
# Usage:
#   bash scripts/process_colabfold_run.sh <src_pdb_dir> <run_name> <status> [note_suffix]
# Example (run2):
#   bash scripts/process_colabfold_run.sh \
#       /mnt/c/Users/Administrator/WorkBuddy/2026-09-23-01-47-50/af3work/results_run2 \
#       run2 \
#       "results on disk (colabfold run2, 2026-10-14)" \
#       " [+colabfold run2 2026-10-14]"
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
SRC="$1"; RUN="$2"; STATUS="$3"; NOTE="${4:-}"
PY=.venv/bin/python

echo "== 1/4 import $SRC -> data/raw/colabfold_predictions $(date +%F\ %T)"
CIFDIR=data/raw/colabfold_predictions
ls "$CIFDIR"/*_model_1.cif 2>/dev/null | sort > /tmp/colabfold_cifs_before.txt
$PY src/meta_acdc/structure/colabfold_import.py \
    --src "$SRC" --out "$CIFDIR" \
    --report "data/processed/colabfold_${RUN}_import_report.tsv" || exit 1
ls "$CIFDIR"/*_model_1.cif 2>/dev/null | sort > /tmp/colabfold_cifs_after.txt

JOBS="/tmp/colabfold_${RUN}_jobs.txt"
comm -13 /tmp/colabfold_cifs_before.txt /tmp/colabfold_cifs_after.txt \
    | xargs -n1 basename | sed 's/_model_1.cif//' | sort > "$JOBS"
echo "new jobs this run: $(wc -l < "$JOBS")"

echo "== 2/4 ledger $(date +%F\ %T)"
ARGS=(--jobs "$JOBS" --status "$STATUS")
[ -n "$NOTE" ] && ARGS+=(--append-note "$NOTE")
$PY scripts/update_pending_ledger.py "${ARGS[@]}" || exit 1

echo "== 3/4 score 5 seeds (rescans all CIFs in colabfold_predictions) $(date +%F\ %T)"
bash scripts/score_colabfold.sh || exit 1

echo "== 4/4 ALLDONE $(date +%F\ %T)"
