#!/usr/bin/env bash
# v11 atom-level electrostatics ablation (README 6.2.1), multi-seed.
#
# Controlled comparison on ONE dataset build (egnn_dataset_v11.pt, 28/18 dims):
# features are sliced at collate time, so all configs share identical graphs,
# decoys, and PDB-grouped split. Multiple seeds per config (model init + batch
# order; the split itself is fixed at seed 42).
#
#   base  26/16  v9.1 feature set (control)
#   node  28/16  + esp, burial node features
#   edge  26/18  + coul, coul_mag edge features
#   full  28/18  v11 (all)
#
# Usage: bash scripts/train_v11_ablation.sh [epochs] [seed0,seed1,...]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

EPOCHS="${1:-30}"
SEEDS="${2:-0,1,2}"
IFS=',' read -ra SEED_LIST <<< "$SEEDS"
DATA="data/processed/egnn_dataset_v11.pt"
OUTDIR="data/processed/v11_ablation"
mkdir -p "$OUTDIR"

run() {
  local name="$1" nd="$2" ed="$3" seed="$4"
  echo "=== [$name] node=$nd edge=$ed seed=$seed epochs=$EPOCHS $(date '+%F %T') ==="
  .venv/bin/python src/meta_acdc/models/train_egnn.py \
    --data "$DATA" --epochs "$EPOCHS" \
    --node-dim "$nd" --edge-dim "$ed" --seed "$seed" \
    --out "$OUTDIR/v11_${name}_s${seed}.model.pt" \
    2>&1 | tee "$OUTDIR/v11_${name}_s${seed}.log"
}

for s in "${SEED_LIST[@]}"; do
  run base 26 16 "$s"
  run node 28 16 "$s"
  run edge 26 18 "$s"
  run full 28 18 "$s"
done

echo "=== per-seed results $(date '+%F %T') ==="
for s in "${SEED_LIST[@]}"; do
  for n in base node edge full; do
    printf "%-4s s%s: " "$n" "$s"
    grep -E "val AUROC=|val \[(graft|displaced)\]" "$OUTDIR/v11_${n}_s${s}.log" \
      | tail -3 | sed 's/.*\(val .*AUROC=[0-9.]*\).*/\1/' | tr '\n' ' '
    echo
  done
done
