#!/usr/bin/env bash
# 234-positive (172 AF3 + 62 crystal) 3-seed DA: stability gate + ladder
# comparison vs the 172-native 5-seed baseline (188 panel, mean r=0.644).
#
# Runs AFTER the 3 checkpoints exist in $CKPT (local retrain or pulled
# from the AutoDL instance). CPU/GTX-only, candidates stay local.
#
# Usage:
#   scripts/run_da234_gate.sh [checkpoint_dir]
# Default checkpoint dir: data/processed/da_v2_crystals
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1
PY=.venv/bin/python
CKPT="${1:-data/processed/da_v2_crystals}"
OUT="data/processed/da234_gate"
mkdir -p "$OUT"

echo "== gate start $(date +%F\ %T) checkpoints: $CKPT"

# 1/4 score the 188-candidate panel with the 3 new seeds
for s in 0 1 2; do
  [ -f "$CKPT/da_seed${s}.model.pt" ] || { echo "missing $CKPT/da_seed${s}.model.pt"; exit 1; }
  if [ -f "$OUT/da234_scores_seed${s}.tsv" ]; then
    echo "seed $s already scored, skip"
    continue
  fi
  echo "== 1/4 score seed $s $(date +%T)"
  $PY src/meta_acdc/structure/score_predictions.py \
      --cifs data/raw/af3_predictions \
      --model "$CKPT/da_seed${s}.model.pt" --mask-plddt \
      --out "$OUT/da234_scores_seed${s}.tsv" --device cuda
done

# 2/4 stability verdict across the 3 new seeds
echo "== 2/4 stability verdict $(date +%T)"
$PY src/meta_acdc/structure/da_stability.py \
    --scores "$OUT/da234_scores_seed0_ensemble.tsv" \
              "$OUT/da234_scores_seed1_ensemble.tsv" \
              "$OUT/da234_scores_seed2_ensemble.tsv" | tee "$OUT/stability.txt"

# 3/4 ensemble of the 3 new seeds + comparison vs 172-native baseline
echo "== 3/4 ensemble + ladder comparison $(date +%T)"
$PY src/meta_acdc/structure/ensemble_mean.py \
    --scores "$OUT/da234_scores_seed0_ensemble.tsv" \
              "$OUT/da234_scores_seed1_ensemble.tsv" \
              "$OUT/da234_scores_seed2_ensemble.tsv" \
    --out "$OUT/da234_ensemble.tsv"

$PY - "$OUT/da234_ensemble.tsv" <<'PY' | tee "$OUT/comparison.txt"
"""Ladder comparison: 3-seed 234-native ensemble vs 5-seed 172-native
ensemble (prediction_scores_ensemble.tsv) on the common 188-panel jobs."""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

new_path = Path(sys.argv[1])
base_path = Path("data/processed/prediction_scores_ensemble.tsv")


def load(p: Path) -> dict[str, float]:
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("score"):
                out[r["job_id"]] = float(r["score"])
    return out


new, base = load(new_path), load(base_path)
common = sorted(set(new) & set(base))
a = np.array([base[j] for j in common])
b = np.array([new[j] for j in common])
r = spearmanr(a, b)[0]
print(f"common jobs: {len(common)} (new {len(new)}, base {len(base)})")
print(f"234-mix 3-seed ensemble vs 172-native 5-seed ensemble: "
      f"Spearman r={r:.3f}")
print(f"234-mix ensemble mean={b.mean():.3f} std={b.std():.3f} | "
      f"172-native ensemble mean={a.mean():.3f} std={a.std():.3f}")
PY

# 4/4 clinical scan with the new ensemble
echo "== 4/4 clinical scan $(date +%T)"
$PY src/meta_acdc/structure/clinical_scan.py \
    --scores "$OUT/da234_ensemble.tsv" \
    --out "$OUT/da234_clinical_scan.tsv" || echo "clinical_scan FAILED (non-fatal)"

echo "== gate done $(date +%F\ %T)"
echo "outputs: $OUT/{stability.txt,comparison.txt,da234_ensemble.tsv,da234_clinical_scan.tsv}"
