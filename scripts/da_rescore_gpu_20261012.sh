#!/usr/bin/env bash
set -euo pipefail
printf "%s\n" "Retired: this driver re-ran pre-QC scoring with legacy checkpoints and the swapped-column submission list, and would overwrite the recovered historical per-seed score tables." "See docs/delivery_audit_20261018/README.md. No structural retraining or clinical scoring is released." >&2
exit 2
