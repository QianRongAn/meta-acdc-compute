#!/usr/bin/env bash
set -euo pipefail
printf "%s\n" "Retired: this driver used unfrozen inputs, pre-QC results and unvalidated clinical interpretation." "See docs/delivery_audit_20261018/README.md. No structural retraining or clinical scoring is released." >&2
exit 2
