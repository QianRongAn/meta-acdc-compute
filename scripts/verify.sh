#!/usr/bin/env bash
# One-command health check for the Meta-ACDC compute pipeline (KN-16).
#
# Runs the regression suite, package smoke imports, key data-integrity
# checks, and (if present) the domain-adaptation stability verdict + clinical
# scan. Safe to run on a fresh checkout: data-dependent checks skip.
#
# Usage:  scripts/verify.sh
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
PY=${PY:-.venv/bin/python}
fail=0

step() { echo; echo "== $* =="; }
ok()   { echo "  OK: $*"; }
bad()  { echo "  FAIL: $*"; fail=1; }

step "package smoke imports"
if $PY -c "import meta_acdc, meta_acdc.models.egnn, meta_acdc.models.train_domain_adapt, \
meta_acdc.structure.graph, meta_acdc.structure.import_af3, \
meta_acdc.structure.native_manifest, meta_acdc.structure.clinical_scan, \
meta_acdc.dashboard.server; print('ok')" >/dev/null 2>&1; then
  ok "all modules import"
else
  bad "import failure"
fi

step "regression tests"
if $PY -m unittest discover -s tests 2>&1 | tail -3; then
  ok "unittest suite"
else
  bad "unittest suite"
fi

step "native manifest == CIF set"
MAN=data/processed/af3_native/manifest.tsv
NAT=data/raw/af3_native_predictions
if [ -f "$MAN" ] && compgen -G "$NAT/*_model_0.cif" >/dev/null; then
  if diff -q <(cut -f1 "$MAN" | tail -n +2 | sort) \
             <(ls "$NAT"/*_model_0.cif | sed 's#.*/##;s/_model_0.cif//' | sort) \
             >/dev/null; then
    n=$(($(wc -l < "$MAN") - 1))
    ok "manifest matches $n native CIFs"
  else
    bad "manifest != CIF set (run native_manifest.py to rebuild)"
  fi
else
  echo "  SKIP: native manifest/CIFs absent"
fi

step "domain-adaptation stability verdict"
if compgen -G "data/processed/da_scores_seed0_ensemble.tsv" >/dev/null; then
  $PY src/meta_acdc/structure/da_stability.py \
      --scores data/processed/da_scores_seed0_ensemble.tsv \
                data/processed/da_scores_seed1_ensemble.tsv \
                data/processed/da_scores_seed2_ensemble.tsv \
      2>/dev/null | grep -E "mean pairwise|verdict"
else
  echo "  SKIP: no DA scores (run scripts/run_da_af3.sh <folds_dir>)"
fi

step "clinical safety scan (KN-11)"
if [ -f data/processed/prediction_scores_ensemble.tsv ] \
   && [ -f data/processed/clinical_gold_standard.tsv ]; then
  $PY src/meta_acdc/structure/clinical_scan.py \
      --out /tmp/kn11_verify.tsv 2>/dev/null | grep -E "FATAL|wrote" || true
else
  echo "  SKIP: scores/clinical table absent"
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "VERIFY: PASS"
else
  echo "VERIFY: FAIL ($fail check(s))"
fi
exit "$fail"