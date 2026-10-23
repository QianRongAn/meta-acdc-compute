# Retrospective shortcut diagnostics

Historical v2 weights used incorrect batch edges. All numbers below diagnose artifacts; none restore a validated gain.

| Representation | Pooled AUROC | Query macro | Family equal (conditional 95% CI) | Same-family macro |
|---|---:|---:|---|---:|
| composition_cosine | 0.9779 | 0.9806 | 0.9872 [0.9826727844276748, 0.9918872891910896] | 0.5494 |
| size_scalar_cosine | 0.5000 | 0.5000 | 0.5000 [0.5, 0.5] | 0.5000 |
| size_log_distance | 0.8157 | 0.8147 | 0.8112 [0.76431116140239, 0.8485206626880396] | 0.5079 |
| radius_distance | 0.7774 | 0.7711 | 0.7636 [0.7107595248876065, 0.8070984886722243] | 0.5208 |
| domain_metadata | 0.6344 | 0.5886 | 0.5668 [0.5052957111245682, 0.6239574714230108] | 0.4417 |
| historical_v2_seed20261021 | 0.6506 | 0.7876 | 0.7850 [0.7230148595882623, 0.8357369655534362] | 0.4962 |

All positives are retained; only negatives are restricted in each stratum. Undefined queries are counted in JSON.
Scalar positive size cosine is degenerate by design; log-size distance is the predeclared sensitivity analysis.
Domain metadata is a confound diagnostic, not evidence the confidence feature varies (it does not in this pool).
Identical composition copies have rho=1 by construction where defined; this is not a trained permutation control.
No audited MHC field: same-MHC/pMHC analysis is unavailable. No fresh binding acceptance data were read.
