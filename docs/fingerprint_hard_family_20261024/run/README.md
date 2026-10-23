# Frozen 2x2 diagnostic results

2x2 controlled diagnostic on historical observed-sequence-key labels; NOT biological receptor truth, NOT new independent validation

NO biological receptor-model promotion: original labels are observed-sequence hashes with confirmed false negatives.

| Arm | Fit hard AUROC | Validation hard AUROC | Validation known-same-excluded AUROC |
|---|---:|---:|---:|
| reference | 0.6515 | 0.4289 | 0.4289 |
| identity_preserved | 0.6036 | 0.4590 | 0.4590 |
| family_batches | 0.9316 | 0.4648 | 0.4648 |
| identity_and_family | 0.9235 | 0.4841 | 0.4841 |

Composition validation hard AUROC: 0.5321428571428571
Family-equal conditional intervals, all seeds retained; eligible validation families: 3.
Known-same exclusion does not establish remaining negatives as true biological differences. No canonical identity labels were imputed.
