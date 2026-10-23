# Corrected fingerprint training: frozen developmental experiment

This panel was exposed previously. No hyperparameter tuning, no early stopping, no best-seed selection; no independent measured acceptance set was read.

| Representation | All-negative family-equal AUROC | Same-family family-equal AUROC |
|---|---:|---:|
| composition | 0.9872 | 0.5475 |
| size_cosine_constant | 0.5000 | 0.5000 |
| size_log_distance | 0.8112 | 0.4650 |
| radius_distance | 0.7636 | 0.4710 |
| true (3 seeds averaged per query) | 0.9379 | 0.4385 |
| shuffled (3 seeds averaged per query) | 0.7368 | 0.4280 |

Primary corrected model minus composition (same-family): -0.10894611523936786, conditional CI95 [-0.37058784234876363, 0.12198534182801046].
Positive developmental evidence criterion met: False. No claim of practical or biological improvement is allowed.
Full replay: {'weights_bitwise_equal': True, 'embeddings_bitwise_equal': True, 'scope': 'full fold0 first seed on same CPU/software; not cross-device guarantee'}
All checkpoints, training-only label mappings, embeddings, similarities, metric fragments and logs are retained; evaluation.json binds hashes.
Only 3 label permutations per fold: a descriptive control, not a calibrated permutation p-value. Same-MHC analysis unavailable.
