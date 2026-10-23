# Cached fingerprint artifacts: corrected retrospective analysis

Historical v2 did not meet the preregistered AUROC >= 0.9 target. The six-change v1/v2 point difference is unverified; determinism replay failed. Training also used incorrect batched edge offsets. These scores are not validated gains.

Family-cluster intervals average seed pairs within query, then queries within family. They are conditional on fixed cached models/gallery, not model-population uncertainty. Procrustes is fitted only on training graphs and evaluated only on test queries. Orthogonal alignment cannot improve within-instance cosine AUROC.

| Fold | Families | Ranking rho | Conditional 95% CI | Held-out aligned cosine |
|---|---:|---:|---|---:|
| 0 | 22 | 0.9550 | ['0.9495', '0.9603'] | 1.0000 |
| 1 | 22 | 0.9752 | ['0.9712', '0.9790'] | 0.9791 |
| 2 | 22 | 0.9361 | ['0.9256', '0.9466'] | 0.9947 |
| 3 | 22 | 0.9271 | ['0.9137', '0.9398'] | 0.9998 |
| 4 | 22 | 0.9452 | ['0.9323', '0.9574'] | 0.9985 |

Historical seed-20261021 mean AUROC: 0.7944
This is a retained historical observation, not a restored result or new independent test.

| Descriptive fold mean | Value |
|---|---:|
| auroc_band_spread | 0.0344 |
| auroc_band_min | 0.7812 |
| auroc_band_max | 0.8156 |
| rank1_band_spread | 0.0487 |
| ranking_stability_rho | 0.9477 |
| rank1_agreement | 0.9812 |
| top5_jaccard | 0.5358 |
| top10_jaccard | 0.5566 |
| pair_flip_rate | 0.0920 |
| pair_flip_rate_same_family | 0.1114 |
| embedding_agreement_same_graph | -0.0090 |
| embedding_agreement_baseline | 0.9816 |
| embedding_agreement_aligned_same_graph | 0.9944 |
| embedding_agreement_aligned_baseline | 0.9816 |
| pooled_aligned_auroc | 0.8048 |
| pooled_auroc | 0.8048 |
| pooled_rank1 | 0.1637 |
