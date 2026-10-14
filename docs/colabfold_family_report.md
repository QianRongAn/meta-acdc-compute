# ColabFold 两批(172)家族排名粗筛报告

> 2026-10-15。数据:`colabfold_scores_ensemble.tsv`(run1 87 + run2 85,5-seed DA 集成)。**QC/粗筛层**:AF2-multimer 结构与 AF3 不可互换用于DA 打分(跨方法 r=0.184, n=28;benchmarks 2026-10-14),下表仅供优先级参考,不作为风险判读;高风险判断须回到 AF3 结构。原生肽基本不在 ColabFold 批次(选的都是变体候选),故本表报各家族**头两名候选与分布**,不报原生排名。完整逐条排名:`data/processed/colabfold_rank_analysis.log`(65 个家族)。

| 家族 | 变体数 | 头名(score) | 次名(score) | 分差(max−min) | 同家族 AF3 候选数 |
|---|---|---|---|---|---|
| 5c07 | 2 | YLGGPDFPTI (0.139) | WGPDPAAA (0.056) | 0.08 | 6 |
| 5c08 | 8 | ALWGPDPAAA (0.907) | MVWGPDPLYV (0.903) | 0.85 | 0 |
| 5c09 | 8 | RQFGPDWIVA (0.987) | ALWGPDPAAA (0.701) | 0.98 | 0 |
| 5c0a | 8 | ALWGPDPAAA (0.882) | YLGGPDFPTI (0.757) | 0.85 | 0 |
| 5c0b | 8 | MVWGPDPLYV (0.861) | RQWGPDPAAV (0.804) | 0.86 | 0 |
| 5c0c | 8 | MVWGPDPLYV (0.970) | RQWGPDPAAV (0.812) | 0.97 | 0 |
| 5eu6 | 1 | YLEPGPVTA (0.347) | — | 0.00 | 5 |
| 5euo | 3 | GLCTLVAML (0.107) | MLAVISCAV (0.013) | 0.11 | 1 |
| 5hhm | 3 | GILGFVFTL (0.004) | GILEFVFTL (0.002) | 0.00 | 1 |
| 5hho | 5 | MLAVISCAV (0.562) | GILGLVFTL (0.062) | 0.56 | 1 |
| 5hyj | 8 | MVWGPDPLYV (0.995) | RQWGPDPAAV (0.989) | 0.99 | 0 |
| 5jzi | 1 | KLVALVINAV (0.049) | — | 0.00 | 1 |
| 5m00 | 5 | KAPYNFATM (0.282) | KAPANFATM (0.222) | 0.27 | 1 |
| 5m01 | 5 | KAVANFATM (0.990) | KAPYDYAPI (0.699) | 0.99 | 1 |
| 5m02 | 5 | KAPANFATM (0.850) | KAVYNFATM (0.779) | 0.47 | 1 |
| 5men | 2 | PKYVKQNTLKLAT (0.238) | NEGVKAAW (0.002) | 0.24 | 1 |
| 5nme | 2 | SLYNTIATL (0.946) | SLFNTIAVL (0.602) | 0.34 | 2 |
| 5nmf | 2 | SLFNTIAVL (0.943) | SLYNTVATL (0.185) | 0.76 | 1 |
| 5nmg | 2 | SLYNTVATL (0.985) | SLYNTIATL (0.960) | 0.03 | 1 |
| 5nqk | 1 | RLARLALVL (0.520) | — | 0.00 | 0 |
| 5til | 5 | KAVYNFATM (0.647) | KAPANFATM (0.625) | 0.63 | 1 |
| 5tje | 5 | KAPANFATM (0.829) | KAPYNFATM (0.802) | 0.82 | 1 |
| 5wkf | 1 | GTSGSPIINR (0.853) | — | 0.00 | 0 |
| 5wkh | 1 | GTSGSPIVNR (0.579) | — | 0.00 | 0 |
| 5xot | 1 | VPLTEDAEL (0.001) | — | 0.00 | 0 |
| 5yxn | 1 | KLVALVINAV (0.485) | — | 0.00 | 0 |
| 5yxu | 1 | KLVALVINAV (0.044) | — | 0.00 | 0 |
| 6am5 | 3 | MMWDRGLGMM (0.982) | AAGIGILTV (0.668) | 0.94 | 1 |
| 6amu | 3 | AAGIGILTV (0.482) | ELAGIGILTV (0.007) | 0.48 | 1 |
| 6bj3 | 1 | VPLTEDAEL (0.599) | — | 0.00 | 0 |
| 6bj8 | 1 | IPLTEEAEL (0.004) | — | 0.00 | 0 |
| 6d78 | 3 | SMLGIGIVPV (0.275) | ELAGIGILTV (0.161) | 0.27 | 1 |
| 6dkp | 3 | SMLGIGIVPV (0.709) | MMWDRGLGMM (0.634) | 0.45 | 0 |
| 6eqb | 4 | ELAGIGILTV (0.168) | EAAGIGILTV (0.005) | 0.17 | 1 |
| 6g9q | 5 | KAVANFATM (0.876) | KAVYNFATM (0.691) | 0.87 | 0 |
| 6p64 | 1 | GILGFVFTL (0.000) | — | 0.00 | 1 |
| 6q3s | 1 | SLLMWITQC (0.042) | — | 0.00 | 1 |
| 6tmo | 4 | ELAGIGILTV (0.822) | ELAAIGILTV (0.557) | 0.82 | 1 |
| 6uk2 | 1 | GILGFVFTL (0.000) | — | 0.00 | 1 |
| 6vm7 | 1 | IMDQVPFSV (0.001) | — | 0.00 | 0 |
| 6vm8 | 1 | ITDQVPFSV (0.007) | — | 0.00 | 0 |
| 6vm9 | 1 | ITDQVPFSV (0.015) | — | 0.00 | 0 |
| 6vma | 2 | ILDQVPFSV (0.097) | IMDQVPFSV (0.091) | 0.01 | 0 |
| 6zkw | 2 | RLPAKAPLLGCG (0.163) | RLPAKAPL (0.004) | 0.16 | 0 |
| 6zkx | 2 | RLPAKAPLL (0.062) | RLPAKAPL (0.046) | 0.02 | 0 |
| 6zky | 2 | RLPAKAPLL (0.720) | RLPAKAPLLGCG (0.186) | 0.53 | 0 |
| 6zkz | 2 | RLPAKAPLL (0.518) | RLPAKAPLLGCG (0.013) | 0.51 | 0 |
| 7n1f | 1 | GILGFVFTL (0.135) | — | 0.00 | 1 |
| 7n2n | 1 | LRVMMLAPF (0.062) | — | 0.00 | 0 |
| 7n2o | 1 | TRLALIAPK (0.473) | — | 0.00 | 0 |
| 7n2p | 2 | TRLALIAPK (0.279) | LRVMMLAPF (0.051) | 0.23 | 0 |
| 7n2q | 2 | TRLALIAPK (0.771) | GQVMVVAPR (0.519) | 0.25 | 0 |
| 7n2r | 2 | LRVMMLAPF (0.100) | GQVMVVAPR (0.000) | 0.10 | 0 |
| 7nme | 1 | LWMRLLPLL (0.000) | — | 0.00 | 0 |
| 7nmf | 1 | LWMRLLPLL (0.531) | — | 0.00 | 0 |
| 7nmg | 1 | QLPRLFPLL (0.152) | — | 0.00 | 0 |
| 7ow5 | 1 | VVVGADGVGK (0.732) | — | 0.00 | 0 |
| 7ow6 | 1 | VVVGAGGVGK (0.532) | — | 0.00 | 0 |
| 7q99 | 4 | LLLGIGILVL (0.122) | ELAGIGILTV (0.103) | 0.12 | 1 |
| 7q9a | 4 | ELAGIGILTV (0.826) | NLSALGIFST (0.497) | 0.63 | 1 |
| 7rtr | 1 | GILGFVFTL (0.408) | — | 0.00 | 1 |
| 8f5a | 1 | NLVPMVATV (0.005) | — | 0.00 | 1 |
| 8gon | 1 | RLQSLQTYV (0.000) | — | 0.00 | 0 |
| 8gvb | 1 | RFPLTFGW (0.040) | — | 0.00 | 0 |
| 8gvg | 1 | RYPLTFGW (0.097) | — | 0.00 | 0 |