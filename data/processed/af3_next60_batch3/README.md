# af3_next60_batch3 — 第三批提交包（60 个，预建）

> **先把 `af3_next60/`（batch 2 剩余 49 个）传完，再传这一批。**
> 判重口径与 batch 2 相同：FASTA 内容哈希 + 台账（`already_uploaded.tsv`）。
> 本批 60 个与任何已提交/已准备的任务内容均无重复。

## 构成（60 个）

| 类别 | 数量 | 说明 |
|---|---|---|
| KAP 变体矩阵 | 24 | KAPANFATM 家族（5m00/5m01/5m02/5til/5tje/6g9q × KAP*变体），batch1/2 已覆盖靶肽，本批补全矩阵（1 个已有结果自动剔除） |
| WGPDPAAA 变体矩阵 | 16 | 3utt 全 8 个 + 5c07 7 个 + 5c08 1 个（3uts 8 个已被你提交并回收结果，自动跳过） |
| RLPAKAPL 家族 | 8 | 6zkw/6zkx/6zky/6zkz × 长度变体 |
| LRVMMLAPF/TRLALIAPK/GQVMVVAPR | 8 | 7n2* 家族 |
| LWMRLLPLL/QLPRLFPLL | 3 | 7nme/7nmf/7nmg |
| NLSALGIFST | 1 | 7q9a |

## 说明

- 生成时点：2026-10-12，基于当时「已有结果 + 已入批文件夹」的全量判重。
- 提交方式与 batch 2 完全一样：每个 `.fasta` 整体粘贴进 alphafoldserver.com。
- 提交后把 `folds_*` 下载放进 Downloads，自动导入；本文件夹里有结果的会移到
  `af3_next60_batch3_submitted/`（机制同 batch 2）。
