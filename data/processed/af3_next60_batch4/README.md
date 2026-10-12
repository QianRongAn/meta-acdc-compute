# af3_next60_batch4 — 收尾批（47 个，队列清空）

> **这是 347 个待传队列的最后一批。** 传完本批，`pending_upload_manifest.tsv`
> 里 **所有内容唯一的 job 全部有了归宿**（batch1 60 已提交归档、batch2 49、
> batch3 60、本批 47；其余为已提交/内容重复/作废旧 ID）。
> 判重口径与前三批相同：`already_uploaded.tsv`（393 条，内容级）+ FASTA
> 内容哈希 + 磁盘已有结果。本批 47 个与任何已提交/已准备的任务 **md5 均不重复**。

## 构成（47 个）

| 类别 | 数量 | 说明 |
|---|---|---|
| 1E6 TCR 变体肽矩阵 | 47 | 5c08 7 + 5c09 8 + 5c0a 8 + 5c0b 8 + 5c0c 8 + 5hyj 8（3uts 8 + 3utt 8 + 5c07 7 + 5c08 1 已在 batch3 或已回收结果） |

本批是 **1E6 TCR（NY-ESO-1 157–165 SLLMWITQV 特异的亲和力成熟 TCR）**
× HLA-A*02:01 的 ALWGPDPAAA 家族变体矩阵收尾：9 条变体肽
（ALWGPDPAAA / AQWGPDPAAA / WGPDPAAA / MVWGPDPLYV / RQFGPDFPTI /
RQFGPDWIVA / YLGGPDFPTI / YQFGPDFPIA / RQWGPDPAAV）× 6 个受体构象。
每个 PDB 给出全部 8 槽位中的剩余部分，按 PDB 分组提交（同组共用 TCR/MHC 链）。

- **为什么值得传**：1E6 是「亲和力成熟增强特异 TCR」的经典体系，
  其变体肽正是「形状相近、化学不同」的交叉反应探针——与 MAGE-A3/titin
  致死案例（KN-11）同一机制。全家族 81 个 job（9 变体肽 × 9 受体结构）
  **至此全部有归宿**：18 个已有结果（3uts 整组 + 每结构 1 个探针）、
  16 个在 batch3、47 个在本批，9×9 矩阵完整闭合，可直接支撑 KN-11 的
  组内/跨结构外推统计（当前外推面板 154 个候选）。

## 说明

- 生成时点：2026-10-12，用 `scripts/build_af3_next_batch.py` 全量判重生成
  （脚本会重算「台账 ∪ 磁盘结果 ∪ 已入批 ∪ 内容哈希」四道闸）。
- 顺序见 `SUBMIT_NOW.tsv`（按 PDB 分组）；每个 `.fasta` 整体粘贴进
  alphafoldserver.com 即可。
- 提交后把 `folds_*` 下载放进 Downloads，自动导入回收；本文件夹里有结果的
  会移到 `af3_next60_batch4_submitted/`（机制同 batch 2/3）。
