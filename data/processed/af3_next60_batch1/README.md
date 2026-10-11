# af3_next60 batch 1 — 已全部提交(2026-10-12 归档)

> ⚠️ 本批 60 个已于 2026-10-12 全部提交到 alphafoldserver.com。
> 事后按「5 链 FASTA 内容哈希」复核:其中 **24 个**与更早提交的任务内容逐字节
> 相同(同一复合物的不同 PDB 别名),属于重复提交;真正的新任务 45 个左右。
> 下一批改判重口径为内容哈希,见 `../af3_next60/README.md`。

以下为归档前的原始说明。

按优先级挑选的 60 个**尚未提交**(无任何结果文件)的 5 链 FASTA,供手动
复制到 [alphafoldserver.com](https://alphafoldserver.com) 提交。
每个 `.fasta` 内含 `>alpha / >beta / >peptide / >mhc / >b2m` 五条链,
直接整体粘贴进 AF3 网页版即可。

判定口径:凡在 `data/raw/af3_predictions/`、`af3_native_predictions/`、
`af3_instances/{A,B1-B4}/` 中已有结果文件的 job 均视为已提交,本批不再包含。

## 构成

| 类别 | 数量 | 说明 |
|---|---|---|
| native | 3 | 天然复合物补全(域适应正样本,manifest 尚缺的最后 3 个) |
| kimmtrak | 5 | KN-11 临床案例 3(gp100 TCR,5eu6):靶肽 + 变体扫描 |
| candidate | 52 | 交叉反应候选(batch1 剩余 2 + batch2 临床重要表位多样化 50) |

候选筛选规则:先按表位临床重要性(GILGFVFTL / ELAGIGILTV / SLLMWITQV /
NLVPMVATV / LLFGYPVYV 等)P1>P2 排序,再保证每条 pdb 结构至少覆盖一个
变体(多样化),最后按 job_id 补足 60 个。详见 `manifest.tsv`。

## 提交后

把 `folds_*` 下载结果放回
`/mnt/c/Users/Administrator/Downloads/`,`orchestrate_da.sh` 的定时任务
会自动导入(`native_manifest.py` 累积匹配 + 域适应重训)。
