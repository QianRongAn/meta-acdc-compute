# af3_next60 — 下一批 AF3 网页版提交（batch 2，待传 49 个）

> **本文件夹现在只放「待提交」的文件**：初始 60 个中已有 11 个提交并回收结果，
> 那 11 个已移出到 `af3_next60_submitted/`。所以本文件夹里 49 个 `.fasta`
> **全部都是没提交过的**，按 `SUBMIT_NOW.tsv` 的顺序挨个传即可，不用担心重复。
> 每提交完一批，把新下载的 `folds_*` 放进 Downloads，我会把对应文件再移出去，
> 文件夹里永远只留待传的。

按优先级挑选、**尚未提交**（且 5 链 FASTA 内容与任何已提交任务都不同）的 job，
供手动复制到 [alphafoldserver.com](https://alphafoldserver.com) 提交。
每个 `.fasta` 内含 `>alpha / >beta / >peptide / >mhc / >b2m` 五条链，
直接整体粘贴进 AF3 网页版即可。

> **上一批（batch 1）** 已全部提交，归档在 `af3_next60_batch1/`。
> 教训：batch 1 只按 job_id 判重，结果 60 个里有 24 个的 5 链 FASTA 与更早提交的
> 任务**逐字节相同**（同一 TCR-pMHC 复合物的不同 PDB 别名，AF3 只看序列），
> 提交时被 AF3 判定为重复。**本批判定口径改为「FASTA 内容哈希」**：
> 凡与已提交任务内容相同的 job 一律排除（见下方判定说明）。

## 判定口径（内容级判重）

「已提交」= 满足以下任一：

1. `af3_pending_upload/already_uploaded.tsv` 记录在案；
2. `data/raw/af3_predictions/`、`af3_native_predictions/`、`af3_instances/{A,B1-B4}/`、
   `af3_predictions_quarantine_20261012/` 里已有结果文件；
3. 5 链 FASTA 内容（md5）与上述任一已提交任务的 FASTA 相同 —— **同一序列
   交给 AF3 只会得到同一预测，属于重复提交**；
4. `obsolete_relabeled_do_not_reupload.tsv` / `af3_relabel_map.tsv` 里的作废旧 ID。

满足条件的候选池：`af3_pending_upload/pending_upload_manifest.tsv`（347 个）中
未提交且内容唯一的 **177 个**，本批取其中 60 个。

## 构成（60 个）

| 类别 | 数量 | 说明 |
|---|---|---|
| 变体矩阵补全（6 个临床表位家族） | 32 | 把上一批没拿全的变体拿全：Melan-A/MART-1 家族 18、HPV16 E7 4、HIV gag 3、Flu M1 3、CMV pp65 2、KRAS 2 |
| 临床相关单点/小组 | 28 | MAGE（MMWDRGLGMM）4、p53 5、EBV（GLCTLVAML/MLAVISCAV/LMP2 LPEPLPQ*）5、PKYVKQNTLKLAT 2、NEGVKAAW 1、GTSGSP* 2、VPLTEDAEL/IPLTEEAEL 3、RLARLALVL 1、RLQSLQTYV 1、R*PTFGW 3、LLLGIGILVL 1 |

筛选规则：

1. **补全变体矩阵**：同一临床表位的剩余变体全部纳入（每条 pdb 结构至少覆盖
   一个变体，同批内不再拆分家族）；
2. **临床重要性排序**：MAGE > p53 > EBV > 其余（与 `clinical_gold_standard.tsv`
   的致命 off-target 案例对齐：MAGE-A3 A3A 是 KN-11 的核心案例）；
3. 候选之间保证 **md5 内容唯一**，与已提交池也保证内容唯一。

## 留给下一批（batch 3，共 117 个内容唯一候选）

- KAP 变体矩阵（KAPANFATM 家族，约 25 个）
- WGPDPAAA 变体矩阵（3uts/3utt/5c07-5c0c/5hyj，约 72 个）
- RLPAKAPL 家族（6zk*，8）、LRVMMLAPF/TRLALIAPK/GQVMVVAPR（7n2*，8）、
  LWMRLLPLL/QLPRLFPLL（3）、NLSALGIFST（1）等

## 提交后

把 `folds_*` 下载结果放回 `/mnt/c/Users/Administrator/Downloads/`，
`orchestrate_da.sh` 的定时任务会自动导入（`native_manifest.py` 累积匹配 +
域适应重训）。
