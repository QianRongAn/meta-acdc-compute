# TCRen 对照协议（预注册，2026-10-21，运行前）

目标：在同一核验语料上运行已发表的结构工具 TCRen（Nat Comput Sci 2024，
DOI 10.1038/s43588-024-00653-0，官方发布，不改源码、不调参数、不选权重），
作为论文 §6.2 的外部结构基线。结果无论方向都按此协议原样报告。

## 输入与子集

- 语料：`data/processed/fingerprint_dataset.pt` 的 331 张图对应结构。
  TCRen 需要结构文件：v7 资料的结构文件路径取自 `structure_records.json`
  （晶体源 PDB 与预测 CIF），候选取 `data/raw/af3_predictions/*_model_0.cif`。
- 受体子集：拥有“已核验天然（cognate）肽结构 + 至少 1 个同受体非同肽候选”
  的受体。cognate 肽来自角色头参考的 peptide 字段；候选肽 = job 后缀。
- 排除：身份未匹配的 9 份候选、图构建失败的 4 份、任何角色头参考缺失的受体。

## 协议

1. 对每个结构运行 `tcren score`（或文档指定的最小评分命令），记录其
   肽-接触能量/评分与命令、版本、哈希。
2. 每个受体单独评估：在该受体全部结构（天然 + 候选）中，按 TCRen 分数排序，
   cognate 是否排第一（rank-1 命中）；同时记录 cognate 的平均百分位排名。
3. 机会水平：逐受体按候选数计算（1/n_candidates）。
4. 与本文嵌入检索的可比口径：仅用"同受体、非同肽"判别，避免跨任务比较；
   明确声明 TCRen 的设计任务是表位排序，本文任务是受体身份检索——只报告
   两个工具在同一受体数据上的判别能力，不宣称任务等价。

## 输出

- `docs/fingerprint_retrieval_20261021/tcren_comparison/`：
  `protocol.json`、逐结构分数 `scores.tsv`、逐受体排序 `ranking.json`、
  汇总 `summary.json`（rank-1、平均百分位、机会水平、n 受体/n 结构）。
- 全部原样保留，不因结果调整子集或口径。
