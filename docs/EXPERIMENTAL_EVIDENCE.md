> **最新状态：** 已扩为 4 项原始研究、23 条分终点观察，并完成保留家族隔离。见 [扩展核查与训练边界](EXPERIMENTAL_EVIDENCE_EXPANSION.md)。下文为首批面板历史记录；其计数和旧版本路径不代表当前版本。

# 实验标签补充与外部保留面板（2026-10-16）

本轮完成的是公开证据获取、条件性标签整理、泄漏排查和一个小型外部核验面板。未训练模型，未生成预测性能。hy4 按用户要求暂不执行；历史 HTTP 402 记录不再作为本轮阻断理由。**这不代表 28 项问题整体闭环，也不代表获得了可训练的大规模真实结合基准。**

## 已实际取得的数据

通过 [IEDB 官方 Query API](https://query-api.iedb.org/docs/swagger/) 查询人源、配对 αβ CDR3、至少关联一个阴性结果的受体组，得到：

| 项目 | 数量及解释 |
|---|---|
| 候选受体组 | 122；数据库的组不能替代具体受体身份，可能含 CDR3 以外的工程突变 |
| 受体导出条目 | 531；使用具体 curated receptor ID 与 assay ID 连接 |
| 不同实验编号 | 1,174 |
| 实验—受体—导出版本证据行 | 1,464；不是 1,464 个独立样本 |
| 数据库阳性／条件性阴性观察 | 1,197／267；所有这些行的训练 label 仍为空，待核对原文 |
| 需要区分条件的正负冲突组 | 105；不做“阳性优先”或多数投票 |
| 通过初步字段筛查的优先原文核查行 | 20；仍不是已验收标签 |

查询刻意选择含阴性结果的受体，存在选择偏倚，不可用这些数据估计总体阳性率或直接声称代表性基准。没有把 VDJdb confidence、结构诱饵、随机配对、未筛中条目当作实验阴性。

每条证据保留实验编号、具体受体编号、受体组、双链序列（有则保留）、肽、MHC、实验终点、原始数值/单位/不等号、PMID、图表位置、实验备注、效应细胞/APC及培养条件和来源链接。完整双链缺失、接头不符、组或实验关联不符、修饰肽、MHC 不明确、截尾测量以及与结构数据重叠均单列原因。不同模态（SPR、多聚体、细胞激活、细胞毒）不混算。

## 首批经原文核实的面板

来源：[Wu 等，2020，PMID 32518267](https://pmc.ncbi.nlm.nih.gov/articles/PMC7283474/)，DOI [10.1038/s41467-020-16755-y](https://doi.org/10.1038/s41467-020-16755-y)。核对 Table 1、Figure 1、Results Par8、Methods Par35，并通过 [RCSB PDB](https://www.rcsb.org/) 的 6VRM、6VRN、6VQO 双链构建体序列确认 CDR3。原始响应均保留 URL、UTC 时间和 SHA256。

| TCR | 突变肽 HMTEVVRHC：测得 KD | 野生型肽 HMTEVVRRC：观察 | 野生型测量最高 TCR 浓度 |
|---|---|---|---|
| 12-6 | 1.1 µM | 未检测到明显 SPR 相互作用 | 328 µM |
| 38-10 | 39.9 µM | 未检测到明显 SPR 相互作用 | 232 µM |
| 1a2 | 16.2 µM | 未检测到明显 SPR 相互作用 | 153 µM |

共 6 行、3 个完整配对受体、1 项研究、1 对野生型/突变肽、1 种 MHC（HLA-A*02:01）。条件：BIAcore T100、25°C、pMHC 固定量 3000 RU、游离生物素空白流道对照、可溶性重组 TCR。构建体常数区含稳定性二硫键工程修饰，序列来源和 1a2 α 链更正均记录，不能直接等同临床原生细胞背景。

这里的 label=0 明确定义为 **该 SPR 条件与测试浓度范围内未检出**，不是 KD=0、KD=无穷大、绝对不结合、细胞不激活或临床安全。阴性记录的 KD 留空。这 6 行只用于外部 sanity check；`training_allowed=false`、`benchmark_allowed=false`。不能用单项研究的 6 行报告广泛泛化 AUROC、交叉反应准确率或临床风险。

## 泄漏排查与隔离

1. 对当前结构版本全部 306 个 native 结构记录（包括预测副本）检查配对可变区、CDR3 与保留受体，未发现重叠；910 张图所属的这些受体来源由结构记录覆盖。
2. VDJdb 关联表发现 12 行精确 CDR3 命中。采用更保守的同链 CDR3 core 相似度 >=0.9、至少 8 个残基的规则后，共隔离 **23 行**，生成 **199,362 行**的无标签新副本。保留原始 199,385 行表和每条隔离理由；这一过滤不创造任何结合标签。
3. 再次逐行检查 199,362 行新副本，保留集命中为 0。不能把这一结果外推为所有公共数据库均无序列暴露；未来新增来源仍必须检查。
4. `holdout_registry.json` 把整个受体及近似链保留用于外部核验。结构构建入口和 release 加载入口已增加强制拒绝检查；单链关联表也有对应的过滤与检查函数。这是防泄漏的保守排除规则，不是用单链相似性推断生物学交叉反应。

## 文件与复现

仓库中可审阅：

- `docs/evidence/p53_primary_curation.json`：6 行原文核查、完整可变区、CDR3/VJ、条件、序列来源哈希。
- `docs/evidence/external_sanity_panel.tsv`：冻结的 6 行面板。
- `docs/evidence/assay_evidence.tsv`：1,464 行待核查证据；label 全空。
- `docs/evidence/condition_conflicts.json`、`primary_review_queue.tsv`：冲突组和优先队列。
- `docs/evidence/evidence_manifest.json`：45 个来源响应的 URL/时间/哈希、依赖代码哈希和冻结产物哈希。
- `docs/evidence/association_reservation.json`、`sequence_exposure.json`、`integration_checks.json`：隔离与复现证据。

本地原始 API/XML/RCSB 响应位于 `data/external_evidence/20261016-repro/`；冻结产物在 `data/external_evidence/20261016-evidence-v3/`；过滤后的 VDJdb 关联表及隔离表在 `data/external_evidence/20261016-training-safe/`。大文件和原文快照不进入 Git；可审阅的证据表、原文核查和来源清单进入 Git。

```bash
cd /home/administrator/docs/meta-acdc-readiness
export PYTHONPATH="$PWD/src" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
PY=/home/administrator/docs/meta-acdc/.venv/bin/python
# 必须选择新目录；带已验证清单的下载可断点复用，不覆盖已有文件。
$PY -m meta_acdc.data.fetch_experimental --out NEW_SNAPSHOT
$PY -m meta_acdc.data.experimental_evidence --snapshot NEW_SNAPSHOT \
  --structural-records data/releases/20261016-structural-v2/structure_records.json \
  --pilot docs/evidence/p53_primary_curation.json --out NEW_EVIDENCE_RELEASE
$PY -m unittest discover -s tests -v
```

若远端原文/序列响应的哈希变化，旧核查文件会拒绝使用，须重新核对，不能只改哈希放行。过滤关联表的函数为 `meta_acdc.data.experimental_evidence.filter_associations(source, output, quarantine, registry_path)`，三个输入/输出路径必须明确且输出必须不存在。

实际做了第二次全新下载：122 组、531 导出、1,174 个实验、45 个响应，与首次发现快照生成的证据表和 6 行面板逐字节一致。更新后的测试为 **56 项，53 项通过，3 项旧夹具缺失跳过**；没有优化器训练。结构数据和历史产物保持原样。

## 放行边界与下一步

- hy4 已按用户要求暂时忽略。发布函数不再把一次性的 HTTP 402 写成所有未来数据版本的固定阻断原因。
- 本轮新增保留集守卫，旧结构 manifest 的代码哈希因此不再匹配；旧 manifest 和验收记录保留为历史，不在原文件上改动。下一次申请结构训练前，必须按当前代码重新校验、签发新版本 manifest，不能关闭哈希/保留集守卫。发布数据本身不自动批准或启动训练。
- 真实结合/激活任务仍缺少足够、跨研究、跨受体且实验终点一致的验收数据。接下来应按优先队列核查原文，扩充多个独立受体和研究；需要保留不等号与检测范围，解决条件差异，不能把当前面板扩写成一个大规模基准。
- 临床索引、单链 CDR3 关联继续只用于索引和候选发现；本轮没有新增临床安全证据。
