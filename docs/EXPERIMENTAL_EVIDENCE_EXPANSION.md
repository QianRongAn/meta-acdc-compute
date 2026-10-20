> **2026-10-17 update:** [Scoped structural retraining is now released](RETRAINING_DECISION_20261017.md); binding and clinical claims remain blocked. The release/guard state described below is historical.

# Experimental evidence expansion and retraining boundary — 2026-10-16

本报告更新 `EXPERIMENTAL_EVIDENCE.md` 的首批单研究面板；旧表、原始下载和历史结构版本均保留。hy4 按用户要求暂缓，不作为当前阻断原因。本轮没有训练、预测提交或新性能指标。

## 已完成的证据核查

`evidence/primary_measurements.json` 保存 **23 条观察、8 个具体受体构建体、4 篇原始研究**：21 条 SPR（10 条测得结合、11 条条件性未检出），另有 2 条 TNFα 细胞激活观察。8 个构建体不代表 8 个独立家族，JDI 亲本与增强版本共享来源。所有 `training_label=null`，只保留为外部核验；不是训练集或总体性能基准。

| 原始研究 | 具体构建体 | 核查结果及界限 |
|---|---|---|
| [p53，PMID 32518267](https://pmc.ncbi.nlm.nih.gov/articles/PMC7283474/) | 12-6、38-10、1a2 | 保留既有 6 条条件明确的 SPR 观察；野生型未检出的测试上限逐受体记录 |
| [KRAS，PMID 36088370](https://pmc.ncbi.nlm.nih.gov/articles/PMC9464187/) | JDI、JDIa41b1 | 亲本 1 条测得结合及 7 条表中 NB；增强版突变肽和野生型均测得结合。不能继承亲本阴性；部分 NB 缺逐肽最高浓度，明确留未知 |
| [NRAS，PMID 39576860](https://pmc.ncbi.nlm.nih.gov/articles/PMC11584006/) | N17.1.2 | 两种突变肽测得结合，野生型至 200 μM 未检出。TRDV1 用于 α 链，不将该 αβ 受体误排为 γδ |
| [PIK3CA，PMID 35484264](https://pmc.ncbi.nlm.nih.gov/articles/PMC9117146/) | TCR3、TCR4 | 两条突变肽 SPR 亲和力；不凭方法提到 WT 就补出 WT SPR 阴性。TCR4 的两条 TNFα 结果另列，不能当成内在结合标签 |

双链均逐实体从 RCSB 沉积构建体核对；保留表达构建体/常数区工程说明、可变区、MHC 等位基因、方法、实验条件、原文图表位置、URL 和 SHA256。原文与序列文件更改会拒绝验证。检测不到不是 KD=0、无穷大或永久不结合；细胞不响应不是临床安全证据。

KRAS 原始科学补充 PDF 和 Source Data.xlsx 经 Europe PMC 官方补充包取得，保留 ZIP/PDF/XLSX 原字节。`supplementary_replicates.json` 精确记录表名、单元格、原始重复测量和单位；重复测量不增加独立样本数，论文四舍五入值与原始重复值分别保存。同行评审信不当作实验数据。JDIa9bwt、JDIa96b35、工作簿仅有 Mutant 编号的变体缺完整双链，未借用其他构建体序列；肽池 NB 不拆成各肽阴性。小鼠 TCR 和 TCR 模拟抗体也未混入人源 αβ 面板。

## 原先 20 条优先队列的逐条裁定

见 `evidence/priority_adjudications.tsv`，每行保留原证据 ID 和理由：

- 9 条 CF34：原文全文/补充的具体条件尚不可取得，保留数据库观察，严格面板不收。SPR 截尾界限与多聚体结合不能强行合并。
- 2 条 DMF4 SPR：两种肽均有结合证据；原始 Figure 7 的本地原文快照不足，暂定而不收入完整溯源面板，后续文献重复引用不计作独立研究。
- 9 条 DMF4 细胞/多聚体：从内在结合面板排除。CD4/CD8 背景、整蛋白/肽刺激和多个构建体合并记录会影响观察；不把 CD4 无应答改写成普遍非结合。

这是完成逐条裁定，不是宣称 20 条都符合纳入条件。完整 1,464 行数据库台账及 105 个条件冲突组仍是待分层核查来源；本轮没有声称全部阅读了其原始研究，也没有投票消除冲突。

## 泄漏隔离及新的结构版本

外部受体注册表扩为 8 个具体构建体，按对应 α/β 链及 CDR3 的保守近似匹配排除；这不是以单链相似推断交叉反应。

原结构 v2 含 6 个 JDI 相关 native 结构：`7pb2_0`、`7pb2_1`、`7ow5`、`7ow6` 及两个预测副本。新建 `data/releases/20261016-structural-v3/`，删除整个命中家族、其派生图以及采用保留家族作为肽供体的嫁接图，共排除 22 张图。保留图的数组和原有 fold 不变，索引仅重排映射；原 v2 及其 manifest 未覆盖。

新版本 **888 张图、300 个 native、300 个 displaced、288 个 graft、86 个家族**；train/validation/calibration/test 为 **537/157/67/127**。机器检查通过，包括图特征、坐标边/RBF、源文件、逐模型 QC、家族切分、供体同 fold 和保留集零命中。发布和加载入口也强制检查保留集。

从原始无标签 VDJdb 关联表 199,385 行重新过滤，隔离 106 行，剩余 199,279 行再次逐行检查保留受体命中为零。隔离表、行号及原因保存在新版本，未给剩余关联添加任何结合标签。现有 4,845,721 条 IEDB bulk 为 MHC ligand 终点，继续单列，不能冒充 TCR 阳/阴性训练集。

`structural_v3_validation.json`、`structural_v3_reservation_report.json`、`structural_v3_release_manifest.json` 记录验证与哈希。manifest 的 `code_commit` 是构建时 Git 基点，实际工作树源文件由 `code_files` SHA256 固定；不能仅凭该基点复现并忽略文件指纹。加载新版本完整检查成功；默认训练入口仍拒绝。原工作区干净，历史六项关键数据哈希不变，见 `final_integrity.json`。

## 可复现步骤

复用原 `.venv`，不安装依赖：

```bash
cd /home/administrator/docs/meta-acdc-readiness
export PYTHONPATH="$PWD/src" OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
PY=/home/administrator/docs/meta-acdc/.venv/bin/python
$PY scripts/curate_primary_panel.py --out NEW_PRIMARY_OUTPUT
$PY -m meta_acdc.data.reserve_release \
  --source data/releases/20261016-structural-v2 --out NEW_STRUCTURAL_OUTPUT \
  --registry docs/evidence/holdout_registry.json
$PY -m meta_acdc.data.validate_release --release-dir NEW_STRUCTURAL_OUTPUT --publish
$PY scripts/prepare_external_structure_inputs.py --out NEW_STRUCTURE_INPUTS
$PY -m unittest discover -s tests -v
```

原文快照须位于各行 `sources`/`sequence_sources` 指定路径并与哈希一致；大文件不进入 Git。整理脚本在新的空目录复现全部 6 个证据产物，逐字节一致。最终补充下载清单也固定在 `expanded_evidence_manifest.json`。测试 **64 项，61 通过，3 项旧夹具缺失跳过**，无优化器训练。

## 实际尚缺的外部条件

已为 21 个 SPR 受体—肽组合生成配对一致的五链输入：α、β、肽、原文对应 MHC 重链、B2M；链源、输入哈希和准备状态见 `external_structure_jobs.json`。正负观察须用同一冻结预测版本、seed 和模板/MSA 条件，经 QC 后才能评分；不能比较晶体阳性与预测阴性来报告性能。作业均 `prepared_not_submitted`，尚未占用外部服务或重复派发。

当前 PATH 无已配置 AF3/ColabFold 入口；检查到 GTX 750 Ti 2 GB。已完成输入准备，**没有模型输出**。还需要可用结构预测后端，以及更多跨受体、跨研究、条件一致的测量；现有 case-selected 外部面板不能支持广泛特异性或临床风险主张。缺完整原文或序列的候选继续隔离，不通过猜测补齐。临床索引与共享 CDR3 关联仍只作索引和候选发现。训练放行保持关闭；代码修复、机器验证和科学性能验证是不同验收事项。
