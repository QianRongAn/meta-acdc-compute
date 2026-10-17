> 当前证据状态：[完整工程DNA与SPR第四轮续查](SEQUENCE_FOLLOWUP_20261018.md)。新增1D3HMCys公开研究质粒6338nt/完整插入1797nt→598aa，临床批次等价性未核实。TCR156新增5条真实SPR，四聚体固定相格式待澄清，单独保留；严格核心仍7条未评分。D1/D6全变量链及原10条条件、2条表值、2条温度仍缺；严格池31条，测试仅2个独立块，训练继续阻断。

> 历史版本记录（已由上方续查更新）：[临床原文、工程构型和结局整理](CLINICAL_PRIMARY_CURATION_20261018.md)。临床资料仍按患者/组合保留；新SPR开发集51条，临床候选队列独立预留，旧版本数据量和结果仅作历史。

> **2026-10-17 补充：结构任务旧放行已撤销，见 [新审计](N_REVIEW_20261017.md)。新增 [TRAIT 临床受体来源核查](CLINICAL_RECEPTOR_FEASIBILITY.md)。本文所述亲和力／功能实验保留原版本和结论。**

# 数据库补齐真实实验数据与验证 — 2026-10-17

本轮实际执行了数据获取、序列/终点核查、来源隔离、从头拟合和一次独立功能面板评分。**现在有真实实验数据可做亲和力研究，并有公开功能验收面板；当前模型仍未获得真实结合或临床风险放行。** 没有把不显著或缺失记录转成阴性，没有安装依赖、申请账户或发送数据请求。

## 取得的数据及各自用途

| 来源 | 本轮实际处理 | 能支持什么 | 不能支持什么 |
|---|---|---|---|
| [ATLAS](https://atlas.ibbr.umd.edu/web/downloads.html) | 原始 1,174 条记录，原值/单位/不等号/误差/PMID/条件/突变保留 | 实验 SPR 亲和力 pKD 回归探索 | `n.d.`、弱结合和 >200 μM 不能自动标成不结合 |
| [TCR3d](https://tcr3d.ibbr.umd.edu/downloads) 与 RCSB | 链注释及 168 个精确 PDB 的配对序列快照 | 具体构建体与肽身份核查 | 没有 KD 的结构不能补出亲和力；模板不能冒充突变体实测结构 |
| [Adaptive 的 IMMREP25 公开发布](https://doi.org/10.6084/m9.figshare.30491582) | 10,000 条、1,000 个配对 TCR、20 种病毒肽、2 种 HLA，MIT 许可 | MIRA 池实验的激活/条件性未激活外部检验 | 不是直接 SPR 非结合，也没有临床毒性结局 |
| [ClinicalTrials.gov](https://clinicaltrials.gov/) | 4 项试验、8 个事件组、690 条“组—事件”观察 | 带人数分母、时间窗和报告条件的回顾性安全审计 | 不是 690 名独立患者；不能标记某条脱靶肽致死或安全 |

ATLAS 是实验数据库整理级证据，本轮没有逐篇人工重读全部 54 篇最终纳入研究，不能将它们写成逐篇原文重新验收完成。数据异常、截尾、模糊物种/MHC、非 SPR 方法、突变 MHC、缺配对/不明确链、肽不符、重复或冲突测量分别隔离，原记录不删除。仅采用 `true_PDB` 所对应的实际构建体，未从模板猜测突变全序列。ATLAS 的 TCR 名称可能是亲本或别名；身份按实际配对序列固定，不按名称继承标签。亲和力范围同时包含弱结合与工程增强构建体，不是天然受体总体抽样。

## 亲和力数据与防泄漏

首次严格筛选得到 129 条精确 SPR 记录（74 篇研究）。再对新增 IMMREP25 面板进行保守近似链排查，命中 38 条，因此本轮新模型只用 **91 条、54 篇研究**。原先 23 条原文核查观察对应的 8 个受体也继续保留。联合保留表 `database_binding/combined_reservation_registry.json` 共 1,008 个具体构建体；它适用于本轮新的亲和力实验，不覆盖旧结构版本原有范围。

冻结训练/验证/测试为 **55/18/18**，连接块 **20/6/6**。同一论文、受体—肽—MHC 对、近似对应 α 或 β 可变区以及相似度 >=0.8 的肽均不能跨 fold。分组不看目标数值。超过检测范围的 KD 保留界限，不参与普通精确值回归；缺失不当成 0。目标 `pKD = -log10(KD[M])`，与结合概率、激活概率和临床危险概率分别对待。

新增外部面板也与 structural-v4 的 **79 个 native 结构记录**命中保守近似链规则。**旧结构 checkpoint 或未来在未排除此重叠的 v4 上训练的模型，不能声称在该功能面板上完成独立验证。** v4 只针对其原先定义的结构诱饵任务放行。本轮使用完全从头拟合的序列/实验条件基线，未加载任何结构或历史模型权重。亲和力 91 条与功能面板 20 种肽之间的显式 >=0.8 相似度复查命中为 0；该肽复查在评分后记录，不冒充预注册前完成。功能受体的联合保留和排除则发生在拟合前。

## 完成的拟合与实际结果

基线为双链及肽的单/双残基频率、TCR×肽交互、MHC、仪器和温度信息上的 Ridge 回归；字典只在训练 fold 拟合。alpha=1/10/100 仅按验证 RMSE 选择，原始实验重复的权重不膨胀为独立样本。

- 129 条试运行：测试 pKD RMSE 0.758，训练均值基线 0.729，未优于均值基线。
- 排除新增功能受体后，91 条版本：测试 RMSE 2.262，均值基线 2.272，提升很小、绝对误差仍大。内部比较条件通过不等于预测任务放行。v1/v2/v3 为来源隔离/依赖指纹/消费者导出修订的连续实验，**不是三次独立盲测复现**；保留各次历史结果。
- 对 10,000 条功能观察进行一次外部排名检查：平均肽级 **AUROC 0.5001**，平均 AP 0.1197。按肽重采样的 AUROC 区间约 **0.4586–0.5431**。该区间只描述固定 TCR/平台下的肽间变化，不是独立受体家族/研究或临床置信区间。

外部评分是在固定假设的 SPR 参考条件下计算的 pKD 分数与 MIRA 激活标签的关联，不是激活概率校准。该面板没有用于拟合或调参；评分后未按结果改模型或阈值。这是当前弱序列/条件基线的阴性结果，**既不验证原 EGNN，也不能推出序列路线已穷尽或结构路线无效**。未报告“真实结合模型验收通过”。

## 临床数据库实际能验证到哪里

下载并冻结 NCT03070392、NCT01273181、NCT01350401、NCT01352286 的当前官方 API 响应。每条事件保留事件组、术语、人数、分母、报告词表、系统性评估说明、报告阈值、时间窗与来源哈希；缺字段保持未知，不补 0。全因死亡不改写为药物致死，多种事件或试验的人数不求和成独立患者。

核查发现需要版本化处理的身份冲突：[titin 原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC6002776/) 曾以 NCT01350401/NCT01352286 描述 MAGE-A3 a3a 试验，而当前注册条目/结果是 NY-ESO-1 c259。不能按同一个编号把当前结果套到历史 a3a/titin 构建体；逐条差异记录在 `clinical_construct_link_audit.json`。相关原始研究通过官方网页核对，但本轮 Europe PMC XML 返回 HTTP500，未伪造本地原文 XML 哈希。

这些表没有把患者的具体构建体双链、脱靶肽、暴露及组织因果结局连起来。它们适合核对已知风险证据和实验设计，**不足以验收分子级临床风险模型**。不同 TCR-T、可溶性 TCR–anti-CD3 和对照药物也不能混成一个机制。未把阴性 AE 行、未报道事件或无病例报告的肽当成“临床安全”。

[FAERS/AEMS](https://www.fda.gov/drugs/drug-approvals-and-databases/fda-adverse-event-reporting-system-faers-database) 可作为药物安全信号来源，不能单独提供发生率、分子脱靶因果或无风险阴性，因此没有用它制造训练真值。

## 其他来源的纳入决定

- [TCRvdb](https://github.com/schumacherlab/TCRvdb)：2026-10-01 版本涉及 5 种表位的再验证，官方要求申请并限制未经许可再分发。本轮只取得公开说明和申请页面；未发送申请、未取得受限数据、未将非显著统计结果一概转成绝对非结合。
- [IMMREP23](https://github.com/justin-barton/IMMREP23)：官方说明测试阴性由交换 TCR 构造，未当作实验非结合纳入。
- [CROSSTCR](https://i.uestc.edu.cn/CROSSTCR/download.html)：取得版本元数据；可继续用于交叉反应候选发现，但没有把其关系条目当作患者结局或非结合真值。

## 对“五关”的当前结论

| 关口 | 本轮解决到的程度 |
|---|---|
| 定义目标 | 分开 SPR 亲和力、MIRA 激活与临床事件，已落实 |
| 可靠标签 | 已取得可建模的实验亲和力及功能观察；广泛、条件充分的物理非结合标签仍不足 |
| 独立验收 | 本轮新模型的来源排除、冻结连接块和未参与训练的外部面板已落实；不能沿用旧结构模型声称独立 |
| 一致输入 | 本轮序列基线采用一致构建体/特征；EGNN 的正负配对结构仍需统一预测与 QC，不能混用晶体阳性和预测阴性 |
| 有效性验收 | 已实际完成比较和外部检查；当前结果不支持泛化放行，临床也不放行 |

取得数据解决了“没有实验依据可研究”的一部分问题，并不保证模型通过。下一次模型开发只能依据训练/验证集进行；这次外部面板已经使用过，应保留为固定回归检查，新模型最终确认仍需另一个未用于方法选择的外部来源。

## 文件与复现

`docs/database_binding/` 包含精确配对数据、冻结切分、联合保留表、来源索引、清单、两次亲和力比较、外部功能评分、临床组摘要和身份冲突审计。大体积原始响应和 10,000 条规范化功能行保存在本地 `data/external_evidence/`，没有覆盖旧数据/日志/checkpoint。源文件、数据文件和消费者代码变更均拒绝使用冻结亲和力版本。测试 **85 项，82 通过，3 个旧夹具缺失跳过**。

```bash
cd /home/administrator/docs/meta-acdc-readiness
export PYTHONPATH="$PWD/src" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=/home/administrator/docs/meta-acdc/.venv/bin/python
# 可选：按来源索引恢复到新的目录。上游版本已变时拒绝，不静默更新哈希。
$PY scripts/fetch_database_snapshot.py --source-index docs/database_binding/source_index.json --out NEW_SNAPSHOT_ROOT
$PY -m meta_acdc.data.functional_benchmark --snapshot data/external_evidence/20261017-public-binding --out NEW_FUNCTIONAL_RELEASE
$PY -m meta_acdc.data.atlas_binding --snapshot data/external_evidence/20261017-databases \
  --registry docs/database_binding/combined_reservation_registry.json --out NEW_AFFINITY_RELEASE
$PY -m meta_acdc.models.affinity_baseline --release NEW_AFFINITY_RELEASE --out NEW_MODEL_RUN
$PY scripts/validate_affinity_on_activation.py --affinity-release NEW_AFFINITY_RELEASE \
  --model-run NEW_MODEL_RUN --functional-release NEW_FUNCTIONAL_RELEASE --out NEW_EXTERNAL_CHECK
$PY -m meta_acdc.data.clinical_database --snapshot data/external_evidence/20261017-databases --out NEW_CLINICAL_RELEASE
$PY -m unittest discover -s tests -v
```

所有输出必须使用新目录；上述命令用于复现，不是反复看测试结果调参的入口。来源许可/引用以 `source_admission_decisions.json` 及官方页面为准。
