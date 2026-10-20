# 新SPR独立预留池：部分完成，未打分

本池保留未来验收，不进入训练。当前只有1个受体家族、1篇研究、1个来源组件和1个连接块；没有满足至少3个来源和3个连接块的原门槛，`evaluation_allowed=false`。未与已消费测试集拼接，也未重洗切分。

## 已核实的新增原文

Chikata等，Nature Communications2025，[原文DOI](https://doi.org/10.1038/s41467-025-64768-2)，PMID41198628。KI-890的TRAV24*01/TRBV27*01受体，完整配对可变区来自[9K2R](https://www.rcsb.org/structure/9K2R)，CDR3与SI Table1对应；不是按胚系推算。alpha/beta分别在明确常量区DIQNP/EDLNK之前截断。

| pMHC（HLA-C*12:02） | KD | 重复 | 分析角色 |
|---|---|---|---|
| ILKEPVHGVYY（WT） | 0.65±0.09µM | 3次独立实验 | primary候选 |
| ILKEPVHGAYY（9A） | 2.77µM | 1次 | sensitivity_only |

Biacore3000，25°C，CAP芯片固定1200RU生物素化pMHC、流过可溶TCR，1:1动力学拟合。Methods把浓度印作0.09–1.4 M，而Fig2d明确标µM；已查看PDF保留冲突，不静默纠正原文。范围判据仅使用Fig2d明确标出的最高1.36µM，所以2.77µM保守退出主指标。动力学拟合超过测试浓度不自动证明数值错误；此处沿用既有保守外推排除标准。

家族审查对历史91条、临床预留、旧独立家族和结构v1–v7均无命中；家族定义为已有可变区/CDR3规则，不能称作V基因完全不重叠。新开发集排除了共享作者的ATLAS:485（3B）后19条维持原fold/block；源、家族和肽相似性重新检查。原alpha1参考配置此前使用过3B参与选择，因此本池不能被宣传为该参考模型的严格来源独立验证；未来必须用新登记的来源隔离协议。

## 无法补齐的候选

- HIV IW9：726_TCR3、726_TCR4、309_TCR3、309_TCR4四条SPR有KD1.53/1.92/6.53/7.98µM、BiacoreT200、25°C，但缺与这些编号一一对应的公开完整配对可变区。旧/新两个预印本全文及旧SI已查，未发现可恢复序列的GSE/PRJ/SRP编号。保留条件记录，不能拟造链对。
- UTP20 CL1/9X7R：有结构与序列，出版社正文不可公开取得；HTTP200接口只给元数据，未核到精确KD和SPR条件。
- 9SOL/9SOK、9D95、9YIR/9Z51/9Y7Y：尚未找到与结构对应的完整定量SPR原文及方法。PDB注释提及SPR不等于测量证据。
- RAC1-5934、AS19.2、JDI、C3-QW9：已有受体或近链家族暴露；新增结构、工程变体或另一篇文章不能恢复独立性。7PB2命中历史7pb2_0/1，7R80命中7r80。
- KRAS9a/9d：临床受体已预留，不能挪作新的独立验收。
- T18A/T11A/T7A：实际Octet QKe BLI，不能因为章节写SPR就纳入严格SPR池。
- crosslinked TCRA–E：已知亲和力与温度，但公开资料仍缺仪器、方向和拟合，不能按同实验室惯例补填。

完整逐项处置见`candidate_disposition.json`，IW9隔离数据见`quarantined_IW9_measurements.json`。原始下载和失败/权限响应保留在本地`data/external_evidence/20261018-fresh-spr-v1`；哈希、URL和体积见`raw_sources.json`，不把响应成功当作原文已核实。主要PDF/XML/补充材料受`required_sources.json`绑定。

校验：`python -m meta_acdc.data.fresh_spr_reserve`。`require_acceptance_ready()`明确拒绝此部分池的评分；未来训练协议必须调用`assert_development_safe`检查预留源、受体和肽。既有冻结训练入口仍因旧测试已消费而拒绝再拟合，未被绕过。新预留校验器并非自动替换所有历史训练脚本的全局入口。
