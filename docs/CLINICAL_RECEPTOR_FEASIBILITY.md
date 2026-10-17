> 2026-10-18 更新：[临床原文、工程构型和结局整理](CLINICAL_PRIMARY_CURATION_20261018.md)。临床资料仍按患者/组合保留；新SPR开发集51条，临床候选队列独立预留，旧版本数据量和结果仅作历史。

# 真实实验训练与临床独立受体：可行性核查

## 训练能否使用现有数据

可以开发新的研究模型，但须保留各自终点及独立验收来源。

- SPR：当前冻结亲和力数据可用于连续 log KD／pKD 回归。已有 91 条排除外部重叠后的配对观测与受体／文献连接块切分；不是 91 个独立受体，也没有自动补足“非结合”二分类。
- MIRA：功能观察可支持细胞激活模型，未激活不能改叫 SPR 不结合。原 10,000 条面板已用于外部打分；若训练或反复方法选择使用它，必须把用途改为开发数据，并获取新的、未参与选择且受体家族隔离的验收来源。不能一边训练一边继续声称原外部验收独立。
- EGNN：需要为实验观测建立身份准确、生成方式一致、QC 可比的结构输入。晶体阳性与预测／嫁接阴性混用会引入来源标签捷径。当前结构原子审查失败，不能直接恢复旧 EGNN 重训。

## 临床来源确实存在

新增查到并下载 [TRAIT Therapeutics](https://pgx.zju.edu.cn/traitdb/therapy/) 的[官方公开临床表](https://pgx.zju.edu.cn/traitdb/download/)，原 ZIP 保存在本地 evidence 目录，来源哈希见 [snapshot](n_review/trait_source_snapshot.json)。这比只有试验事件总表更接近“受体—抗原—临床”连接。

对工作簿做了只读检查：72 个有受体名称的行；51 行含可直接解析的双链 CDR3，去重为 48 个 CDR3 配对候选。与现有结构、亲和力及保留面板作保守的任一链筛查时，51 行中有 8 行触发已见链重叠；未触发并不证明家族独立。未把含括号工程突变说明的字符串自动当成实测最终序列。表内没有完整 alpha/beta 可变域序列列。48 不是独立临床受体数，更不是 48 个可独立归因的毒性样本。[可复现逐行身份／暴露筛查](n_review/trait_identity_exposure_audit.json)。

当前已经发现需要隔离的例子：Excel 第 6 行 gp100 项标为 DMF4，并复用了 MART-1 行的 CDR3 对；[PMID19451549](https://pubmed.ncbi.nlm.nih.gov/19451549/)明确区分 MART-1 的 DMF5 与小鼠来源 gp100(154) 受体。这条不能原样入库。第 9 行 MAGE-A3 A118T 神经毒性案例缺 beta CDR3，不能凭它声称已有完整临床双链。

有多于 titin/a3a 一个受体来源的线索，例如 MART-1 DMF5、HPV16 E7、个体化新抗原 TCR、KRAS G12D。是否与训练集“真正独立”必须逐条回到原文和补充表核对完整工程序列、亲本谱系、同患者多受体共给药、同队列重复报道，以及是否已在 ATLAS／结构／功能开发数据中暴露。

[个体化受体替换临床研究](https://pubmed.ncbi.nlm.nih.gov/36356599/)是可继续追溯的原始来源，但多个受体与同一患者结局不能拆成互相独立的标签。MAGE-A3/MAGE-A12 神经毒性的[原始研究](https://pubmed.ncbi.nlm.nih.gov/23377668/)与 HLA-A1 titin/a3a 机制不同，可以构成另一来源线索；仍需补全工程构建体及来源隔离。

## 放行需要形成的受体级证据表

每条保留：完整 alpha/beta 变量域、工程突变、构建体／亲本谱系 ID、肽和 HLA、真实测量终点与条件、原文定位、试验／队列／患者映射、共给药、毒性类型和因果等级、训练来源暴露，以及预先冻结的训练／验证／外部用途。

临床毒性须区分靶上组织毒性、分子脱靶、CRS、治疗预处理毒性和未归因事件。未报告事件或数据库写 None 不自动生成安全阴性。先建立受体级可追溯案例面板，可以做有限范围盲测；现有公开数据尚未完成足以训练和验收临床危险概率的证据连接。


复现：`PYTHONPATH=src python scripts/audit_trait_clinical.py --zip data/external_evidence/20261017-trait-clinical/Therapeutics.zip --structures data/releases/20261017-structural-v4/structure_records.json --affinity docs/database_binding/affinity_dataset.json --registry docs/database_binding/combined_reservation_registry.json --out NEW_TRAIT_AUDIT.json`。当前结果是来源初筛，不是已经验收的数据集。来源工作簿未改写或重导出。
