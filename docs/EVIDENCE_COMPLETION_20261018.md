# SPR原文续查、工程参考序列和外部验收资料（2026-10-18）

本轮从剩余26条SPR候选继续追溯，新增7条达到当前严格池条件，确认1条肽—数值错配并隔离。**当前严格池23条；训练15、验证8、原测试0，原有fold和连通块均保留。没有重训、评分验收集或启用调度。** 旧16条、51条及更早版本均保留。

## SPR逐条处置

[51条更新清单](evidence_completion_20261018/review_51.tsv)与[完整证据账本](evidence_completion_20261018/review_ledger.json)保留原始行和新判定。新追回5项研究，覆盖原来8条未核查记录；其余18条、13项研究仍未拿到足够原文，其中D30补充表已核实136±7uM，但仍缺完整Methods/温度。

|本轮记录|结果与来源|
|---|---|
|ATLAS:201，SB27，18.7uM|通过；[PMID22343629](https://pmc.ncbi.nlm.nih.gov/articles/PMC3320977/)，Results与SPR Methods，25°C|
|ATLAS:717，HCV1406，16uM|通过；[PMID28572406](https://pmc.ncbi.nlm.nih.gov/articles/PMC5474831/)，Figure1C及SPR/DSF Methods，25°C|
|ATLAS:502/505，TK3，93.2/25.4uM|通过；[PMID24759101](https://pmc.ncbi.nlm.nih.gov/articles/PMC4059114/)，Figure3A/D及SPR Methods，25°C；本轮未独立核实数据库误差条，派生数据不沿用它们|
|ATLAS:485，3B，81.8uM|通过；[PMID24600035](https://pmc.ncbi.nlm.nih.gov/articles/PMC3962895/)，Figure2A、Results和SPR Methods，25°C；误差值同样留空|
|ATLAS:644/645，C25/C7，4.7/5.1uM|通过序列/实测亲和力池；[PMID26429912](https://pmc.ncbi.nlm.nih.gov/articles/PMC4661422/)，Figure1B/D及SPR Methods，25°C。原文测HLA-A*0201，但数据库把沉积MHC注为A*0202；已禁止把这些KD直接用于对应结构训练，待核对MHC序列|
|ATLAS:914，TK3–Q5，2.2uM|隔离：原文2.2uM明确指HPVGEADYFEY，并引用旧研究20566715；本行却是HPVGQADYFEY/4PRP。没有替Q5编造纠正值，也不把重复引用算新测量|

当前51条处置：23条可溶性SPR池、18条原文条件仍不足、4条温度不明、2条肽—数值错配、1条完整TCR-CD3格式分存、1条最初测量待追溯、2条方法/引用冲突。

**4条原有温度缺口仍未解决：NYBR1、WT1-alpha7beta2、LS01、LS10。** 数值已有来源；这次补查JCI正文PDF、补充材料及LS论文可见资料，仍未找到这些SPR实验的明确温度。不能用MD的298K、细胞培养温度或仪器惯例补填。逐条见[温度续查](evidence_completion_20261018/temperature_followup.json)。

## 临床工程资料

**WT1-C4新增完整公开参考编码盒：1809nt，602aa，beta–P2A–alpha。** 来源是[US20160083449A1](https://patents.google.com/patent/US20160083449A1/en)的SEQ96（PDF第80页/印刷第62页）。逐行转录后的DNA翻译与另外列出的蛋白SEQ14（PDF第45–47页）完全一致，核对了CDR3及全链alphaT177C/betaS188C，对应临床文献恒定区编号T48C/S57C。保存[序列与证据](evidence_completion_20261018/C4_published_reference_cassette.json)及[FASTA](evidence_completion_20261018/C4_reference_insert.fasta)。

它与临床原文的克隆、链顺序、半胱氨酸对和pRRLSIN设计相符，**仍不是临床批次序列证书，也不是完整载体序列**；不使用专利中的C4-DLT亲和力增强变体替代临床C4-WT。专利早期病例与后来的12人队列存在重叠，不额外累计独立患者。

新追回DMF4原始试验全文（PMID16946036），确认当时是**MSGV1AIB alpha-IRES-beta**，不能替换成后来P2A构型。试验17人、2例客观缓解包含淋巴清除及IL2背景，不转成受体单独疗效标签。还追回[SPEARHEAD1机构库全文](https://discovery.ucl.ac.uk/10213399/1/nihms-1984642.pdf)（PMID38554725），保留52人、19个客观缓解、37例CRS；它延续同一afami-cel受体，不能算新的独立受体。

E7新取得[临床原文的83页补充材料/方案](https://www.nature.com/articles/s41591-020-01225-1)，确认MSGV1完整载体7310bp、PG13包装、鼠恒定区、工程半胱氨酸和alpha跨膜区疏水改造、beta–furin2A–alpha示意。此前的1800nt公开插入片段因此有更明确的临床设计联系，但没有补造整个7310bp载体或批次证明。方案对胸腺选择的安全性推理没有当作“无脱靶”的实验结论。

仍缺：Foy最终工程盒（TCR919序列冲突仍隔离）、p53最终MSGV1/鼠恒定区/P2A编码盒、1D3HMCys最终插入序列、KRAS9mer/10mer临床构型与精确序列的完整对应、各临床批次序列证明。另有13篇原临床来源全文未取得。见[明确缺口](evidence_completion_20261018/clinical_remaining_gaps.json)，获取失败不等于作者从未公开。

## 外部实测验收资料

没有移动原训练/验证行来生成测试集，也没有评价这些面板。

|分层|规模|独立性与用途|
|---|---|---|
|常规外部SPR核心|4项测量、3个受体家族、2项研究|N17.1.2的两项KD及TCR3/TCR4各一项；来自此前明确预留、已记录亲和力模型未使用的外部面板，**不是本轮新发现数据**。保留原值、25°C、配对序列和来源，冻结供后续锁定模型验收|
|新收集的格式迁移挑战|6项配对测量、1个工程受体家族、1项研究|[HBV/HLA-E原文](https://www.nature.com/articles/s41467-024-54378-9)与SourceData Figure2c原始重复测量；a09b08 ImmTAV、3种肽、2种HLA-E。以M为单位的3次重复正确转换为uM，不把重复或六种配对计成独立受体|

新HBV数据不在已有ATLAS1174行中；已检查本地v1–v7结构历史、原91条亲和力记录、MIRA/原外部保留和临床保留家族，没有检出对应链/近CDR重叠。检查使用alpha/beta到历史保留格式的显式转换，缺失配对身份即拒绝检查，避免空序列产生虚假的“无重叠”。[逐源独立性审计](evidence_completion_20261018/external_independence_audit.json)保留输入哈希。

HBV数据的HLA-E、亲和力增强和抗CD3融合格式不同于常规可溶性TCR，单独验收；半衰期20小时可能受测量窗限制，保留作者拟合值与局限，不宣称超高亲和力被精确确证。PDB只覆盖配对TCR部分，不包含完整抗CD3融合域；E0101配对没有冒用E0103作为其复合物结构。

**当前新增资料仍不足以宣称广泛的独立验证，更不能放行临床风险模型。** 常规核心只有2项研究、4个观测，新面板也只有1个受体家族。冻结[外部验收协议](evidence_completion_20261018/external_test_protocol.json)要求以后锁定模型后分别报告pKD误差、训练均值对照和研究/受体相关性；禁止混池凑大样本、基于面板调参或把未检测到当临床安全阴性。原23条发布仍为inspection_only，训练入口拒跑。

## 复现与保存

当前本地发布：`data/external_evidence/20261018-primary-reviewed-affinity-v4`。活动代码要求本轮清单、C4/HBV新增保留依赖和原始来源哈希；旧manifest不能绕过更新后的审计。

```sh
PYTHONPATH=src python -m meta_acdc.data.spr_primary --snapshot data/external_evidence --out NEW_RELEASE_DIRECTORY
PYTHONPATH=src python -m unittest discover -s tests -v
```

本轮177项测试：174通过、3项因既有结构fixture缺失跳过。覆盖错误单位、TK3错配、C4翻译/链顺序/批次声明、历史身份字段转换、新受体保留、原始工作簿单元格和旧冻结文件不变。第三方全文/补充材料保留本地原始快照，Git保存派生事实、定位、URL、SHA256与失败原因；原始工作区、历史结果/日志、原临床18个验收候选及旧冻结目录保持不变。
