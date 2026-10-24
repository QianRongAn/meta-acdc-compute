# 完整αβ受体跨PDB、跨研究来源核验（2026-10-24）

本轮找到3组完整沉积αβ变量链完全一致的跨论文受体，其中1组作者不重叠；找到3种同V基因组合的不同变量链受体候选对照。**仍没有同时满足来源平衡、困难对照及原子QC的评估单元，因此没有冻结新训练/评估协议，没有训练模型，也没有产生新的AUROC。** 独立SPR、功能与临床验收集均未访问；所有历史发布保持原样。

## 搜索范围和身份准入

从TCR3d当前快照的799个PDB获取RCSB完整聚合物序列、主引文、作者及结构元数据。以旧24个完整双链身份的48条唯一链进行精确序列搜索，补15次档案搜索（含仅用于发现候选的80%相似度搜索）；逐项核验176个候选PDB，覆盖原107个天然PDB。检索并不保证覆盖世界上所有合格结构，也不把局部100%比对等同于全长一致。

最终89个链对单元、44个PDB、29个完整双链变量序列ID通过身份层核验。一个PDB的多个晶体副本不能增加受体或研究来源数。身份层通过不等于原子QC通过，更不等于结合/非结合标签成立。

沿用此前完整变量域准入条件：STCRDab显式双链配对及编号、原始SEQRES/聚合物完整构建体、IMGT V/J端部证据必须一致；保留工程突变，不补齐缺失端部，不凭名称、单链CDR3或相似度合并受体。仅当同一RCSB polymer entity的完整构建体逐字相同时，允许从该实体另一晶体副本继承编号注释，并记录证明。

修正一个来源解析假阴性：原始COMPND写“42F3 ALPHA”时，旧解析器因缺少TCR字样返回unknown。新核验器结合独立TRAV/TRBV证据确认显式ALPHA/BETA角色；若原文明确相反或歧义，仍拒收。旧解析结果及修正结果同时保留。

## 跨论文同受体证据与原子QC

| 完整变量链身份 | PDB与研究来源 | 来源层结论 | 当前结构阻断 |
| --- | --- | --- | --- |
| YLQ7 / YLQ-SG3 | 7N1F：[Wu等，2022](https://doi.org/10.1038/s41467-021-27669-8)；7RTR：[Szeto等，2021](https://doi.org/10.3390/cells10102646) | 完整αβ变量链完全相同，作者不重叠；两篇Methods分别说明表达、结晶与独立衍射采集 | 7N1F通过；7RTR链E第28位GLU实际缺CG/CD/OE1/OE2，不能组成完整QC通过对 |
| 工程α24β17 | 4JFE/4JFH：[2013](https://doi.org/10.1074/jbc.M113.464560)；6EQB：[2019](https://doi.org/10.1002/eji.201948085) | 跨论文相同完整变量链，作者有重叠；不能称独立团队复现 | 6EQB通过；4JFE被固定粗筛拒收；4JFH是游离TCR，缺完整pMHC |
| 42F3 | 3TFK/3TJH/3TPU：[2011](https://doi.org/10.1016/j.immuni.2011.09.013)；4MS8/4MVB/4MXQ/4N0C/4N5E：[2015](https://doi.org/10.1038/ni.3310) | 跨论文相同完整变量链，作者有重叠；全文Methods未取得，不能声称完成其核读 | 紧凑MHC构建体缺β2m，超出现有完整class-I图的范围；不等于生物学无效 |

2019年原文Methods明确：**6EQA是native MEL5，6EQB是工程α24β17**。不能沿用泛化“MEL5”名称合并天然与工程受体。完整变量链ID反映沉积构建体，不推断供体克隆或临床构建体身份。

对13个原始PDB坐标逐项重建并运行现有固定策略，只有6EQB、7N1F各1图通过。7RTR的缺失已直接检查ATOM记录，不是altloc或解析遗漏。4JFE的筛查触发点为肽C1 N与MHC A7 OH间距2.231154 Å、估算重叠0.838846 Å；这是保守粗筛的短供受体接触判据，**不能称作物理不可能，也不是MolProbity结论**。没有调整阈值、补原子或松弛结构使其通过。

来源与QC的联合状态见[cross_study_pair_readiness.json](cross_study_pair_readiness.json)，逐原子结果见[native_atomic_qc.json](native_atomic_qc.json)。2张候选图保存在native_graph_candidates.pt，明确trainable=false、无结合标签；不是训练release。

## 同V基因组合的不同受体对照

| V基因组合 | 完整变量链不同的候选PDB | 尚缺什么 |
| --- | --- | --- |
| human TRAV12-2 / TRBV6-5 | 5E9D、6RPB | 两个研究来源，但该组合中尚无跨研究同完整受体正对照 |
| human TRAV38-2/DV8 / TRBV25-1 | 5JZI、5YXU | 5YXU主引文DOI/PMID未核实；尚无跨研究同完整受体正对照 |
| mouse TRAV21/DV12 / TRBV29 | 7N5C、7N5P | 同一研究；尚无跨研究同完整受体正对照 |

这里的“负对照”只指身份检索中的不同受体，不是实验不结合。完整变量链不同本身也不证明来源、谱系或功能独立。

来源平衡要求：对研究A中的查询，同完整受体正例与同V组合的不同受体对照必须**共同来自另一研究B**。不能把B的正例和C的负例混合，让研究风格直接预测身份。当前3组跨论文正对照与3种同V对照的组合没有交集，合格查询块为0。排除图联合连接完整双链身份、共享单链、V组合与研究来源；22个连接分量中合格平衡分量为0。最少3个分量仅是“1个拟合、2个开发验证分量”的可计算门槛，不是统计功效或生物学放行标准。门槛与阻断均见[source_balanced_design.json](source_balanced_design.json)。

## 隔离项复核和仍未知事项

原285个隔离项的源文件SHA256全部重新核对一致。新增原始注释解决**1项8i5c_8（索引138）**的完整双链身份对应：同一聚合物实体的完整构建体一致，另一晶体副本提供编号证明。它仅为同PDB重复副本，不新增受体、研究或独立样本，未自动进入训练。其余284项（109天然、175预测）继续隔离，逐项原因见[isolation_followup.json](isolation_followup.json)。旧发布没有覆写。

8QFY/a42b20的主文确认其工程受体背景，但完整N端工程边界与变量域证据尚未核实；沉积序列以MAKEVEQNSG开头，不能改写为胚系QKEVE。DDD/EEE链不采用旧PDB单字母猜测。附件请求结果为403或HTML，**没有取得真实可核读的SI PDF，不知道其序列表内容**；原始请求见[8qfy_attachment_access.json](8qfy_attachment_access.json)。其HLA-E与YLQ的HLA-A2还构成MHC混淆，不能直接补成困难对照。

还需补足：完整原子且跨研究的同受体结构、同一gallery研究中的同V不同受体、更多互不连接的来源/家族分量，以及未知工程构建体的原始序列证明。42F3若要纳入紧凑MHC任务，需要另行定义并验证统一图范围；本轮没有为单个候选放宽范围。天然结构通过几何QC也不能验证任何合成诱饵为真实阴性。未解决项是当前证据不足，不代表永久无法获得。

## 验证与复现

19项相关测试全部通过，无跳过；候选身份、来源门、原285项复核JSON逐字重演一致，2图数组和记录重演一致。前3次冻结交付的211、233、592文件哈希均未改变。见[tests.log](tests.log)、[verification.json](verification.json)、[verification_final.log](verification_final.log)。失败请求、初版解析及修正产物一并保留，不选择有利版本隐藏问题。

复用现有环境，无安装。仓库根目录执行（请使用新的输出路径；身份解析和联合验证脚本拒绝覆写已有目标）：

```bash
export PYTHONPATH=src:.
PY=/home/administrator/docs/meta-acdc/.venv/bin/python
$PY -m unittest tests.test_cross_source_retrieval tests.test_paired_variable_identity -v
$PY scripts/curate_cross_source_receptors.py --out /tmp/new-cross-curation.json
$PY scripts/qc_cross_source_structures.py --out-dir /tmp/new-cross-native-qc
$PY scripts/verify_cross_source_delivery.py --replay-dir /tmp/new-cross-replay
```

原始候选元数据、编号快照、13个坐标和原始URL/哈希已入库，候选与原子QC可离线重演。原285项来源复核还依赖既有本地历史源文件归档；不能宣称新克隆无需历史归档即可重演这一部分。全文XML缓存保留本地；仓库提供原文定位、URL、哈希和事实摘要，不将无法取得的全文或附件当作已核读证据。

本轮仅完成数据取证、解析修正和放行条件核验。没有重洗划分、读取独立验收集、挑种子或运行新训练；结果不足以支持结构指纹增益、真实结合预测或临床风险结论。
