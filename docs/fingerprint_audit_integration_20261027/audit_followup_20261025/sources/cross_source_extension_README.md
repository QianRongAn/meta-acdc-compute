# 跨来源结构取证扩展（2026-10-24）

**本轮补上1组原子QC通过的同研究、同V工程构建体对照，但仍没有原子QC通过的跨研究同完整αβ变量链正对照；来源平衡块仍为0。未冻结新协议、未生成训练release、未训练新模型。** 独立SPR/功能/临床验收集未访问，未重洗划分或放宽准入。

## 新增搜索与准入

复用此前799-PDB TCR3d/RCSB快照，先通过功能性IMGT V/J端部匹配发现候选，再要求STCRDab显式配对及编号和原始沉积构建体完全一致。发现阶段的相似性、名称、V基因多解不能充当标签。随后覆盖快照中**全部307个CLASSI PDB**，连同游离/溯源线索共331个PDB下载编号和FASTA；联合旧176项形成扩展证据。45个序列层候选身份的两条链均已有RCSB精确搜索结果（本轮补42条链查询），未发现该799条快照之外的新双链交集PDB。精确局部比对不等于完整变量链，搜索命中仍要重新核验。

严格原始头文件核验后：113个链对单元、61个PDB、41个完整变量链ID通过身份层；这些不是41个相互独立生物学受体。新加的primary-header拒收项含5E6I、6VRM、6VRN、8CX4，COMPND身份角色未得到充分一致支持，未按序列相似猜标签。缺编号、端部证据不足和角色冲突均逐项留在[sequence_curation_all_class_I.json](sequence_curation_all_class_I.json)、[source_balanced_design.json](source_balanced_design.json)。

对61个身份层通过PDB的原始坐标全部运行既有完整class-I原子规则：**30个PDB、42个链对单元、28个完整变量链ID通过QC**。同PDB多个副本不增加研究/受体数；42不是独立样本量。未加原子、未优化坐标、未放宽阈值。见[native_atomic_qc.json](native_atomic_qc.json)。候选PT明确trainable=false、无结合标签，不是可训练release。

## 跨研究同受体：为什么仍不能放行

| 身份组 | 结果 |
| --- | --- |
| YLQ7 / YLQ-SG3，7N1F–7RTR | 作者不重叠；7N1F通过，7RTR E28缺CG/CD/OE1/OE2，原子规则仍拒收 |
| 工程α24β17，4JFE/4JFH–6EQB | 6EQB通过；4JFE固定粗筛拒收，4JFH游离无pMHC；作者有重叠 |
| 42F3，2011与2015两篇 | 完整变量链一致；紧凑MHC缺β2m，超出现有完整class-I范围；不是认定生物学无效 |
| 新发现XV19，4EI6–4ELK | 完整变量链一致、两篇论文作者不重叠；但为游离/CD1背景，不能补入class-I任务 |

新增XV19提升的是溯源覆盖，不是合格class-I正对照数量。4JFE保守粗筛不是物理金标准，没有据此声称原结构不可能。各组具体证据继承[上一轮报告](../fingerprint_cross_source_20261024/README.md)并在[joint_readiness.json](joint_readiness.json)联合重算。所有身份级跨论文组在原子QC后均没有完整正对照对。

## 已补齐：同一研究中的同V不同构建体

[JDI研究原文](https://doi.org/10.1038/s41467-022-32811-1)说明7PB2是亲本JDI，7OW5/7OW6是工程JDIa41b1。两种完整变量链ID不同、V组合均为TRAV19/TRBV6-2；3个PDB均有通过原子QC的单元。原文Table 4的五条CDR序列逐项与沉积变量链匹配，证明见[jdi_construct_receipt.json](jdi_construct_receipt.json)。

7OW6与7PB2还具有相同沉积pMHC构建体与G12D肽，可保留为同研究、同V、同pMHC的不同工程构建体对照。**它们属于同一个工程谱系，不能称不同独立天然克隆、不能算跨研究复现、更不是实验阴性。** 7PB2的另一个晶体副本被固定原子规则拒收，拒收记录也保留；不把副本选择包装为新增受体。

来源平衡仍要求研究A的查询同时检索研究B的同受体正例及同V不同受体对照。JDI没有跨研究同完整受体正例，不能拿另一研究的负例拼接凑块。当前QC后仅1个同研究同V对照组、0跨研究正对照对、0来源平衡块、0严格作者独立平衡块。最低分量门槛是计算可行性要求，不是样本功效保证。

## 隔离项与8QFY原文补证

旧285隔离项再次逐一核对源文件SHA一致，结果仍是上一轮仅8i5c_8的身份对应已补足、其余284隔离；本轮没有新增解隔离。该1项是同PDB副本，仍未自动准入训练。见[isolation_followup.json](isolation_followup.json)。

8QFY官方Europe PMC补充附件接口本次成功返回ZIP，保存原始包、提取23页SI PDF并渲染核对关键页。补证确认a42b20与8QFY的构建体对应及结晶/SPR方法，但没有补足完整工程序列/端部边界证据，继续隔离；不把分子置换搜索模板4WW2认作其受体身份。证据与页码见[8qfy_primary_supplement_receipt.json](8qfy_primary_supplement_receipt.json)。原文和SI由作者2024发表于PNAS，DOI[10.1073/pnas.2318003121](https://doi.org/10.1073/pnas.2318003121)，按[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)保留来源并归属作者；原始下载和失败请求全部保留。没有将该研究实测读数追加到训练或验收集。

## 旧v2外部复核回应与可交换元数据

已用精确逐图索引重演15个旧实例的原指标，恢复0.7943805及160旧标签/55多图标签；补出一个已核实完整变量链身份跨旧折实例。实际输入的晶体和预测置信度列均为0.5、归一化均为0，不能把假想AF3列变化套在这份数据上。详见[逐项回应](EXTERNAL_REVIEW_RESPONSE.md)。

提供[不含图的配套包](legacy_review_package.zip)：原30个sim/emb矩阵、331行精确索引映射、折分配、feature_schema、15实例原指标重演绑定及文件哈希。旧records缺MHC，明确标null。此包只供已消费历史面板诊断，不是干净训练数据或独立验收数据。

## 检查与重演

29项相关测试通过，无跳过。候选身份、来源门、42图数组、隔离项、15旧实例指标与逐图元数据均独立脚本重演一致；前4批冻结交付保持哈希不变。见[verification.json](verification.json)、[tests.log](tests.log)、[replay/replay.log](replay/replay.log)。重演复用现有环境，不安装依赖：

```bash
export PYTHONPATH=src:.
PY=/home/administrator/docs/meta-acdc/.venv/bin/python
$PY scripts/verify_cross_source_extension.py --replay-dir /tmp/new-cross-source-extension-replay
```

候选解析及原子QC从仓库快照离线重演；旧285源文件重核仍需要本地历史结构归档，不能宣称纯新克隆完全重演该部分。上述检查不解除任何生物学/临床证据阻断。

仍缺：固定规则下原子完整的跨研究同完整变量链受体、其共同gallery来源里的同V对照、足够独立来源/家族分量和缺失工程端部的原始证明。当前检索覆盖有限快照及指定精确查询，并非证明所有公共库永远没有这些证据。保留未知，不通过重新洗牌、补原子或改阈值制造可用性。

## 原始附件的分段传输

Git上传连接反复重置，因此18,093,497字节的官方原始ZIP以5个分段入库，原包及原始完整提交仍在本地保留；提取的SI PDF和其他证据不变。`sources/attachment_probe/multipart_manifest.json`记录每段及原包的长度和SHA256。在仓库根运行`python scripts/reassemble_cross_source_attachment.py --out /tmp/supplement-original.zip`可得到逐字节相同的原包；已有输出会拒绝覆盖。该调整仅为传输包装，不改变数据准入或实验结果。
