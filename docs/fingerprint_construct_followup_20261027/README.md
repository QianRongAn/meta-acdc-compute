# 6AT6、7N5C、7N5P原文链身份与天然困难对照续查

2026-10-27。**三项指定链角色冲突已用原文序列证据裁定；补到AS3.1／AS8.4不同患者天然克隆的同来源同V依据，但尚不能组成合格结构对。合格跨研究同完整αβ正对照仍为0，来源平衡块仍为0；不冻结协议、不生成训练release、不训练，独立验收集未访问。** 根README同步更新；历史冻结交付保持原字节。

## 指定条目的原文构建体与链身份

| PDB | 证据支持的身份与作者链 | 构建体核验 | 当前准入 |
|---|---|---|---|
| 6AT6 | KFJ5；A为α、B为β；TRAV4/TRBV28 | 原文Table 1的六条CDR均在唯一对应沉积链中匹配；方法记载恒定区工程二硫键。不能将沉积SP3.4等模板名当作变量链身份 | 完整配对变量链身份恢复；该结构未结合pMHC，不能补class-I复合物正对照 |
| 7N5C | 6218αC；D为α、E为β；鼠变量区、人恒定区 | 原文Fig.1a、合成方法及SI Table 1（PDF第12页）共同支持；α链Ser110→Cys，工程替换完整保留 | 通过既有完整性和原子QC，新增恢复1张候选图；不是新独立天然受体 |
| 7N5P | 6218亲本；D为α、E为β；鼠变量区、人恒定区 | Fig.1a及方法、结构编号对应均支持；与7N5C的β变量链相同，α变量链只差论文指定的S/C | 身份恢复，但现有原子规则拒收；不作为合格结构对照 |

2018原文：[Chan等，Divergent T-cell receptor recognition modes of a HLA-I restricted extended tumour-associated peptide](https://doi.org/10.1038/s41467-018-03321-w)。2022原文：[Covalent TCR-peptide-MHC interactions induce T cell activation and redirect T cell fate in the thymus](https://doi.org/10.1038/s41467-022-32692-4)。两篇CC BY原始XML、SI及关键图像快照与请求回执在sources/中保留。

本轮没有把通用解析器改为“优先相信molecule、忽略synonym”。原始冲突保持记录；只对有明确论文序列锚点的指定条目建立附加裁定。逐链SEQRES必须与原API构建体完全一致；完整变量域仍重新运行原有STCRDab编号、IMGT功能性V首端及J/FR4末端精确验证。[最终裁定与全部序列](construct_adjudication_final.json)保留作者链、原声明、文献定位、序列及源哈希。

## 同来源同V的独立天然受体：找到来源依据，但结构不合格

[Yang等，Autoimmunity-associated T cell receptors recognize HLA-B*27-bound peptides](https://doi.org/10.1038/s41586-022-05501-7)的SI Table 2（PDF第4页）提供：

- 患者P的AS3.1：TRAV21/TRBV9；α CDR3为`CAVSLGTGAGSYQLTF`，β为`CASSVGLYSTDTQYF`。
- 患者K的AS8.4：同V组合；α CDR3为`CAVNSPGSGAGSYQLTF`，β为`CASSVGTYSTDTQYF`。

这是两个明确不同患者的单细胞αβ配对，不是由名称或相似度推测的工程前后版本。SI Table 6（PDF第8页）将其对应到AS3.1的7N2S和AS8.4的8CX4。[事实转录、原文定位及SI哈希](natural_control_primary_receipt.json)保留核验依据。

**仍不能将它们纳入合格结构对照：** 7N2S实际沉积β链起始为`GVTQTP...`，没有通过原有完整N端规则；同论文7N2N/O/P/Q/R也有真实沉积β端部缺口，其中N/O还有α FR4缺口。没有从原文全序列补字母。8CX4则已用SI Table 2的两条CDR3唯一匹配作者链D（α）与F（β），并重新验证完整变量链，但完整坐标触发现有原子筛查，仍不准入图。

原文天然克隆来源证明与合格结构证明是两层证据。当前准入后的同来源同V组仍只有JDI亲本与工程变体一组，不能称独立天然克隆或实验不结合对照。6218亲本／S110C也不能凑成独立天然组。

## 跨研究同完整αβ正对照补查

本轮对KFJ5（以6AVF完整变量链为查询）、6218、6218αC、YLQ7和6VRM逐条α／β进行RCSB 100%比对身份检索：10条查询，2条首次500后重试成功，最终10/10有可读结果。原始参数和返回均归档；**100%局部比对命中只用于发现，必须再核对全长与链配对**。

| 线索 | 本轮结果与不准入原因 |
|---|---|
| KFJ5 | 双链命中6AT6、6AVF；同一2018论文，前者未结合pMHC，不能作为跨研究正对照 |
| 6218亲本 | 双链检索交集含3PQY、7N4K、7N5P。2011年[Day等原文](https://doi.org/10.1073/pnas.1106851108)确认6218谱系；但3PQY及7N4K实际沉积α链从`KTTQP...`起始，不含当前完整变量链起始`DAKTTQP...`中的`DA`。β完整匹配不弥补α缺口；不补DA、不用同克隆名称当完整沉积变量域相等 |
| 6218αC | α链只命中7N5C；与亲本保留工程差异，不合并身份 |
| YLQ7 | 双链命中7N1D、7N1F、7RTR；7N1D未结合pMHC，7RTR既有界面侧链缺原子问题未解决；不构成QC通过的跨研究对 |
| 6VRM | 双链仍命中6VRM、同研究未结合pMHC的6VTH；没有新增合格跨研究来源 |

[最终联合报告](joint_readiness_final.json)包含各查询原始命中、PDB交集、6218实际端部比较、完整资格判定与可追溯哈希。检索范围是上述锚点及已有候选，不是对所有公共结构已作穷尽证明。

## 原子检查的具体阻断

沿用相同图范围、完整重原子要求及固定原子policy，没有改坐标、补原子或调阈值：

- 7N5P：肽C1的OG与MHC A163的OE1距离2.391434 Å，估算重叠0.648566 Å，触发现有短供受体接触筛查。
- 8CX4：肽C1的O与MHC A159的OH距离2.373012 Å，估算重叠0.666988 Å，触发同一筛查。

这两项是保守工程准入规则的拒收，**不是MolProbity结论，也不宣称原晶体物理无效**。本轮没有为使条目通过而修改规则。[最终逐原子结果](native_qc_final/native_atomic_qc.json)同时保留通过与拒收记录；候选PT明确trainable=false、无结合标签。

## 当前计数与冻结条件

| 层级 | 当前数量 |
|---|---:|
| 身份支持的PDB／完整配对变量序列ID | 64／44 |
| 固定原子QC通过图／PDB／配对变量ID | 43／31／29 |
| QC通过的跨研究同完整αβ正对照 | 0 |
| QC通过的同来源同V不同序列组 | 1（JDI工程谱系） |
| QC通过的独立天然受体对照组 | 0 |
| 来源平衡块／作者不重叠平衡块／可用联合分量 | 0／0／0 |

相对上一批60个身份PDB，恢复三项指定条目并补足8CX4角色证据，得到64；新增QC图只有7N5C一项，42变43。多个晶体副本／变量序列键不是独立样本量。

冻结前仍需同一完整配对变量身份的合格跨研究正伙伴，以及共同gallery研究中的同V独立天然受体对照；训练、验证还须满足既有联合来源／家族／单链隔离。最低三个可用联合分量只表示可计算，不保证统计功效。当前条件未满足，协议保持未冻结，独立SPR／功能／临床验收集未访问。

## 验证、原始资料与复现

[最终验收](verification_final.json)：40项相关测试，无跳过；四项角色裁定、原子JSON及联合报告逐字节重放一致；7N5C候选图数组与元数据重建一致。七份既有冻结交付3622个文件逐哈希未变，通用角色、变量域、图和原子／来源平衡规则代码哈希未变。[测试日志](tests.log)、replay_transport/replay.log及DELIVERY_MANIFEST.json记录检查范围。

复用既有环境，不安装依赖。仓库根目录运行，输出目录须不存在：

```sh
PYTHONPATH=src:. /home/administrator/docs/meta-acdc/.venv/bin/python scripts/verify_construct_followup.py --replay-dir /tmp/construct-followup-replay-UNIQUE
```

2018／2022两篇CC BY资料保留全文与归属；2011 PNAS与AS研究的未确认CC再分发许可全文／SI保留本地下载，Git仅含事实摘要、序列转录、源链接与SHA256。AS事实转录是人工核读SI Table 2/6后的科学事实，不声称已公开其完整附件；新克隆可重演计算链，但重新核读这两份完整原文需由原URL获取或使用本地归档。见[SOURCE_ATTRIBUTION.md](SOURCE_ATTRIBUTION.md)。

失败请求、初版输出均保留。`construct_adjudication_final.json`、`native_qc_final/`、`joint_readiness_final.json`和`verification_final.json`为本轮最终版本。初版8CX4曾把旧头文件用作原子检查输入，得到空单元；已重新下载完整坐标并加“无蛋白坐标即报错”检查，最终按实际原子筛查拒收。初版空单元不作为生物学结论；[SUPERSEDED.json](SUPERSEDED.json)列明被替代文件。未覆写任何旧冻结结果。


## CC BY原始附件的传输包装

一次性推送连接重置，两个CC BY补充PDF改以2MB以内的原字节分段入库；本地原始PDF和首次完整提交另存保留。不是重新压缩或编辑附件，原始请求回执与SHA256不变。

```sh
python scripts/reassemble_construct_supplements.py --name 6at6 --out /tmp/6at6-original.pdf
python scripts/reassemble_construct_supplements.py --name 7n5 --out /tmp/7n5-original.pdf
```

[分段清单](sources/supplement_parts/manifest.json)绑定各段和两个原PDF的长度与SHA256。重组脚本拒绝覆盖已有输出，最终验收会重新组装两份并逐字节校验；新克隆可核读相同CC BY附件。该调整不改变角色裁定、原子规则、候选图或联合放行结论。Git按小批次发布，每批次仅是传输进度，完整结果以最终树为准。
