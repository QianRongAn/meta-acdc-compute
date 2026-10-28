# NP链身份、参考编号及天然配对裁定（2026-10-29）

本轮补齐NP105的公开同细胞天然配对来源，并把生物学α/δ身份与HMM分类区分开；核心四条β链的编号差异已定位并由IMGT参考裁定。V/等位基因、完整天然变量区及固定外部编号门仍有未解决项，未升级身份或冻结协议。

## 已解决：表达链身份不由最高V型HMM决定

[IMGT鼠TRAV等位基因表](https://imgt.org/IMGTrepertoire/Proteins/taballeles/mouse/TRA/TRAV/Mu_TRAVall.html)明确记录TRAV/DV基因可重排到TRAJ或TRD的(D)J基因。因此，共享V/DV特征和mouse_D略高的HMM分数不能决定成熟表达链是δ。

对最终版结构构建体的配对CDR3，回查公开GEO单细胞注释：核心NP202/NP101/NP207/NP203的α对应TRA、TRAJ53、TRAC，并在同一细胞配对TRB/TRBJ2-7；原文构建体也明确标α/β。证据支持这些是αβ受体的α链，经典工具的D型HMM首位仅是保留的算法输出。**生物学链角色问题已澄清，原始HMM输出没有被改写为A。**

NP1、NP2来自同一10鼠混池的不同条码化处理，不是两个独立来源。其核心四组β常数基因列为None：本轮没有补写TRBC。报告区分“TRAJ/TRAC与配对TRB/TRBJ支持链角色”和“双方常数基因均记录”的更严格证据层；后者对这些核心组仍为0。它们都不等于完整域身份准入。

## 已解决：所谓插入差异实际是FR1缺口位置

四条核心β链只有三个沉积残基的标签不同：

| 沉积索引（0起点） | 氨基酸 | ANARCII | 经典ANARCI | IMGT TRBV5带缺口参考 |
|---|---|---|---|
| 13 | G | 14 | 15 | 15 |
| 14 | R | 15 | 16 | 16 |
| 15 | A | 16 | 17 | 17 |

两工具把同一个FR1缺口分别放在17与14；没有新增或删除氨基酸，也没有CDR3插入分歧。TRBV5*01/*02的官方带缺口参考支持缺口14，故这些争议位置支持经典ANARCI。其余位点逐一相同。[逐位裁定及原始参考头](adjudication_source_tiers.json)。

参考坐标使用了保护条件：只对首/次保守Cys分别落在23/104且完整框架精确匹配的TRBV参考直接使用整数槽位。TRAV参考可带额外对齐列，不能把每个原始三联体槽位直接当作IMGT整数编号；本轮没有用该捷径强行“统一”α全域编号。

## V/J归属：从首尾10残基缩小到完整框架与较长J片段

下表以整个V框架至保守Cys的连续氨基酸完全一致筛选，J侧比较完整参考的最长精确后缀。它们是证据相容集合，不是唯一基因组定位或历史重排切点。

| 链组 | 完整V框架相容等位基因 | 最长J后缀支持 |
|---|---|---|
| 核心四组α | TRAV16*06, TRAV16D/DV11*03 | TRAJ53*01 (F, 19 aa) |
| 核心四组β | TRBV5*01, TRBV5*02 | TRBJ2-7*01 (F, 14 aa); TRBJ2-7*02 ((F), 14 aa) |
| NP204_alpha | TRAV7-4*02, TRAV7D-4*01, TRAV7D-4*02, TRAV7N-4*01 | TRAJ56*01 (F, 20 aa) |
| NP204_beta | TRBV1*01 | TRBJ1-6*01 (ORF, 17 aa) |
| NP105_alpha | TRAV7-3*04, TRAV7D-3*01, TRAV7D-3*03 | TRAJ23*01 (F, 19 aa) |
| NP105_beta | TRBV2*01 | TRBJ1-6*01 (ORF, 17 aa) |

核心α框架无法用蛋白序列区分TRAV16*06与TRAV16D/DV11*03；核心β无法区分TRBV5*01/*02。较长J片段将核心β归到与原文一致的TRBJ2-7相容集合，超出了原先J末10残基的多解，但仍不能区分同蛋白等位基因或反推出实际核酸切点。

NP204/NP105 β的17 aa后缀精确支持原文TRBJ1-6*01，而非经典工具约0.692相似度的TRBJ2-4近邻调用。原始IMGT参考将TRBJ1-6标为ORF，单细胞productive标记属于另一个判据，不能把该参考改成F。**既定功能性J准入门保持。**

## NP105天然配对来源：已补到92个严格同细胞条码

前一轮最终Data S2的242克隆主表没有NP105完整配对，是因为NP105来自另一批再挑战样本。[GEO GSE302680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302680)的设计和[GSM9108725](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9108725)记录DN4_LN_TCR，来自2只B6鼠混池、X31异源再挑战后第4天。其公开clonotypes和filtered_contig_annotations文件已下载并交叉核验。

98个条码同时出现目标双方CDR3；92个通过严格检查：is_cell、high_confidence、full_length、productive均为True，只有一条活跃TRA和一条TRB，同一细胞、同一克隆号，CDR3核酸翻译一致，V/J与原文相符，常数调用分别为TRAC和TRBC1，克隆表也同时记录该双链。

- α：TRAV7D-3 / TRAJ23，CDR3 `CAVMNYNQGKLIF`。
- β：TRBV2 / TRBJ1-6，CDR3 `CASSQDRRSSYNSPLYF`。
- 原始克隆号clonotype1。92个条码支持一个天然配对受体的来源，不是92个独立受体或个体。6个未达严格层的候选不混入92；另一个仅缺β常数调用的条码只计入链角色层。

[同细胞配对与逐条注释收据](adjudication_source_tiers.json)保留严格条码、双方contig ID、CDR3核酸、V/J/C调用及拒收计数。NP1/NP2没有NP105双链命中，不把其他α链共享的β命中认作NP105。

## 仍未解决：完整原生变量区与外部完整编号

公开注释CSV只有18列，不含完整contig核酸或框架序列。full_length=True是来源QC标记，不是公开了完整序列。因此配对来源已补齐，但不能用V/J标签与CDR3拼接成“已观测完整天然变量链”，也不能证明工程结构构建体每个变量区残基均与天然细胞完全相同。

[ENA原始读段收据](sources/read_metadata/NP105_ENA_reads.tsv)进一步定位到SRX29686421 / SRR34543498，两个FASTQ压缩文件共6,277,171,658字节。原始reads有公开入口；本轮未下载并重建完整V(D)J，不能声称源数据不存在。后续完整序列结论需要从这些读段重建并验证细胞、UMI、配对和组装，不能以参考补齐代替。

α的生物学角色与β的FR1争议已澄清，但外部STCRDab完整配对编号仍缺，部分等位基因仍有不可由现有蛋白框架区分的解，TRBJ1-6的固定功能性参考门仍不满足。没有修改既有角色、边界或原子门代码，没有自动准入NP结构或释放训练数据。

## 跨研究可用正对照：重取坐标后仍未放行

本轮针对上一轮377个命中中最接近可用的两组完整身份，重新下载5个PDB的当前原始坐标，并用原固定原子与完整class-I单元规则重跑。不是以近邻命中率代替泄漏对照。

| PDB | 当前固定规则结果 |
|---|---|
| 7N1F | 通过；图数组哈希与既有库存一致 |
| 7RTR | E28仍缺CD/CG/OE1/OE2界面重原子 |
| 6EQB | 通过；图数组哈希与既有库存一致 |
| 4JFE | 仍被固定原子筛查拒收 |
| 4JFH | 无要求的完整pMHC单元 |

7N1F–7RTR、6EQB–4JFE、6EQB–4JFH三对均未达到双端通过。两张通过图是已有开发库存的精确重演，不增加受体数。固定原子粗筛的拒收不等于物理无效；未补原子或调阈值。[新鲜坐标复核](fresh_cross_summary.json)、[完整原子报告](fresh_cross_QC/native_atomic_qc.json)。

当前身份库存仍64 PDB/44配对变量ID，固定QC 43图/31 PDB/29 ID；合格跨研究同完整αβ正对照、合格同来源同V独立天然结构对照组、来源平衡块均为0。协议未冻结、未训练，独立验收集未访问。

## 复验与文件

67项测试通过；十一份既有冻结交付4474个文件哈希保持。天然配对/参考裁定逐字节重演，新鲜坐标QC与两张图数组精确重演；原始门代码保持。[验证](verification_final.json)、[测试日志](tests_final.log)、[交付清单](DELIVERY_MANIFEST.json)。

```bash
PYTHONPATH=src:. python scripts/verify_np_adjudication.py --replay-dir /tmp/np-adjudication-replay-unique
```

原始GEO/IMGT网页、CSV.gz、坐标和HTTP失败均留收据。PDB按gzip保存并保留未压缩SHA256；验证会在新目录恢复，不重写旧结果。下载失败不当作阴性；[归档布局](archive_layout.json)列出原始长度、哈希与包装。两个PT图明确trainable=false、无结合标签，是公开开发数据复验产物。
