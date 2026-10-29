# NP105整run复验与TCR1414来源核验（2026-10-30）

原计划现已推进到**完整原始run的序列重建和逐细胞复核**。SRR34543498全文件通过官方大小/MD5验证；全部50,653,512个spot按原筛选协议处理，获得1,849,447对reads、92个预先登记条码。60个细胞组装出完整αβ编号区间，54个具有唯一一致的候选核苷酸双链且每一碱基有同细胞Q30支持。**这是同一个NP105来源的序列证据，不是54个独立受体，也不是模型阳性。**

新增合格跨研究同完整αβ正对照、同来源同V独立天然结构对照组均为0；固定QC仍为43图／31 PDB／29 ID。不降低标准、不冻结协议、不启动训练，独立验收集未访问。旧的“整run尚未完成”描述只适用于前两轮历史快照。

## 完整数据与复验范围

- [下载和官方MD5](whole_sra_download.json)、[存储恢复后再次校验](post_storage_recovery_sra_check.json)：6,040,866,339字节，MD5 `bb8d681b918233fada8e0dbe875e4887`，SHA256 `b7dc3b09865eb381d102f7be3c7ff3e129a0297897e42cfc98e8ca4935628bbb`。
- [完整筛选收据](whole_selection_receipt.json)：92个原白名单条码不变；原先546,555个spot处全部186个筛选输出的字节数和SHA256精确复现。全量186个文件也在磁盘上逐个复核，[结果](whole_selection_inputs_verified.json)。SRA最初3个spot的6条生物read及逐碱基质量与ENA前缀一致，[跨格式核对](first3_cross_format.json)。8nt技术read未冒充生物read。
- 完全沿用原始barcode/UMI/Q20处理；无barcode纠错、无参考补缺、无生物碱基替换。每个细胞使用RNA SPAdes 3.15.5、2线程、4GB、`--only-assembler`，未使用参考序列或read纠错。[工具版本/包哈希](sra_tool_manifest.json)和[编号依赖恢复核对](restored_numbering_dependencies.json)。
- [输入范围和失败项](whole_run_scope.json)：92次退出码均为0；88个有非空transcripts.fasta，4个没有转录本文件。四个无输出细胞的日志均记录插入长度估计失败、无正确成对比对及跳过repeat resolution；不能从正常退出推出成功重建，也不据此断定具体生物学原因。原始日志完整归档。
- 未改动旧驱动，因此`whole_assemblies.json`中继承的`same bounded prefix`文字是过时描述。实际输入范围由完整筛选收据、186个输入哈希、逐细胞命令及本轮scope文件界定；原始报告逐字保留，没有事后改写。

## 逐细胞结果

| 项目 | 本轮完整run |
|---|---:|
| 原白名单细胞 | 92 |
| 非空转录本输出 | 88 |
| CDR3-bearing ORF / 完整编号ORF | 443 / 394 |
| αβ均完整 / 仅β完整 / 仅α完整 / 两者均不完整 | 60 / 17 / 7 / 8 |
| 60个中，旧候选每一碱基均有同细胞Q30支持 | 57 |
| 同时满足每条链唯一完整核苷酸序列、与旧观察候选一致、每碱基Q30 | 54 |
| 上述54个中，全部完整ORF两种编号均同角色且逐位一致 | 26 |
| 其余54个中的编号失败或分歧细胞 | 28 |

[全部92细胞摘要](cell_summary.json)。Q30以BWA primary、MAPQ≥30、baseQ≥30、同barcode+UMI去重；某一位同UMI有冲突即不计其无歧义支持。比对使用前轮**观察到的组装ORF**，只用于验证；没有把比对模板输入组装或拿它替换碱基。这一覆盖检查仍可能受比对偏倚、barcode/UMI碰撞或扩增影响，不能独立证明天然基因组单倍型。

3个质量缺口细胞：CCTATTAGTCAATGTC、CGCCAAGGTCACTGGC、GAACCTACACCTGGTG。另3个存在替代核苷酸组装：ACTTACTAGTGTCCAT、GCTGCAGTCTAACGGT、TAAACCGAGGGAAACA；不因同时能比对到旧候选而把这些组装自动改回。它们与质量缺口集合不重叠。

前缀37个完整双链细胞到全量后的状态为23个仍双链、9个仅β、5个两链均不完整；[完整状态转移及变体](whole_variant_inventory.json)。深度增加不保证组装结果单调改善；不跨深度择优合并，不将旧失败项删除。

原前缀阶段预先按read量选定的3个细胞用全量输入再次运行，均与各自首轮全量transcripts.fasta**逐字节相同**，[收据](whole_predeclared3_replay.json)。其中包括有替代组装的TAAACCGAGGGAAACA。只重复这3个；没有宣称92个都重复运行。

## 序列、域编号与V/J：哪些已补齐，哪些没有

[观察核苷酸FASTA](observed_variable.nt.fasta)和[蛋白FASTA](observed_variable.aa.fasta)直接导出自全量独立组装、Q30及两方法均通过的ACAGCCGAGAAAGTGG；[转录本和哈希出处](observed_variable_provenance.json)。α编号观察区间333nt/111aa，β348nt/116aa。两者等于前轮三个预选细胞观察序列，未用IMGT补碱基。

固定IMGT参考快照下，天然核苷酸框架精确匹配支持α **TRAV7D-3*03 / TRAJ23*01**、β **TRBV2*01 / TRBJ1-6*01**；蛋白层α V仍有TRAV7-3*04、TRAV7D-3*01/*03多解。β J全53nt匹配但数据库功能标记仍是**ORF**，不能改写成F；原HMM近似J首位TRBJ2-4也保留。详细参考、原始标签及坐标见[前轮核苷酸取证](../fingerprint_np105_raw_20261029/README.md)。本轮只增加原始reads支持，不更换参考快照或重定义准入。

V精确匹配区间为α[0,272)、β[0,287)，完整J片段为α[274,334)、β[296,349)，均为0起始半开区间。这些是参考匹配区间，**不是唯一裁定的V/J切除、D片段或N插入边界**；基因组相位仍缺。α原构建体112aa在编号128P之后还有观察到的N，β显式128缺口均按原样保留；FASTA只覆盖上述编号区间。本轮54个的逐碱基检查不等于54个均已验证额外N密码子，不能外推前轮3个细胞的末端检查。

ANARCI与ANARCII对394个classic完整ORF有339个同角色、相同非gap位置映射，**55个分歧继续隔离**：[逐项列表及细胞交集](whole_numbering_disagreement_summary.json)。其中45个神经编号失败（44个分数不足、1个非插入残基不足），5个角色不符，5个同角色但映射不同。保留完整ORF输入和原始输出，不裁剪成更易通过的输入，不强改α/δ/β标签，也不推断失败机制已查明。26个联合一致细胞也不能绕过β J ORF及固定结构准入门。

## 替代组装的原始read复查

全量4个替代ORF对应3组barcode/链核苷酸候选；再加原先隔离的K82N，逐一检查所有差异位点的Q30同UMI相位、同read相位及精确CDR3关联。[原始支持计数](whole_variant_raw_support.json)：

| 候选 | 差异位置（编号区间nt，1起始） | 候选无歧义UMI覆盖全部差异 | 同时关联精确CDR3的UMI |
|---|---|---:|---:|
| ACTTACTAGTGTCCAT TRA | 39,132 | 0 | 0 |
| GCTGCAGTCTAACGGT TRB | 39 | 0 | 0 |
| TAAACCGAGGGAAACA TRA | 4,12,55,56,60,71 | 1 | 0 |
| TCGGGACTCGTGGGAA TRA（旧K82N） | 201 | 0 | 0 |

不能据此确认独立天然受体，也不把它们定性为已证实的测序/PCR错误。候选、低质量和冲突证据保留，不修复、不新增身份。[重构辅助函数后整SAM重查一致](variant_refactor_recheck.json)。

## 跨研究与天然对照续查

接续前轮41个本周RCSB公开命中及282目标搜索，原3条失败请求已全部成功；TCR1414的4条完整链查询现在均有成功回执。逐条精确子串验证后，双链交集仍仅TCR1414_1的9X87/9XC3和TCR1414_2的9X9I。[原始请求](retry_requests.json)、[来源裁定](tcr1414_source_audit.json)。这轮补的是这几条未决线索，未声称重跑了整个数据库。

找到了[2023受体生成原文](https://doi.org/10.1136/jitc-2022-006001)：其方法分别取得α/β候选链并测试15个组合，筛出TCR1414_1/_2；这没有建立同一细胞的天然αβ配对。图3的V/J/CDR3与结构链名相容，但不能替代跨论文完整变量序列证据。9X87是无pMHC的未结合受体，9XC3/9X9I的结构论文DOI仍未补到；不能将受体生成论文DOI填作结构论文DOI来制造跨研究正对照。两个受体的V组合也不同，不能用作同V天然对照。新增合格数继续为0，原因逐项记录。

原文和图像因CC BY-NC许可保留本地；版本化内容为事实摘录、请求、SHA256和段落定位哈希，未把完整原文放进Git。

## 可复核交付与复现

[evidence](evidence)保存四份原始JSON的逐字gzip，以及92个组装任务的原始transcripts/params/log/input YAML归档（88份transcripts，4份缺失不伪造空FASTA）。[归档布局](archive_layout.json)逐文件记录长度和SHA256；原文件不删除。6GB SRA、186个FASTQ、SAM与大型临时目录留在外部存储，不推送Git。全量输入文件名、字节数和SHA均在筛选收据内。

轻量复核（解压核对归档、重算主要计数、逐文件复核14份旧交付4625个文件，并运行99项测试）：

```sh
python scripts/verify_np105_recovery.py
```

[验证结果](verification_final.json)、[99项测试日志](tests_final.log)、[既有交付/20个原代码/5个固定门未变](frozen_integrity_after_storage_recovery.json)。轻量复核不声称重新下载整run或重复92次组装。

完整重演需先按下载收据取得MD5一致的完整SRA，以及固定版本SRA Toolkit 3.2.1、SPAdes 3.15.5、BWA 0.7.19、ANARCI 2026.2.13.2、ANARCII 2.0.8及各自依赖。当前驱动保留本机固定工具位置；安装包和恢复映射见工具收据与脚本，不静默改用最新版。以下示例输出目录必须尚不存在，D盘需足够空间，`PYTHONPATH=src:.`：

```sh
python scripts/select_np105_full_sra.py --sra /path/to/verified/SRR34543498.sra --out /mnt/d/np105-replay/selected
python scripts/assemble_np105_all_cells.py --selected /mnt/d/np105-replay/selected --out /mnt/d/np105-replay/assemblies --report /mnt/d/np105-replay/assemblies.json
python scripts/number_np105_all_cells.py --assemblies /mnt/d/np105-replay/assemblies.json --out /mnt/d/np105-replay/classic.json
python scripts/validate_np105_whole_pairs.py --selected /mnt/d/np105-replay/selected --numbering /mnt/d/np105-replay/classic.json --out /mnt/d/np105-replay/validation
python scripts/number_np105_whole_neural.py --classic /mnt/d/np105-replay/classic.json --out /mnt/d/np105-replay/neural.json
python scripts/audit_np105_whole_variants.py --numbering /mnt/d/np105-replay/classic.json --sam /mnt/d/np105-replay/validation/whole_reads.sam --out /mnt/d/np105-replay/variants.json
```

存储中断后重新核验完整SRA、既有交付与编号依赖。新的大型筛选/临时/组装/比对文件均在D盘，导出增加D盘路径与空间检查。旧579,862,528字节部分下载及原临时文件保留；未删除用户数据，未移动WSL发行版。此前错误地只检查WSL虚拟空间、漏检C盘物理空间的问题已在本轮操作中纠正。

## 尚未解决与下一门槛

1. NP105已有完整观察编号区间与多细胞Q30证据，但NN分歧、β J ORF、额外末端及唯一生物学重组边界仍未全部闭环；原结构完整身份门不改。
2. 替代组装缺乏无歧义且精确CDR3关联的候选UMI支持，保留隔离，不用参考纠正。
3. TCR1414缺跨论文完整变量序列/天然同细胞配对/结构原文来源；已知重复PDB不能补齐这些证据。
4. 合格跨研究同完整αβ结构正对照、同来源同V独立天然结构对照仍缺。取得原文完整序列并通过原固定QC后，才允许讨论冻结及后续身份修正评估。
5. 本轮没有运行模型或获得新的性能证据；现有模型尚未显示超过残基组成基线的受体信息。序列恢复成功不能改变这一结论。
