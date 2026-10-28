# NP最终版、域编号与更广跨研究检索（2026-10-29）

本轮取到最终版所需补充材料，公开逐残基编号候选与精确沉积区间，并扩大外部精确搜索。**尚未补齐可准入的NP完整域配对或跨研究正对照；标准、原子门和独立验收集保持不变。** 未冻结协议、未训练。

## 最终版附件：已补到

出版社API元数据确认DOI [10.1016/j.celrep.2026.117908](https://doi.org/10.1016/j.celrep.2026.117908)、PII及最终版题名。API只返回元数据，不是全文；FULL请求401，网页403，Europe PMC搜索503，NCBI转换接口报告未收录PMC。上述访问状态都保留，不能当作没有补充材料。

沿同一PII的出版社CDN取到22页补充PDF及6份工作簿。PDF封面已渲染核对，题名和作者相符；第18页为Table S1，第19页为结构统计Table S2。[取证收据](final_supplement_receipt.json)。本次没有声称穷尽所有附件。

- [最终构建体表（mmc7）](https://ars.els-cdn.com/content/image/1-s2.0-S2211124726009861-mmc7.xlsx)：6个结构构建体、12条原文变量区序列逐格核验，与预印本相同；CDR3、V/J标签及常数域标签也相同，新增2H11_mAb工作表。
- [最终克隆表（mmc3）](https://ars.els-cdn.com/content/image/1-s2.0-S2211124726009861-mmc3.xlsx)：仅用配对CDR3及双方V/J标签复核来源，不用功能值或记忆表型选模。
- [公开事实单元格](final_primary_cells.json)、[最终克隆来源收据](final_repertoire_receipt.json)。

| 结构构建体 | 最终克隆表配对CDR3 + V/J完全一致的行 | 解释 |
|---|---|---|
| NP202 | 1 | 支持天然克隆来源；同氨基酸配对的多个核酸克隆不增加独立受体数 |
| NP101 | 2 | 支持天然克隆来源；同氨基酸配对的多个核酸克隆不增加独立受体数 |
| NP207 | 1 | 支持天然克隆来源；同氨基酸配对的多个核酸克隆不增加独立受体数 |
| NP203 | 3 | 支持天然克隆来源；同氨基酸配对的多个核酸克隆不增加独立受体数 |
| NP204 | 1 | 支持天然克隆来源；同氨基酸配对的多个核酸克隆不增加独立受体数 |
| NP105 | 0 | 未找到完整双链配对；β命中的α不同，不能补天然配对来源 |

核心NP202/NP101/NP207/NP203仍有相同原文TRAV16/TRAJ53、TRBV5/TRBJ2-7标签，配对CDR3不同。混池样本不能推出独立个体，原文标签也不能自动升级成唯一IMGT等位基因调用。NP105的最终配对来源仍未闭环。

## 完整域编号与精确区间：候选已补，分歧未掩盖

对原始沉积完整构建体运行固定版本ANARCII 2.0.8与经典ANARCI 2026.2.13.2/HMMER 3.4。均为独立注释，未训练本项目模型；没有修剪输入、补残基、调低bit score阈值或修改准入代码。安装在临时独立目录，版本、包/模型文件SHA256与安装收据均保留。

ANARCII给出11/12链的编号；经典ANARCI给出8/12链。两者7条链的完整域序列和端部区间一致，但只有3条的每个IMGT位置/插入编号完全一致。四条核心β链的域端部一致，内部插入位置编号有分歧。不能把端部一致写成完整编号完全一致。

| 链 | 沉积域区间（0起点、右端不含） | 原文V-Region额外尾部 | 两工具端部 / 逐位一致 |
|---|---|---|
| NP202_alpha | [0, 118) | 无 | 否 / 否 |
| NP202_beta | [0, 113) | 无 | 是 / 否 |
| NP101_alpha | [0, 118) | N | 否 / 否 |
| NP101_beta | [0, 113) | 无 | 是 / 否 |
| NP207_alpha | [0, 118) | N | 否 / 否 |
| NP207_beta | [0, 113) | 无 | 是 / 否 |
| NP203_alpha | [0, 118) | N | 否 / 否 |
| NP203_beta | [0, 113) | 无 | 是 / 否 |
| NP204_alpha | [0, 113) | D | 是 / 是 |
| NP204_beta | [0, 116) | 无 | 否 / 否 |
| NP105_alpha | [0, 111) | N | 是 / 是 |
| NP105_beta | [0, 116) | 无 | 是 / 是 |

[完整逐残基映射及参考匹配](boundary_audit_both_methods.json)保留两套结果，含沉积实体/作者链号、每个SEQRES索引、IMGT位置和氨基酸。[ANARCII原始输出](numbering/numbering_raw.json)、[经典ANARCI原始输出与全部HMM命中](numbering/classic_raw.json)。表中的区间是至少一个工具给出的精确连续子串，不能当作两方法共识或准入身份。

具体未闭环项：

- NP202/NP101/NP207/NP203的α链：经典ANARCI的mouse_D得分略高于mouse_A，未输出允许的α编号；ANARCII给α编号。保留分歧，未强制重标或改HMM选择规则。
- NP204 β：ANARCII报“Less than 50 non insertion residues numbered.”；经典工具虽给出编号，最近J调用为TRBJ2-4、相似度约0.692，不能覆盖原文TRBJ1-6标签。
- NP204/NP105 β：现有IMGT原始参考中的TRBJ1-6*01为ORF；序列后缀能匹配ORF，不满足本项目既定功能性J参考门，未改成F。
- 部分α链原文V-Region末尾多1个N或D；工具把它置于编号域之后。这里只报告边界差异，原始构建体和输入始终保留。
- 刷新TRAV/TRBV/TRAJ/TRBJ四份IMGT参考，逐字节与旧冻结参考一致。V前10残基匹配有多个等位基因/相关基因解，β的J末10残基也有多解，不能据此宣称精确等位基因或历史V(D)J重排切点。
- 六个NP的STCRDab页面仍均无有效显式链对和完整编号。10/12链能通过独立编号后的功能性端部匹配，也不等于通过既定完整身份门；独立工具结果没有冒充STCRDab来源。

因此：可以核对候选域的精确沉积区间和逐残基编号，但完整配对编号共识、唯一V/J归属以及部分天然来源仍有缺口。历史重排的V/J切点不能只从同源氨基酸片段唯一反推。

## 跨研究正对照：搜索由44扩至282个目标

复用旧799项公共档案中281个功能性端部发现候选，并联合44个已核实身份和NP编号候选，去重后为282个目标：44个已有完整身份，238个仍仅为发现假设。385条去重链查询全部成功，双链交集377个PDB，逐个刷新论文、完整构建体和作者链号；这些PDB都已在旧799项档案内，不能称新增377个结构。

原44身份仍只有11条跨论文完整子串线索，双端固定QC后新增合格对为0。发现层另外得到197条目标—伙伴线索，涉及27个变量区假设、11个锚PDB簇；其中存在同一锚的不同候选端部，它们不是27个独立受体或197个正对照。没有把发现层标记改为完整身份。

| 发现层锚PDB | 跨论文伙伴PDB | 阻断或范围原因 |
|---|---|---|
| 1KB5 | 1KJ2 | existing alpha N-terminus quarantine and single-chain fusion context; no paired-chain admission |
| 2CKB | 1G6R, 1MWA, 1TCR | existing full numbered alpha N-terminus quarantine; endpoint-only substring variants cannot repair it |
| 2UWE | 1LP9 | existing full numbered alpha N-terminus quarantine; endpoint-only substring variants cannot repair it |
| 3HUJ | 3VWJ, 3VWK | CD1d/NKT complex; outside complete classical class-I peptide scope |
| 3TFK | 4MS8, 4MVB, 4MXQ, 4N0C, 4N5E | 42F3 compact MHC without beta2m; outside complete class-I scope |
| 5FKA | 4UDT, 4UDU | existing alpha complete J/FR4 terminal quarantine; endpoint-only alternative boundary is not verified numbering |
| 5L2K | 4G8F | GEM42 CD1b/apo context; outside complete classical class-I peptide scope |
| 5VCJ | 3O8X, 3O9W, 3QUX, 3QUY, 3QUZ, 3RTQ, 3RZC, 3TA3, 3TVM, 4IRJ, 4IRS, 4Y16, 4ZAK, 6CW6, 6CWB, 6CX5, 6CX7, 6CX9, 6CXA, 6CXE, 6CXF, 6MIV, 6MIY, 6MJ4, 6MJ6, 6MJA, 6MJI, 6MJJ, 6MJQ | CD1d/iNKT context; outside complete classical class-I peptide scope |
| 6PVC | 26SJ, 26SK, 26SL, 4IIQ, 4L4T, 4L4V, 4LCC, 4LCW, 4NQC, 4NQD, 4NQE, 4PJ5, 5U16, 5U17, 5U1R, 5U2V, 5U6Q, 5U72, 6PUC, 6PUD, 6PUE, 6PUF, 6PUG, 6PUH, 6PUI, 6PUJ, 6PUK, 6PUL, 6PUM, 6W9U, 6W9V, 7UFJ, 8Y6X, 9BTX, 9BTY, 9BTZ, 9BU0, 9C42, 9CGR, 9CGS, 9EK6, 9EK7, 9O05, 9O06, 9O07, 9O08, 9Y0T, 9Y0U, 9Y0V, 9Y0W, 9Z4R | MR1/MAIT ligand context; outside classical class-I peptide scope |
| 7N1D | 7RTR | existing full numbered beta N-terminus quarantine; 7RTR also lacks interface heavy atoms |
| 9BBC | 8TW4, 8TW6 | TCR-CD3 detergent/nanodisc complexes without required pMHC; no qualifying paired class-I unit |

[完整分组原因](broad_cross_study_summary.json)、[385条请求与查询计划](broad_search/search_plan.json)、[原始搜索结果](broad_search/search_results.json)。含所有链对和原始核验原因的2.44MB JSON按[gzip](broad_search/primary_validation_final.json.gz)保存；[包装收据](broad_search/report_packaging.json)记录原文件SHA256，可逐字恢复：

```bash
python scripts/reassemble_np_search_report.py --out /tmp/np-broad-primary-validation.json
```

已知class-I线索继续保留原N端/J端、单元范围或界面原子缺失隔离；MR1、CD1和未结合TCR-CD3背景没有改成class-I任务。完整原文/编号/双端原子证据缺失时，不能因外部搜索命中就冻结协议。当前仍为身份层64 PDB/44配对ID，固定QC 43图/31 PDB/29 ID，合格跨研究同完整αβ正对照、合格独立天然对照组、来源平衡块均为0。

## 复验、原件与限制

59项测试；十份既有冻结交付3992文件哈希保持。原始输入、身份/角色/原子门代码不变，未运行新原子QC来跳过身份门，未训练，独立验收集未访问。[复验收据](verification_final.json)、[测试日志](tests_final.log)、[交付清单](DELIVERY_MANIFEST.json)。

```bash
PYTHONPATH=src:. python scripts/verify_np_boundaries.py --replay-dir /tmp/np-boundary-replay-unique
```

最终版PDF/XLSX按出版社元数据CC BY-NC-ND 4.0本地保留；仓库提供事实收据、原始URL和SHA256，不重新发布完整原件。在线入口失败与重复安装包读取超时全部记录，均不视作阴性结果。安装包和编号权重也以版本/来源哈希追踪，不混入训练release。公开结果可离线重放；重新提取原件单元格及重新运行编号另需所列原件/固定软件。 [归档布局](archive_layout.json)逐项说明哪些原件未打包。

本轮有限范围检索不证明公共库不存在合格对照；仍缺的证据在上表逐项保留，不通过改阈值、填原子或猜测边界制造准入。
