# TCR 端部隔离项续查与 GPC3 天然克隆对照（2026-11-16）

本轮定位到同研究、同 V 组合、同 pMHC 的两个天然克隆。两结构均有界面重原子缺失，固定图准入拒收。跨研究完整双链正对照、合格天然困难对照组和来源平衡块仍各为 0。

## GPC3 天然克隆与原始坐标

| 官方坐标 | 克隆与原始实验 | 缺失界面重原子的残基 |
| --- | --- | --- |
| [9TMM](https://www.rcsb.org/structure/9TMM) | TCR-A，1H4，实验 1 | HLA A:ARG74、TCR β D:GLN29 |
| [9TML](https://www.rcsb.org/structure/9TML) | TCR-B，3E9，实验 2 | HLA A:GLU59、A:GLU155、TCR β D:ASP55 |

[克隆原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12266797/)及[补表 3、4](https://cdn-links.lww.com/permalink/hep/j/hep_2024_11_22_stubbs_hep-24-1558r1_sdc2.pdf)第 3 页将两条受体追溯到分别分选的 HHD-DR1 小鼠 CD8 T 细胞。两者同为 TRAV9-2 / TRBV1，α CDR3 和 J 基因不同，β CDR3 氨基酸相同，连接区核苷酸、D 基因不同。原文 CDR3 分别唯一出现在 PDB 构建体的 C（α）和 D（β）链。两 PDB 的 β 构建体完全相同，α 构建体不同，肽 FLAELAYDL、HLA 和 β2m 沉积序列逐一相同。结构来自[同一项研究](https://doi.org/10.64898/2026.09.15.750028)。

TCR3d 的 α V 注释为 TRAV9D-4，沉积 SOURCE 字段标为 Homo sapiens。原文提供受体小鼠来源。保留 V 名称差异，精确等位基因同义关系仍待核。STCRDab 两页均未给出显式链对和编号。补表可核实 V/J、CDR3 和 β 连接区。完整天然重排序列和现用外部编号证据仍有缺口。结构预印本全文请求限流，克隆原文与补表已读取。

固定 class-I 接触分配在两端都选中 C/D、肽 P、重链 A 和 β2m B。几何短接触子筛通过，图构建的界面重原子门两端失败。所有缺失位点均与原始 PDB 的 REMARK 470 相符，逐原子清单见 [audit.json](audit.json)。即使补齐身份编号，这两份坐标仍会被固定门拒收。

以 A/B 两条沉积构建体和 A/B/C 原文 CDR3 查找替代坐标，A、B 只找到上述对应 PDB，C 无双链命中。两者识别同一 pMHC，属于受体身份检索的困难对照候选，未添加实验非结合标签或把功能读数纳入 SPR。

## 跨研究正对照续查

从旧设计 267 条端部隔离记录提取 121 个唯一观察双链片段对。215 条 RCSB 查询全部成功，其中复用 46 条同参数回执，新增 169 条，70 次初始网络超时均重试成功。218 个 PDB 有双链查询命中。核验沉积实体后，216 个在不同实体上有精确双片段。片段命中只作为发现线索。

其中 187 个已在原设计，29 个在设计之外，包括曾独立核过的 A6 Fab 和 TK3。逐项理由见[候选裁定](candidate_decisions.json)，原始证据见[查询结果](search_results.json)。

- **3-2E，8I5D / 9WFS。** 双片段准确命中。两端主要引文缺 DOI/PMID，共有五位作者，9WFS 缺有效外部配对编号。固定原子筛分别发现 1 对与 2 对超阈值重叠。9WFS 同一沉积研究的 9WFT 使用另一 V 组合。[8I5D](https://www.rcsb.org/structure/8I5D) · [9WFS](https://www.rcsb.org/structure/9WFS) · [9WFT](https://www.rcsb.org/structure/9WFT)
- **SN45，7BYD / 9WNT。** 两个 DOI 共享三位作者，新结构使用 7BYD 作为初始精修模型。两端均为四残基脂肽，固定完整肽复合物范围拒收。9WNT 同文另外六项沉积无 TCR。[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12804144/) · [精修来源](https://www.rcsb.org/experimental/9WNT)
- **A6c134，4FTV / 9NDS。** 9NDS 使用替代 β2m 的小型化 MHC，缺独立 β2m。[官方沉积和原文摘要](https://www.rcsb.org/structure/9NDS)
- **11OE、9ZEP/9ZEQ/9ZER/9ZES、9ZET。** 固定范围拒收单链融合 pMHC。germline-like DMF5 与 TCR156 的 CDR3 也经过改造。[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13648438/)

新取 TCR3d 目录有 317 行，10 项标记于上一轮后发布。本轮未下载其中两条新工程 A3A 的结构，也未用其做临床验收。

## 交付与复验

固定库存仍为 64 PDB / 44 个完整配对身份，43 图 / 31 PDB / 29 个身份通过固定 QC。本轮新增准入 0，未冻结协议或重训，未访问独立验收集。公开结构/API 快照随交付归档。完整论文、专利页面和补件保存在本地，仓库提供回执、摘要事实、源定位与 SHA256。样本缺口是具备跨研究独立来源、完整编号和原子完整性的结构，当前限定查询不能排除其他公共资源中有合格对照。

    PYTHONPATH=src:. python scripts/verify_terminal_rescue.py --out /tmp/terminal-rescue-replay-unique

本轮重放的检索和坐标审计结果逐字节一致，38 项既有相关测试通过，并复核了 18 批历史交付的 4,795 个文件。历史 NP105 验证链依赖的 Biopython 和 ANARCI 当前环境未提供，本轮未重跑完整历史测试链。详见[复验记录](verification.json)。
