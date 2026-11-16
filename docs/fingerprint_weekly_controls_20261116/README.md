# 新发布与修订结构复核（2026-11-16）

找到一组两端均通过固定原子 QC 的同研究、同 V 天然克隆候选。完整变量区身份规则仍阻断这组结构。合格的跨研究正对照、同 V 天然困难对照组和来源平衡块仍为 0。

## 同研究同 V 替代候选

| 候选 | 原文和坐标核对 | 裁定 |
| --- | --- | --- |
| **6ULN / 6ULR，KRAS TCR9d / 9a** | 原文说明两者来自两位受试者的肿瘤浸润淋巴细胞，Vα4/Vβ5 相同且 β CDR3 不同；沉积注释为 TRAV4*01/TRBV5-6*01。两端在官方 assembly 1 坐标通过固定图重原子和原子筛，生成 63 和 64 节点图。 | 两端 α、β 各自均在固定 V 首 10 残基核验失败。α 观察起始为 LAKTTQPISV，参考为 LAKTTQPISM；β 为 MAGVTQSPTH，参考为 DAGVTQSPTH。差异已在原始沉积构建体中出现。保留原始字母，未假定它们的天然或工程成因。暂不准入。[Sim 等原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7293613/) · [6ULN](https://www.rcsb.org/structure/6ULN) · [6ULR](https://www.rcsb.org/structure/6ULR) |
| **5SWS / 5SWZ，流感 NP2-B17 / NP1-B17** | 原文补表 1 记录同 V、不同 α/β CDR3。5SWS 通过身份端部校验和坐标 QC，属于既有入库结构。5SWZ 官方坐标内的四个配对单元均跑完筛查。 | 5SWZ 四条肽 C/H/M/R 的首位 Asn 均缺主链 N。α 的观察 FR4 末 10 也与现用功能性 J 参考不符。PDB 的 5SWS 实体说明把两条 TCR 链写作 NP1-B17，原文 NP2 的双链 CDR3 和 PDB 标题支持 NP2，实体名不能代替构建体序列。[Gras 等原文及补表](https://orca.cardiff.ac.uk/id/eprint/100009/1/Reversed%20Tcell%20receptor%20docking%20Rossjohn-Gras.pdf) · [5SWS](https://www.rcsb.org/structure/5SWS) · [5SWZ](https://www.rcsb.org/structure/5SWZ) |
| **8GVB / 8GVG / 8GVI，HIV TD08 / H25-11** | 原文记录受体源自两名受试者，均为 TRAV8-3/TRBV7-9。8GVG 是与突变肽结合的另一个 TD08 结构。 | 三份坐标分别有 3/3/1 对超阈值短接触，每份还有 α J 端和 β N 端身份阻断。[Ma 等原文](https://doi.org/10.4049/jimmunol.2200191) · [8GVI](https://www.rcsb.org/structure/8GVI) |
| **9PBG / 9PBH，AS19.2** | 两份坐标都通过原子 QC，但发表文的 Data availability 将两份都明确标为 AS19.2。两条完整沉积链逐字母相同。 | 同一受体、同一论文、两个肽；目录里 19.2 和 AS19.2 的不同显示名会误导候选分组。两个条目均缺当前有效外部身份编号。[原文](https://pmc.ncbi.nlm.nih.gov/articles/pmid/42129507/) |

对 6ULN/6ULR 的原始非对称单元，固定接触分配找不到 TCR–肽接触。RCSB 官方 assembly 1 对 D/E 使用沉积的晶体平移。本轮以官方组装提供的坐标核验，同原始文件逐项比对原子名称、残基、链及 alternate location，集合完全相同且无重复。没有补任何原子，固定图参数和原子阈值保持。非对称单元与组装结果均保留。

## 跨研究同受体核对

- **5SWZ / 7JWI** 的完整沉积 α、β 工程链分别逐字母相同，两篇 DOI 共享 4 位作者。7JWI 原始坐标通过图和原子筛，α J 端身份未通过；5SWZ 四份配对的肽主链均缺 N。7JWI 同论文的 7JWJ 通过单体身份和坐标检查，但其 α V 是 TRAV4-4/DV10；不能和该 TRAV14 谱系配成同 V 对照。[2016 原文](https://doi.org/10.1016/j.immuni.2016.09.007) · [2021 原文](https://doi.org/10.1126/science.abe9124) · [7JWI 精修来源](https://www.rcsb.org/experimental/7JWI)
- **6ULR / 7SU9** 的完整沉积 α、β 工程链分别逐字母相同，两篇 DOI 共享 4 位作者。两端有相同的 V 首段差异。7SU9 用 RCSB 官方 assembly 1 仍不能在固定规则下唯一分配完整接触单元，未标为通过原子 QC。[2020 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7293613/) · [2022 原文](https://doi.org/10.7554/eLife.75670) · [7SU9](https://www.rcsb.org/structure/7SU9)

15 条 RCSB 序列查询全部完成，覆盖 10 个双链锚。查询使用沉积完整构建体或观察编号段，结果按 PDB 双链交集记录，从沉积实体再核对精确序列。没有出现表外的新双链伙伴。3 个新修订的目录内 class-I 条目中，9YIR/9Z51 仍缺界面重原子，9Z50 缺当前外部配对编号且未通过固定链范围识别。9TMM/9TML、7N1F/7RTR、4JFE/6EQB 和四个 NP 同 V 候选的官方版本均未在本轮窗口修订。

## 覆盖范围与复验

RCSB 的 11 月 5–16 日窗口返回 274 个新发布条目和 587 个修订条目。274 个条目的标题和沉积聚合物实体逐项抓取，未检出 αβ TCR 候选；12 个有 class-I 相关实体的条目包括 10 个 KIR 复合物与 2 个无 TCR 的 pMHC。317 行 class-I 目录没有新条目；最近发布目录补出的 26JD 是 CD1c 脂质复合物，超出当前范围。合计抓取 592 条发布与目录实体记录，另取 7JWI/7JWJ 两条。15 个同研究、同 V 且不同显示名的目录组用于发现线索，包含此前已知的工程变体或待核组。

[机器可复验审计](audit.json)保留 HTTP 回执、序列、组装原子核对及失败原因。此次只处理公开结构和研究来源。尚未达到冻结协议的条件，未训练，未访问独立验收集。限定日期和序列检索无法覆盖未发布结构或元数据不足的条目。
