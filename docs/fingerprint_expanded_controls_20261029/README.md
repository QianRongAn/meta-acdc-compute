# 扩大完整双链与天然对照检索（2026-10-29）

本轮扩大了检索覆盖并补到天然克隆原文序列证据，尚无新增满足全部固定门槛的对照。身份、原子、来源要求不变；未冻结协议、未训练、未访问独立验收集。

## 实际扩大的范围

1. 完整身份锚点由上一轮29个QC身份扩大到全部44个已核实配对变量ID，包含15个没有QC通过锚的身份。87条去重链查询全部HTTP 200成功，双链PDB交集由49扩至77；逐个刷新完整构建体、作者链号和论文来源。
2. RCSB描述字段检索同时使用两个T cell receptor短语变体与TCR单词。三条查询分别返回363、363、380项，合并693个PDB；与旧799项比较，增加14个旧快照外的线索，均已核实实际分子组成。693是描述检索命中数，不能称完整TCR结构库规模。
3. 补查NP系列6个天然克隆的原文Supplementary Data 5：12条被原文标为变量区的序列均在各自沉积实体中唯一连续匹配，CDR3一致；再做11条去重链发现查询，全部成功，仅命中该研究8个已知结合/游离PDB。
4. 重新读取TCR3d/TCRdb的307项class-I表，本次与旧799项快照相比无新PDB。网页更新时间不能当作新增结构。

完整链交集77项、NP研究8项、描述新增14项合计99个不同PDB有本轮新鲜RCSB元数据。两条最初的宽泛全文查询各返回83611项，精度不足，已排除出候选选择；原始响应按gzip保存，不能当作83611个TCR。

## 完整双链跨研究线索：两端都要通过

11条完整子串跨论文线索分为以下几类。这里的子串匹配仍是发现线索，不自动证明新条目的完整域编号或接触配对。

| 组 | PDB线索 | 未准入原因 |
|---|---|
| YLQ | 7N1F→7RTR | 7RTR E28缺界面重原子；固定原子QC继续拒收。 |
| α24β17工程谱系 | 4JFE→6EQB | 6EQB通过，但4JFE未通过固定原子QC；单端通过不能组成正对照。4JFH游离项的既有阻断也保持。 |
| 42F3 | 3TFK→4MS8/4MVB/4MXQ/4N0C/4N5E | 紧凑MHC缺β2m，超出完整class-I协议范围；不补原子或改范围。 |
| BC8B | 6CUH→8DV4 | 游离/CD1b-脂质背景，不是协议要求的完整class-I肽复合物。 |
| MAIT | 4DZB→4PJ8 | 游离/MR1-5-OP-RU背景，超出协议范围。 |
| XV19/Hy19.3 | 4EI6→4ELK/4ELM | 游离/CD1d-脂质背景，超出协议范围；名称不充当身份标签。 |

在扩到无QC锚点之后，必须显式检查正对照两端。初始通用报告只标记伙伴6EQB通过，不能据此把4JFE–6EQB计为合格。权威输出为[双端核验](complete_search/primary_validation_both_ends.json)；[版本说明](complete_search/SUPERSEDED.json)保留初始报告及解释，未覆盖历史结果。补了双端检查的回归测试。

## 同来源同V天然克隆：补到了什么，仍缺什么

[Akitsu等预印本v1](https://doi.org/10.1101/2025.11.05.686789)的结构构建体表列出NP202、NP101、NP207、NP203，原文V/J标签均为TRAV16/TRAJ53与TRBV5/TRBJ2-7，四组配对CDR3不同。研究来自混池小鼠自然TCR谱系；不同克隆不等于独立个体或独立研究。表中还列NP204与NP105，它们的V组合不同，不能凑入同V组。

| 克隆 | 结合态PDB | 原文常数域构建 | 本轮结果 |
|---|---|---|
| NP202 | 9NU7 | Mouse | 原文双链序列/CDR3唯一匹配；编号与固定身份门仍未闭环 |
| NP101 | 9NV6 | Human | 原文双链序列/CDR3唯一匹配；编号与固定身份门仍未闭环 |
| NP207 | 9NV2 | Mouse and Human | 原文双链序列/CDR3唯一匹配；编号与固定身份门仍未闭环 |
| NP203 | 9NVE | Human | 原文双链序列/CDR3唯一匹配；编号与固定身份门仍未闭环 |
| NP204 | 9NW2 | Mouse | 原文双链序列/CDR3唯一匹配；编号与固定身份门仍未闭环 |
| NP105 | 9NVA | Human | 原文双链序列/CDR3唯一匹配；编号与固定身份门仍未闭环 |

这些受体不是同一亲本的CDR工程变体，但结构实验使用小鼠或人常数域嵌合构建体，应保留工程信息。不能把人常数域误写成人来源天然受体。原文TRAV16与TCR3d的TRAV16D/DV11标签没有凭名称合并为精确等位基因。

21个STCRDab请求中20个HTTP 200页面均没有有效编号/显式链对；7K0Z两次HTTP 500，作为访问失败记录。六个NP结合态均缺有效外部编号配对，继续隔离；没有跳过身份门去声称原子QC通过或失败。[逐链原文与沉积匹配、功能性端部参考命中](audit_final.json)公开保留。

RCSB引用的最终论文为[Cell Reports 2026，117908](https://doi.org/10.1016/j.celrep.2026.117908)。本次逐格序列核验针对可取到的预印本v1补充表；最终论文元数据接口两次503，不能宣称已核对最终版SI，也不能据此认定最终版没有附件。原文v1的Data Availability只预告发表时释放PDB，没有把其日期当作第二项研究。

## 旧快照外14项：逐项原因

| PDB | 实际分子 / 原因 |
|---|---|
| 1TCE | SHC plus CD3 zeta phosphopeptide; no paired alpha/beta variable domains |
| 1YGR | CD45 phosphatase plus CD3 zeta ITAM; no paired alpha/beta variable domains |
| 2E5K | Sts-1 signaling suppressor SH3 domain; not a paired TCR |
| 2H0Q | Sts-1 signaling suppressor PGM domain; not a paired TCR |
| 2IKQ | Sts-1 signaling suppressor PGM domain; not a paired TCR |
| 3RI4 | Ets1 plus TCR alpha promoter DNA; not a paired receptor protein |
| 6MF8 | TCR alpha transmembrane fragment; no complete paired variable domains |
| 6U07 | engineered alpha/beta constant domains only; no paired variable domains |
| 7K0Z | marsupial mu/gamma TCR; outside alpha/beta class-I scope |
| 9CI8 | gamma/delta-CD3 complex; outside alpha/beta class-I scope |
| 9CIA | gamma/delta-CD3 complex; outside alpha/beta class-I scope |
| 9JY2 | gamma/delta-CD3-Fab complex; outside alpha/beta class-I scope |
| 9JY3 | gamma/delta-CD3-Fab complex; outside alpha/beta class-I scope |
| 9ZZV | engineered orthogonal alpha/beta constant domains only; no variable domains |

## 当前准入和下一步证据需求

当前仍为身份层64 PDB / 44完整配对变量ID；固定QC通过43图 / 31 PDB / 29 ID。合格跨研究同完整αβ正对照=0，合格同来源同V独立天然受体对照组=0，来源平衡块=0。以上新增原文证据不自动改变准入计数。

缺口具体是：合格class-I跨研究结构对的两端原子完整性；NP系列可核验的完整域编号及精确V/J边界；天然克隆来源与当前同V定义的对应。缺编号不是认定原始结构无效，固定粗筛拒收也不等于物理无效。当前有限关键词和精确链检索不能证明其他公共资料中不存在合格对照。只有达到原门槛才可冻结协议。

## 证据、归档和复现

- [汇总裁定](audit_final.json)、[完整身份查询](complete_search/search_plan.json)、[原始命中](complete_search/search_results.json)、[双端准入核验](complete_search/primary_validation_both_ends.json)。
- [原文事实单元格](natural_primary_cells.json)、[天然克隆发现查询](natural_search/search_plan.json)、[命中](natural_search/search_results.json)。
- sources及complete_search/sources、natural_search/sources保存查询响应、HTTP状态、时间与SHA256；[归档布局](archive_layout.json)说明gzip和本地原件。
- 预印本v1正文、28.6MB附件ZIP与原始XLSX受CC BY-NC 4.0约束，仅本地保留，不随仓库重新发布；公开事实收据记录原始URL、成员名、工作表/单元格与SHA256。可按原请求获取后核对，不声称新克隆仓库包含这些原件。
- [离线复验](verification_final.json)、[测试日志](tests_final.log)、[交付清单](DELIVERY_MANIFEST.json)。

```bash
PYTHONPATH=src:. python scripts/verify_expanded_controls.py --replay-dir /tmp/expanded-controls-replay-unique
```

复验检查51项测试、原门代码与九份既有交付3814文件；不重写旧结果。原件未随仓库发布的部分仅报告本地哈希核验，公开事实收据和其余公开API响应可离线重放。
