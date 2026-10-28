# 完整双链与天然对照续查（2026-10-29）

本轮完成有限范围的公开结构开发资料续查，没有新增满足原标准的跨研究正对照或独立天然受体结构对照。协议仍未冻结，未训练，未访问独立验收集。

## 搜索范围与复核

以现有固定QC通过的29个完整配对变量ID为锚，加上AS8.4完整双链及AS3.1、AS4.2、AS4.3沉积片段，共33个目标、66条去重链查询。RCSB序列服务采用100%局部一致、返回全部命中；66条均HTTP 200成功。双链按PDB求交只用于发现线索，再对49个PDB刷新完整沉积序列、作者链号、论文DOI和作者。所有可与既有裁定比较的构建体均未改变。

局部一致不能替代完整变量域、实际链配对或原子QC。AS片段不补残基、不升级身份标签；实体/链的组合只是待审线索，不能据此宣称接触配对。

## 逐项结论

| 线索 | 本轮结果 / 未准入原因 |
|---|---|
| AS8.4 | 仅8CX4；已有固定原子筛查失败保持。 |
| AS3.1 | 仅7N2S；沉积端部缺口未补齐。 |
| AS4.2 / AS4.3 | 仅7N2N/O与7N2P/Q/R；片段命中未找到完整端部补救。 |
| 6VRM / 6VTH | 同研究；6VTH新鲜元数据仅两个TCR实体，为未结合结构，不能组成所需pMHC正对照。 |
| 7N1D、5YXN、7DZM | 局部链命中，但不含完整查询β序列；不能当作同完整αβ。 |
| 6AT6 / 6AVF | 同研究，6AT6为未结合KFJ5；Aα/Bβ裁定保持。 |
| 7N5C / 7N5P | 原文6218αC与亲本6218工程关系保持，不能计作独立天然受体；7N5P固定QC失败保持。 |

完整双链跨研究线索仍只有以下3个，均是旧库存中已知阻断项：

| 锚PDB | 命中PDB | 固定原因 |
|---|---|
| 7n1f | 7rtr | incomplete interface heavy atoms: [('E', 28, ['CD', 'CG', 'OE1', 'OE2'])] |
| 6eqb | 4jfe | existing fixed native atomic screen rejected unit |
| 6eqb | 4jfh | outside_complete_class_I_alpha_beta_scope |

7RTR来源作者不重叠，但界面重原子缺失；4JFE/4JFH还与6EQB有作者重叠。不同论文不能自动视作独立来源。原子筛查是本项目固定准入规则，其拒绝不等于晶体物理无效或MolProbity结论。

## 当前状态与边界

身份层仍为64 PDB / 44配对变量ID；固定QC通过43图 / 31 PDB / 29 ID。合格跨研究同完整αβ正对照=0，独立天然受体同来源同V合格组=0，来源平衡块=0。AS3.1/AS8.4的原文患者与同V证据保留，但尚不构成通过结构QC的对照组。

这次搜索覆盖上述33个目标，不能推断全世界不存在合格对照。需要新公开坐标或原始完整序列及来源证据解决当前缺口，达到原标准后才可冻结评估协议；不据此寻找模型正结果。

## 文件与复现

- [查询计划](search_plan.json)、[原始请求](requests.json)、[搜索结果](search_results.json)。
- [最终元数据裁定](primary_validation_final.json)；primary_validation.json为中间版本，以final为准。
- [原文构建体与天然来源证据](../fingerprint_construct_followup_20261027/README.md)。
- sources/search与sources/metadata保存原始响应及时间、状态、SHA256收据。
- [离线验证](verification_final.json)、[测试日志](tests_final.log)、[交付清单](DELIVERY_MANIFEST.json)。

```bash
PYTHONPATH=src:. python scripts/verify_native_search.py --replay-dir /tmp/native-search-replay-unique
```

47项测试通过；八份既有冻结交付3731个文件哈希未变。验证只重放公开开发资料，不运行模型或读取独立验收数据。
