# N2／N8 原子数据修复与重训阻断处置（2026-10-17）

**旧污染控制已停用，天然结构资料已重新构建并独立核验；结构分类重训仍不放行。** 当前交付物是无结合标签的结构资料库和很小的构造控制试运行，不是合格的结合／非结合训练集。

## 实际处置

| 原对象 | 本轮处理 | 结果 |
|---|---|---|
| 9 个存在界面重原子缺失的天然结构 | 完整隔离相应结构，保留原始文件、18 个缺失残基明细及来源；不删除残基伪装完整，也不猜补原子冒充实测 | 全部未进入新结构资料库 |
| 300 个固定方向位移控制 | 从当前构建与训练入口退役，旧函数只保留作事故重建 | 不再作为当前训练样本 |
| 288 个逐残基嫁接控制 | 整个拼接协议退役，不只筛去其中 146 个极短距离案例 | 不再作为当前训练样本 |
| 300 份历史天然结构记录 | 重新核对来源、链身份、保留界面原子、肽连续性及预设接触筛查 | 168 份进入无标签资料库：145 晶体、23 预测；其余 132 份保留隔离原因 |
| 新的整体天然肽构象转移 | 同分区、同域、相同观测 MHC/B2M 序列、同肽长度且肽序列不同；一次正规旋转／平移，拒绝反射、断链和新碰撞 | 仅 5 对通过工程检查，涉及 3 个受体家族 |

132 份未纳入记录不等于 132 份都被证明“结构错误”：97 份触发当前保守接触筛查，另有 9 份界面原子缺失、25 份预测身份／参照 QC 未满足、1 份会改变历史参照归属。筛查不能代替实验密度、力场或完整结构验证。

历史 888 张图的逐项去向：588 张旧合成控制退役、168 份天然资料保留、132 份天然记录隔离。见 [完整处置账本](structure_repair_20261017/historical_graph_disposition.json)与[原始 N8 缺失明细](structure_repair_20261017/missing_interface_quarantine.json)。原 v4 与之前版本、源结构、历史结果、日志和 checkpoint 均未覆盖。

## 新控制的范围与限制

新方法整体转移供体天然肽的所有重原子。单个正规刚体变换保持供体内部原子距离、键角、二面角和手性，避免逐残基侧链拼接引入新的内部几何。对接后的肽与受体环境仍可能不相容，须通过相同的原子接触检查。保留天然构象不保证新环境中的构象正确，更不能证明它“不结合”。

原生与控制使用同一套完整性和接触规则。两者使用各自界面节点的并集，并要求并集内两边都没有缺失原子；图节点身份、肽 mask 完全相同。两类所有 pLDDT／缺省置信度特征统一屏蔽，避免供体或数据域的置信度模式泄漏类别。

接触检查保留原有 1.2 Å 极短距离拒收，并增加按元素半径计算的非供体／受体原子范德华重叠检查（大于 0.4 Å 拒收）。潜在氢键仅按标准残基的化学供体／受体识别，且重原子间距至少 2.4 Å；没有加氢和角度验证，所以不称为真实氢键或 MolProbity clashscore。[当前方法声明和固定参数](structure_repair_20261017/policy_v2.json)。

方法依据：[MolProbity 原始方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC2803126/)区分非供体／受体碰撞与氢键接触，并依赖加氢等检查；半径来自 [Bondi 原始工作](https://doi.org/10.1021/j100785a001)。本实现只是重原子工程筛查，未运行 MolProbity、旋转异构体库验证或能量优化，不能宣称通过完整物理验收。

## 失败过程也保留

第一轮 v5 把所有重原子对都零容忍套用同一重叠条件，错误地把可能的紧密氢键一起判为碰撞，最终 0 对控制；[原参数](structure_repair_20261017/policy.json)和[失败统计](structure_repair_20261017/failed_v1_build_stats.json)原样保留。

v2 修正了接触化学分类，1.2 Å 和非供体／受体 0.4 Å 数值阈值未放宽。随后 v6 的可选 QC 错误列导致导出失败；修复为显式列并集，再在新目录完成 v7。所有失败目录和日志保留，未覆盖旧输出；没有根据 AUROC、种子相关性或模型表现调规则，也没有偷偷重分测试集。

## 独立复核与放行判决

本地资料包：`data/releases/20261017-structural-v7-atomic-repair/`。

- `native_inventory_graphs.pt`：168 份源原子可复现的无标签图，不含可用于结合训练的 response labels。[独立核验结果](structure_repair_20261017/native_inventory_validation.json)。
- `graphs.pt`：5 个天然／构造控制对，仅作构造协议试运行。独立验收从源原子重新生成图，并检查完整性、接触、刚体变换、源链、哈希、节点并集、置信度屏蔽、供体条件及家族归属；工程检查通过。[成对复核](structure_repair_20261017/paired_validation.json)。
- 冻结分区：训练 1 对、验证 4 对、校准 0 对、测试 0 对。没有校准／测试样本且只有 3 个受体家族，**不能恢复结构分类重训**。未挪动已有样本去填空分区。
- `asset_manifest.json` 状态固定为 `inspection_only`，所有训练许可均为 false；普通 release loader 拒绝它。没有发布合格训练 manifest，也没有执行 optimizer 或启动调度。
- 122 项测试：119 通过，3 项旧本地夹具缺失跳过。[测试日志](structure_repair_20261017/unit-tests-final.log)。6 个原始关键产物、旧 v4 全部 manifest 文件、全部历史源结构哈希均核对不变。[保全记录](structure_repair_20261017/preservation.json)。

完整性核验命令（复用现有环境）：

```bash
PYTHONPATH=src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  /home/administrator/docs/meta-acdc/.venv/bin/python scripts/verify_atomic_repair.py \
  --bundle data/releases/20261017-structural-v7-atomic-repair \
  --out NEW_ATOMIC_VERIFICATION.json
```

源结构可用时可重建到新目录：

```bash
PYTHONPATH=src /home/administrator/docs/meta-acdc/.venv/bin/python \
  -m meta_acdc.data.build_release \
  --previous-release-dir data/releases/20261017-structural-v4 \
  --out NEW_ATOMIC_REPAIR_DIRECTORY
```

## 仍须解决的研究条件

旧“困难阴性”假设没有通过这次处置变真。整体转移控制也可能真实结合，只能检查构造协议，不能充当实验不结合标签。要训练真实结合模型，应优先为已有 SPR／功能观测建立受体身份准确、生成方式一致、几何 QC 可比的结构输入，并另留未参与开发的独立受体／家族验收来源。当前清理完成的是污染处置和可复核资产，未完成这项新的生物学数据集。
