# Meta-ACDC(计算管线)

**Deciphering the Proteome-wide TCR Cross-reactivity Landscape via Active
Meta-learning and Synthetic Immune Cells**

PhD research project — applicant Nan Wang, supervisor Prof. Sai T. Reddy,
Department of Biosystems Science and Engineering (D-BSSE), ETH Zurich.

> 本仓库为研究计划的**干实验(计算)执行仓库**:由 Claude Code 自主推进,
> 用户确认关键节点。湿实验技术路线保留于 `docs/milestones/KEY-NODES.md`,
> 不在本仓库执行范围内。所有实验数字见 `docs/benchmarks.md`。

---

## 工作内容总览

### 1. 数据管线(`src/meta_acdc/data/`, `src/meta_acdc/structure/`)

| 数据集 | 规模 | 说明 |
|---|---|---|
| VDJdb(2026-06-03 版) | 199,289 条去重记录;178,974 独特 CDR3-β;2,068 表位 | 稀疏性:17.9 万 TCR 只覆盖 2,068 表位——"泛化缺口"的实证 |
| IEDB mhc_ligand_full | 2,302,095 条递呈/结合记录(9.2GB 流式清洗) | 递呈级数据,与激活标签区分 |
| TCR-pMHC 共晶结构 | 291 个复合物(STCRDab 634 条目自动分类) | 内容型链角色分类(YFC 基序 + 长度规则),非链字母 |
| 结构 ↔ VDJdb 映射 | 2,251 个候选对(1,487 交叉反应 + 764 对照) | CDR3 窗口匹配(Cys 起点,IMGT 定义) |

关键工程决策:PDB 链字母不可靠(1AO7 的 TCR 在 D/E 链),必须按序列
内容分类;水分子/HETATM/备选构象过滤;多复合物非对称单元自动拆分
(训练集 149 → 423 阳性,×2.8)。

### 2. 模型(EGNN + 界面图)

- **界面图**:10Å 界面选择、8Å 边、RBF 距离(12 维)+ **侧链范德华接触
  直方图(4 维,C-C/C-异质/异质-异质/总数,4.5Å)**、节点特征 = 理化性质
  + pLDDT + 残基身份 + 肽链标志 + 侧链伸展度
- **EGNN**(Satorras 2021):质心居中、归一化位移、肽感知池化;E(3)
  不变性验证通过(旋转误差 <0.002)
- **诱饵构造协议**(可复现性核心):
  - 位移诱饵:肽链整体平移 6Å(几何破坏,强信号)
  - 嫁接诱饵(结构模拟):外来肽侧链经 **Kabsch 骨干对齐移植**到天然
    骨干上——骨架坐标/节点集/肽标志完全相同,只有侧链化学改变
  - 防泄漏:节点集 pin(forced_nodes)、供体序列去重(否则诱饵=天然)、
    PDB 分组切分(同复合物所有变体同折)
- **元学习**:ProtoNet + MAML++(functional_call 链式适应、逐参数可学习
  内环学习率)

### 3. 主动学习(`src/meta_acdc/active_learning/`)

- 采集函数:EIG(BALD)、预测熵、MC 方差、ε-greedy 批量采样
- 模拟闭环:未见表位池 49,696 对,4 轮 × 3,000/轮,多种子统计

### 4. AF3 交叉反应实验(`src/meta_acdc/structure/`)

- 候选清单生成(KN-5):2,251 个 TCR × 肽结构预测对,按交叉反应证据排序
- 提交批次:FASTA/JSON/纯序列三种格式(适配 AF3 网页版粘贴)
- mmCIF 解析器(AF3 输出格式,与 PDB 管线输出一致)
- 导入(import_af3.py)→ 打分(score_predictions.py)→ 每 TCR 排名
  (rank_analysis.py)

---

## 主要结果(详见 `docs/benchmarks.md`)

### ✅ 正面结果

1. **侧链接触特征解锁化学相容性学习(v9.1)**:
   嫁接诱饵任务(几何不变、化学改变)AUROC **0.500 → 0.933**;
   位移诱饵 0.946 → 0.980。消融:去掉接触特征回落 0.501(因果确认)。
2. **主动学习**:EIG 采集在 6% 采样处 2.0× 阳性召回增益(10.4% vs 随机
   5.9%);EIG > MC 方差 > 熵 ≈ 随机。
3. **AF3 预测结构上的交叉反应排名**:A6 TCR 的 8 个候选肽中,5 个
   VDJdb 验证的交叉反应肽得分 0.96-0.99;B7 TCR × Tax 同源肽 0.999。

### ⚠️ 有意思的阴性结果(全部记录在案)

1. **序列方法天花板**:ProtoNet/MAML++/逻辑回归全部聚在 AUROC
   0.62-0.64(留出表位),更大模型/更多训练回合无法突破——与 TITAN
   (0.62)/ ERGO TPP-III (0.669) 文献一致。**结构信息是唯一出路。**
2. **Cα 分辨率的信息论不可分**:嫁接诱饵在 Cα 级图上 AUROC 0.500,
   **3 倍数据(149→423)不改善**——判别信息根本不在 Cα 表征里。
3. **三连排除**:肽感知池化(0.498)、残基身份嵌入(0.499)、序列+图
   混合(0.905 < 纯 EGNN 0.933)——常见改进全部无效,**侧链接触才是
   关键变量**。
4. **单对绝对打分失败(MAG-IC3/titin 图灵测试)**:5BRZ(MAGE-A3)0.605
   vs 5BS0(titin)0.223——模型未在单对比较中检出致死结构模拟。由此
   定义正确的问题设定:**每 TCR 的候选肽相对排名**,而非绝对分。
5. **熵采样在弱模型上无效**(≈随机)——均值预测的熵不含模型不确定性
   的差异信息,支持选 EIG。
6. **构造陷阱(诚实性记录)**:v9 首训 0.988 系诱饵 atoms 字段泄漏
   (诱饵残基无原子 → 接触特征暴露身份),已回退并修复;**AF3 首次
   17 个任务因五链被拼成单链全部作废**(网页版粘贴问题,已用多链
   FASTA 解决)。这些教训写入了 `docs/benchmarks.md` 供论文
   Methods 部分引用。

---

## 目录结构

```
docs/
  proposal/      解码后的计划书全文
  literature/    文献综述×3 + BibTeX 引用库
  milestones/    关键节点与里程碑追踪(干实验执行计划 + 完整技术路线)
  benchmarks.md  全部实验数字(正面+阴性结果)
  paper1_outline.md 论文1大纲(Meta-TCR-GNN)
  paper2_outline.md 论文2大纲(主动学习安全评估)
  dataset-report.md 数据集统计报告
src/meta_acdc/
  data/          VDJdb/IEDB 下载清洗、统一 schema、防泄漏切分、临床金标准案例
  structure/     界面图构建、链分类、mmCIF 解析、复合物拆分、AF3 导入/打分/排名
  models/        EGNN、混合模型、元学习(ProtoNet/MAML++)、硬负样本、训练脚本
  active_learning/ 采集函数、模拟闭环
scripts/         环境复现(setup_env.sh)
data/processed/af3_starter|af3_json|af3_plain/  AF3 提交任务(启动包/JSON/纯序列)
```

## 复现

```bash
bash scripts/setup_env.sh          # venv + torch 2.9.1+cu126 + 科学栈
.venv/bin/pip install -e . --no-deps
.venv/bin/python src/meta_acdc/models/dataset.py      # 构建界面图数据集
.venv/bin/python src/meta_acdc/models/train_egnn.py   # 训练 EGNN(约10分钟/150epoch, 4GB GPU)
.venv/bin/python src/meta_acdc/structure/score_predictions.py  # 给 AF3 结构打分
```

## 环境

GTX 1050 Ti 4GB / 8GB RAM / Python 3.14(torch 2.9.1+cu126)。结构预测走
AF3 官方网页(用户手动提交,条款合规);TCRmodel2 为条款安全的替代方案
(TCR-pMHC 基准上 DockQ 0.566 优于 AF3 0.499)。
