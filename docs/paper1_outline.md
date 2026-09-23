# 论文 1 大纲:Meta-TCR-GNN(方法学论文)

> 目标期刊:Nature Machine Intelligence / Bioinformatics
> 状态:骨架 v0(2026-09-21),实验数据实时填充中
> 本文件随实验推进持续更新,是论文写作的工作底稿

## 题目(候选)

"Geometric deep learning of TCR-pMHC interfaces reveals that side-chain
contact information — invisible at Cα resolution — is required for
chemical-compatibility learning"

或简洁版:"Atom-contact graph networks for TCR-pMHC interface compatibility
prediction"

## 故事线(已有实验支撑)

1. **Motivation 1 — 序列天花板**:序列级方法(逻辑回归 0.609;ProtoNet
   0.635;MAML++ 0.625;TITAN 0.62 / ERGO TPP-III 0.669 文献一致)聚在
   0.62-0.64,信息瓶颈明确。
2. **Motivation 2 — 结构预测现状**:AF3 对 TCR-pMHC DockQ 0.499、TCRmodel2
   0.566(STCRDab-22)——现成结构预测不足以直接作交叉反应判据,需要
   学习型几何评分。
3. **核心发现(本论文的增量贡献)**:
   a. **Cα 分辨率的信息论诊断**:嫁接诱饵(几何不变、化学改变)在 Cα 级
      图上不可分(AUROC 0.50),且 3 倍数据(149→423)不改善——判别信息
      不在 Cα 表征中。
   b. **侧链接触特征是钥匙**:加入侧链 vdW 接触直方图(4.5Å,4 维)后,
      嫁接任务 0.500 → 0.933;位移任务 0.946 → 0.980。
   c. **对照谱系**:纯序列 0.495 → 序列+图统计 0.649 → Cα-EGNN 0.500(嫁接)
      → 侧链-EGNN 0.933。
4. **方法**:EGNN(Satorras 2021)+ 界面图(10Å 界面/8Å 边)+ 侧链接触
   特征 + 严格的诱饵构造(Kabsch 骨干对齐侧链移植、节点集 pin、
   供体序列去重、PDB 分组切分)。
5. **与临床交叉反应的联系**:MAGE-A3/titin(Raman 2016,骨架 RMSD 0.285Å
   的"直接分子模拟")——侧链互补性正是我们特征建模的对象;金标准
   验证集已建(clinical_gold_standard.tsv)。

## 结构草案

- **Title / Abstract**
- **Introduction**:TCR 交叉反应安全危机(MAGE-A3 死亡病例)→ 为什么
  结构方法必要(序列天花板 + 结构预测现状)→ 我们的问题设定
- **Results**
  1. 数据与表征:界面图、侧链接触特征的定义与验证
  2. Cα 分辨率的不可分性诊断(三连实验 + 数据量消融)
  3. 侧链接触特征解锁化学相容性学习(0.50→0.93)
  4. 对照谱系与消融(特征消融:去接触直方图、去残基身份、去侧链伸展度)
  5. AF3 预测结构上的真实交叉反应排序:域适应扩样(63→89→144 天然)
     后跨实例判决由 UNSTABLE 转为 STABLE(mean Spearman 0.729);
     临床金标准致死 off-target 被正确标记(§2.7)
  6. (待补)与序列 SOTA(pMTnet-omni/PanPep)在真实标签上的对比
- **Discussion**:表征层级 vs 数据规模的争论;对"结构预测+学习评分"
  范式的含义;局限(Cα 之外的缺失:静电表面、主链动力学)
- **Methods**:数据管线、EGNN 细节、诱饵构造协议(可复现性核心)、
  切分协议、指标(AUPRC 为主)

## 已就绪的实验资产

| 资产 | 位置 |
|---|---|
| 数据管线(VDJdb/IEDB/结构集) | src/meta_acdc/data/, structure/ |
| EGNN + 侧链接触特征 | src/meta_acdc/models/egnn.py, structure/graph.py |
| 诱饵构造协议(Kabsch 移植等) | src/meta_acdc/models/dataset.py |
| 基准记录(全部数字) | docs/benchmarks.md |
| 金标准临床案例 | data/processed/clinical_gold_standard.tsv |
| AF3 预测评分管线 | src/meta_acdc/structure/score_predictions.py |

## 待补实验(优先级排序)

- [x] 特征消融(去接触直方图 → 回落到 0.501;去侧链伸展度 → 小降;
      见 benchmarks 特征消融谱系)
- [x] AF3 结构上的真实交叉反应排序(2026-09-23:144-native 域适应后
      跨实例判决 STABLE,mean Spearman 0.729;论文 §2.4)
- [ ] 与 pMTnet-omni/PanPep 在同数据上的公平对比(需其代码复现;
      当前以内部基线逻辑回归 0.609 / ProtoNet 0.635 / MAML++ 0.625 对照)
- [x] 结构预测器对比(AF3 vs TCRmodel2:DockQ 0.499 vs 0.566;
      TCRmodel2 路线因服务器队列废弃,见 README 第 15 步)
- [x] 图灵测试式验证(2026-09-23:临床扫描把致死性 MAGE-A12 模拟肽
      KVAKELVHFL 排到靶肽之上、titin ESDPIVAQY 与靶肽分差 0.02;
      论文 §2.7 / KN-11)
