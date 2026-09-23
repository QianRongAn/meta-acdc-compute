# 论文 2 大纲:AI 驱动的蛋白组尺度 TCR 交叉反应安全评估(干实验版)

> 目标期刊:Bioinformatics / Nature Methods 级别(取决于 AF3 数据规模)
> 状态:骨架 v0(2026-09-21)。范围:纯干实验(湿实验技术路线保留但不执行,
> 见 KEY-NODES.md)

## 题目(候选)

"An AI-driven dry-wet loop simulation for proteome-wide safety assessment of
therapeutic TCRs: how few experiments suffice?"

或:"Active learning on structural compatibility models: ranking cross-reactive
off-targets for therapeutic TCRs"

## 故事线

1. **Motivation**:TCR-T 疗法的致死交叉反应(MAGE-A3/titin、MAGE-A12)在
   体外常规筛选不可预见;蛋白组尺度实验筛选不可能(10^7-10^15 组合)。
2. **方法**:
   - 递呈预筛(MHCflurry 本地,NetMHCpan 在线版已废弃):蛋白组 → ~5M 候选压缩
   - 结构预测(AF3/TCRmodel2)+ EGNN 相容性打分(论文 1 的模型)
   - 主动学习:EIG/熵/方差采集函数 + ε-greedy,逐轮"验证"(以公开数据
     作 oracle 的模拟闭环)
3. **已有结果**:
   - 模拟闭环:6% 采样 2.0× 阳性召回增益(v0,序列模型)——待结构模型升级后
     重跑
   - 采集函数对比(进行中):EIG vs 熵 vs 方差 vs 随机
4. **进展(AF3 数据到达后)**:
   - [x] 每 TCR × 多候选肽的结构预测 → 稳定 DA EGNN 打分 → 排名
     (144-native 域适应后跨实例 STABLE,mean Spearman 0.729)
   - [x] 金标准检查:致死性 MAGE-A12 模拟肽 KVAKELVHFL 排到靶肽之上
     (FATAL-mimicry flagged);titin ESDPIVAQY 与靶肽分差 0.02(§2.8)
   - [ ] 模拟闭环升级:oracle = EGNN 自身 vs 金标准标签的对比
   - [ ] 阈值标定:风险分与临床安全边界的映射(计划书 Module 3 的
     "Enhanced Safety Margin":高危器官表达蛋白阈值 0.5→0.3)
5. **交付物**:TCR-Safety-Radar 看板(KN-13)、打分管线、排名报告

## 结构草案

- **Introduction**:安全评估的搜索空间问题 + 主动学习文献支撑
  (ALDE/ALSEBO/EVOLVEpro:少量实验达到数倍效率)
- **Results**
  1. 模拟闭环的效率曲线(采样比例 vs 阳性召回,多采集函数)
  2. 结构打分对排名的贡献(序列模型 vs 结构模型在闭环中的增益)
  3. 金标准案例:titin/MAGE-A12 的排名(KN-5 数据)
  4. 阈值与假阳性率的权衡(安全边界分析)
- **Discussion**:干实验模拟的适用边界 vs 真实 ACDC 闭环(技术路线保留);
  对临床前安全评估工作流的含义
- **Methods**:模拟器协议(oracle、批量、ε-greedy、多种子统计)

## 实验资产

| 资产 | 位置 |
|---|---|
| AL 模拟器(多模式+多种子) | src/meta_acdc/active_learning/simulate_al.py |
| 采集函数库 | src/meta_acdc/active_learning/acquisition.py |
| EGNN 打分管线 | src/meta_acdc/structure/score_predictions.py |
| 金标准临床案例 | data/processed/clinical_gold_standard.tsv |
| KN-5 候选清单(2,251 对) | data/processed/kn5_submission_list.tsv |
| AF3 批次(50 任务) | data/processed/af3_batch1/, af3_plain/ |
