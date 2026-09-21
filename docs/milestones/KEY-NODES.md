# 关键节点与里程碑追踪

> 项目启动:2026-09-20(计划书 Year 1, Month 1)。计算侧由 Claude Code 自动推进;
> 标注 🔴 的节点 = **需要用户确认/操作**,到点停下等确认。
>
> **范围(2026-09-21 用户定):只做干实验。** 湿实验环节(寡核苷酸池订购、
> ACDC 实验、原代 T 细胞、NSG 小鼠)全部移除;干-湿闭环改为**主动学习
> 模拟器**(公开数据验证)作为独立交付物。

## 当前状态

- [x] KN-0:任务范围确认(执行研究计划)+ 起点确认(从零开始)— 2026-09-20 ✅
- [x] KN-1:计算环境搭建 ✅ 2026-09-21 — torch 2.9.1+cu126 GPU 冒烟通过(GTX 1050 Ti)
- [x] KN-2:文献综述 v0 ✅(3 篇综述 + 3 个 bib;🔴 待用户审阅)
- [x] KN-3:VDJdb ✅(199,289 条)+ IEDB ✅(2,302,095 条,1,394,080 肽,544 等位基因)— 2026-09-21
- [x] KN-4:EGNN 原型首训 ✅ 2026-09-21 — 诱饵判别 AUROC 0.771 / AUPRC 0.556(PDB 分组切分,无泄漏);同任务对照谱系:纯序列 0.556 < 序列+图统计 0.649 < EGNN 0.771
- [ ] KN-5:AlphaFold Server API(🔴 需用户提供)+ 首批结构批量预测
- [x] KN-6:v0 ✅ 2026-09-21 — 硬负样本三重筛选模块 + ProtoNet 少样本 AUROC 0.635(留出表位,超序列基线 0.609)
- [x] KN-7:v0 ✅ 2026-09-21 — 主动学习模拟器:EIG+ε-greedy 在 6% 采样达 2.0× 阳性召回增益
- [ ] KN-7+:AL 模拟器升级(EGNN 打分 + ACDC 预筛肽库分布)
- [ ] KN-4+:VDJdb 标签 ↔ 结构映射(289 复合物 TCR 序列与 VDJdb 交叉匹配)
- 项目仓库:`/home/administrator/docs/meta-acdc/`(18 commits)

### 文献综述关键发现(2026-09-20)
1. ACDC 原文确认:bioRxiv 2025.07.23.666264(Hong K.-L. 等 14 人,Reddy 通讯)。**缩写展开与计划书不一致**:原文是 "antigen-presenting cell detecting cytokine",计划书写 "Artificial Cell Display of pMHC" → 需与实验室核对
2. NetMHCpan 最新是 **4.2**(2025,Front Immunol),计划书引用的 4.1 已过时 → 预筛管线改用 4.2
3. 直接方法学基线:**MCGLPPI++**(JCIM 2025,几何 GNN 少样本 PPI,含 TCR-pMHC 簇验证)、**ZeroBind**(MAML++)
4. vanilla MAML 有任务冲突/过度记忆陷阱 → 实现 MAML++ / ProtoNet
5. 主动学习文献规则"批数多、每批少、收敛快" → 模拟器需做批大小消融

## 关键节点总表

| # | 节点 | 计划时间(计划书) | 实际日期 | 谁操作 | 通过标准 |
|---|------|-----------------|---------|--------|---------|
| KN-1 | 🔴 计算环境:装 Python 3.12 venv(uv)+ PyTorch | 2026-09(Month 1) | 待确认 | 用户批准安装 | torch 在 1050 Ti 上跑通 CUDA 冒烟测试 |
| KN-2 | 文献综述 v0 + 参考文献库 | Month 1-3 | 进行中 | AI | 3 篇综述文档 + bibtex 入库 |
| KN-3 | 数据管线 v0:VDJdb/IEDB 下载、清洗、统一格式 | Month 1-3 | — | AI | 标准格式数据集 + 统计报告 |
| KN-4 | 🔴 EGNN 原型跑通首个基准(VDJdb vs NetMHCpan/ERGO) | Month 4-6 | — | AI 训练/用户验收 | AUPRC 达到文献可比水平 |
| KN-5 | 🔴 AlphaFold Server API 接入(需用户注册/提供 key)+ 首批结构批量预测 | Month 4-6 | — | 用户提供 API | 首批 100 条 TCR-pMHC 结构入库 |
| KN-6 | 硬负样本挖掘模块 + 元学习基线(MAML/原型网络) | Month 4-6 | — | AI | 三重筛选管线 + few-shot 评测 |
| KN-7 | 主动学习模拟器(公开数据干跑,3-5 轮 EIG 采样) | Month 6-9 | AI | 模拟显示 ≤10% 实验量达 >90% 覆盖 |
| KN-8 | ~~5 万寡核苷酸池订购~~ **已移除(纯干实验)** | — | — | — |
| KN-9 | ~~ACDC 湿实验~~ **已移除** → 替代:模拟器升级(结构打分 + 多种采集函数对比) | Month 9-15 | AI | 模拟闭环 3-5 轮完整报告 |
| KN-10 | ~~NGS 管线~~ **已移除** → 替代:公开多数据集统一评测(VDJdb/IEDB/结构集) | Month 13-18 | AI | 跨数据集评测报告 |
| KN-11 | 🔴 首批 TCRmodel2/AF3 结构批量预测提交(KN-5 候选清单 2,251 对,优先前 ~100) | Month 4-6 | 用户协助提交 | 首批结构入库 |
| KN-12 | ~~原代 T 细胞验证~~ **已移除** → 替代:临床 TCR 深扫描的 in silico 部分(Kimmtrak 等,公开数据) | Month 25-27 | AI | off-target 风险排序清单 |
| KN-13 | 风险看板 v1(Web 交互式"安全雷达") | Month 28-30 | AI | 部署可用 |
| KN-14 | 第一篇论文(Meta-TCR-GNN 方法学,含 Cα 分辨率诊断实验) | Month 22-24 | AI 起草/用户投稿 | 投 Nature MI / Bioinformatics |
| KN-15 | 第二篇论文(主动学习模拟器 + in silico 安全评估) | Month 34-36 | AI 起草/用户投稿 | 投 Bioinformatics / Nat Methods 级别期刊 |
| KN-16 | 交付物打包:软件套件 + Docker + 数据集 | Month 37-39 | AI | 开源仓库发布 |

## 各阶段执行清单

### Year 1 Q1(2026-09 ~ 2026-11,Months 1-3):整合与文献综述
- [ ] 文献综述:平台(进行中)/ ML 模型(进行中)/ 数据库与主动学习(进行中)
- [ ] 环境搭建(KN-1)
- [ ] 数据管线 v0(KN-3)
- [ ] ACDC 原始数据盘点:向实验室要历史数据清单(🔴 需用户对接实验室)

### Year 1 Q2(Months 4-6):Meta-TCR-GNN 原型
- [ ] EGNN 骨架实现 + 接口图构建(10Å 截断)
- [ ] AlphaFold 3 结构预测管线(KN-5)
- [ ] 硬负样本三重筛选(Levenshtein ≤3 / 静电势 Pearson >0.8 / 不确定性 0.4-0.6)
- [ ] 元学习框架(MAML / Prototypical)
- [ ] 基准评测(VDJdb、IEDB;对照 NetMHCpan、ERGO、NetTCR)

### Year 1 Q3(Months 7-9):基准测试 + 首批文库下单
- [ ] 公开数据基准评测报告
- [ ] 首批 50,000 候选的 EIG 排序(全蛋白组 500 万池)
- [ ] 🔴 KN-8:下单 50k 寡核苷酸池(用户)

### Year 1 Q4(Months 10-12):年度答辩 + 文库接收
- [ ] 正式开题报告材料(D-BSSE 委员会)
- [ ] 慢病毒质粒构建启动(用户湿实验)
