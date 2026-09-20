# 数据库、基准与主动学习/元学习方法综述

> Year 1 Q1 文献综述 · 2026-09-20 · 与计划书 Section IV Module 2 / VI 对应。

## 1. 数据库与基准(含本项目实测)

### VDJdb(实测,2026-06-03 版)
- 下载:GitHub 镜像 antigenomics/vdjdb-db release(41 MB zip,含 vdjdb.txt 等 11 文件)
- 本次清洗(管线 `src/meta_acdc/data/`):**199,289 条去重记录、178,974 个独特 CDR3-β、2,068 个表位肽、170 个 HLA 等位基因**
- 稀疏性证据:17.9 万 TCR 只覆盖 2,068 个表位 → 表位覆盖率极低,直接印证计划书"泛化缺口"
- 标签:vdjdb.score 0–3(0=无应答 … 3=强应答),全部 199,289 条带分;0 分条目可作天然阴性(需谨慎:多数未系统验证)

### IEDB
- mhc_ligand_full 单文件导出(CSV):MHC 配体结合/递呈数据,含 qualitative_measure(Positive/Negative)与 IC50 定量
- 用途:NetMHCpan 类**递呈预筛**的训练/基准数据(与激活标签区分——见 schema 注释)
- 下载进行中(`download_iedb.py`,legacy downloader URL)

### NetMHCpan:计划书引用 4.1,最新已是 4.2
- **NetMHCpan-4.2**(Nilsson et al., Front Immunol 2025, DOI 10.3389/fimmu.2025.1616113):基于 4.1 的 NNAlign_MA 框架 + 新 EL 数据集(41 个新集)+ 结构特征(删除位点的 BLOSUM50 组成、PDB 推导的 MHC 接触频率)+ IEDB 表位迁移学习
- 覆盖 **163 个 MHC 分子**(4.1 为 130),人群覆盖率 ~96%(A)/~93%(B)/~97.3%(C)
- 结论要点:迁移学习增益**有限**;EL 模型 > BA 模型(独立基准 bioRxiv 2025.04.10.648169 亦证实)
- 🔴 **行动项**:预筛管线应改用 4.2(计划书写的是 4.1);`>10^7 HLA 限制性肽` 的说法需按新工具复算

## 2. 主动学习:已有实验量压缩证据

| 工作 | 领域 | 关键数字 |
|---|---|---|
| ALDE(bioRxiv 2026) | 远红光荧光蛋白定向进化 | **5 轮仅 120 个变体**筛出 12 个改良,最佳 4.35× 亮度;**批数多、每批样本少 → 收敛更快**(通用设计规则) |
| ALSEBO(bioRxiv 2026) | GFP 设计(协同进化+BO) | 20^L 空间压缩到 2,143 候选,40–45 次评估命中 oracle 最优 |
| EVOLVEpro(Science 2025) | pLM + 少样本主动学习 | 6 个蛋白达 **100 倍**性能提升 |
| DeepDE | avGFP 进化 | 4 轮 74.3× 活性提升(对照:sfGFP 40.2× 用多年) |
| ML-DE 综述性结论 | — | 累计重训练 > 单轮;模型类别常非瓶颈;Spearman 0.098–0.286 → **迭代实验反馈 > 纯预测** |

**对干-湿闭环的启示**:
1. 计划书"0.05% 采样覆盖 100 倍效率"的方向与文献一致,但文献更倾向**更多小批次**(每轮 5 万条肽 vs 分成多批) —— 每轮 5 万是文库合成经济性决定的批量上限,可考虑**池内子批次**设计
2. ε-greedy(10–20% 随机探索)有文献支撑:纯贪心会陷入局部最优(计划书 Risk C)
3. 冷启动风险(Risk A)是真实痛点:ALDE 证明 **pLM 零样本先验与真实适应度相关性弱**,结构预筛(NetMHCpan/AF3 富集)比随机采样更稳妥——与计划书缓解方案一致
4. 噪声:文献确认实验噪声会经正反馈放大——计划书 Risk B 的 GCE 噪声鲁棒损失 + 校准锚点方案合理

## 3. 元学习:方法谱系与陷阱

- **ZeroBind**(Nat Commun 2023):MAML + MAML++ 稳定化(多步损失、逐层学习率),ESM-2 特征 + AlphaFold2 图,少/零样本 DTI 预测——我们的直接方法论前身
- **MCGLPPI++**(JCIM 2025):几何 GNN + 簇切分 + Lennard-Jones/Coulomb 显式能量建模的少样本 PPI 亲和力预测,**在 TCR-pMHC 测试簇上验证**——与 EGNN 设计最接近的已发表工作,🔴 **应作为我们的直接对照基线**
- **BioBridge**(Adv Sci 2025):归纳-联想元学习,5-shot BindingDB 全面超越 MAML++/原型网络/ANIL
- **RR-ADS**(J Med Chem 2025):报告 vanilla MAML 的**任务冲突**与**过度记忆**问题;解法:内部衰减器 + 元正则化 + 难度感知调整
- **MetaMBP**(JCIM 2025):深度度量元学习,少样本多标签肽功能预测
- 趋势综述(Adv Drug Deliv Rev 2025):少样本/元学习已是数据受限药物发现的标配策略

**对元学习模块的启示**:
1. 任务定义:计划书"一个 TCR 对多样肽的亲和力预测 = 一个 Task"与 ZeroBind/MCGLPPI++ 的簇级任务切分吻合;任务切分应**按 TCR/抗原簇**(CDR3 聚类或结构簇),避免随机切分导致的任务泄漏
2. 支持集构造(计划书 1:4 硬负样本)有先例可循;难度感知采样(decision boundary 0.4–0.6)对应 RR-ADS 的 difficulty-aware 设计
3. MAML 训练不稳定是已知问题 → 直接用 MAML++ 稳定化变体或原型网络 + 度量的组合,别裸写 vanilla MAML
4. 评估:报告 seen-TCR / unseen-epitope / unseen-both 三档(CV 切分),与 TITAN 的 0.87/0.62 模式可比

## 4. 汇总行动项

- [ ] 预筛工具从 NetMHCpan-4.1 → **4.2**
- [ ] 基线清单加入 MCGLPPI++、ZeroBind、ERGO-II、pMTnet-omni、PanPep、TITAN(见 02 综述)
- [ ] 主动学习模拟器先跑**批大小消融**(验证"多批小量"规则在 5 万/批约束下的最优子批次)
- [ ] 元学习框架优先实现 **MAML++ / ProtoNet**(避免 vanilla MAML 陷阱)
