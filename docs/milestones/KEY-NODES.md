# 关键节点与里程碑追踪

> 项目启动:2026-09-20(计划书 Year 1, Month 1)。计算侧由 Claude Code 自动推进;
> 标注 🔴 的节点 = **需要用户确认/操作**,到点停下等确认。
>
> **执行范围(2026-09-21 用户定):本会话只执行干实验。** 湿实验技术路线
> **完整保留**(见下文"技术路线"表,简化表述),由用户/实验室另行执行,
> 不在本会话执行范围内。干-湿闭环在本会话内以**主动学习模拟器**
> (公开数据验证)形式实现。

## 当前状态

- [x] KN-0:任务范围确认(执行研究计划)+ 起点确认(从零开始)— 2026-09-20 ✅
- [x] KN-1:计算环境搭建 ✅ 2026-09-21 — torch 2.9.1+cu126 GPU 冒烟通过(GTX 1050 Ti)
- [x] KN-2:文献综述 v0 ✅(3 篇综述 + 3 个 bib;🔴 待用户审阅)
- [x] KN-3:VDJdb ✅(199,289 条)+ IEDB ✅(2,302,095 条,1,394,080 肽,544 等位基因)— 2026-09-21
- [x] KN-4:EGNN 原型首训 ✅ 2026-09-21 — 诱饵判别 AUROC 0.771 / AUPRC 0.556(PDB 分组切分,无泄漏);同任务对照谱系:纯序列 0.556 < 序列+图统计 0.649 < EGNN 0.771
- [ ] KN-5:AlphaFold Server API(🔴 需用户提供)+ 结构批量预测 — 进行中:AF3 网页版已收 6 批;重提感知导入(链对匹配)+ 晶体链指纹回贴(`relabel_instances_by_chains.py`)+ 固定实例集成打分 + 可靠性闸门(ipTM)就位;⚠️ 2026-09-23 修复 pdb 标签错误(48/73,kn5 tcr_a/tcr_b 列写反);**2026-09-23 真实 native DA:危机解除**——天然 AF3 正样本 63→89→**144**(`af3_native/manifest.tsv`),`native_manifest.py` 链指纹累积匹配;3-seed val AUROC **0.959/0.952/0.964**;候选判决 mean pairwise Spearman **0.447(代理)→ 0.479(89)→ 0.729(144)=STABLE**、seed 过滤后全 73 实例 **0.802**(CI 0.72-0.86);——Rashomon 根因 = 训练数据规模(弱 seed 须过滤),已确认;排名声明可升级回可报告。剩余天然 29 + 交叉反应候选 358,🔴 待用户提交。**2026-10-12 更新**:native 台账 172;候选面板 21→154(5-seed 判决 STABLE r=0.609,CI 0.529–0.677,组级分化 5/17 组内 r<0.5);AF3 待传队列已**全部入批**——batch1 60(已传)/ batch2 49 / batch3 60 / batch4 47(收尾,内容唯一 job 清零;判重建包脚本 `scripts/build_af3_next_batch.py`),🔴 待用户提交 batch2/3/4
- [x] KN-6:v0 ✅ 2026-09-21 — 硬负样本三重筛选模块 + ProtoNet 少样本 AUROC 0.635(留出表位,超序列基线 0.609)
- [x] KN-7:v0 ✅ 2026-09-21 — 主动学习模拟器:EIG+ε-greedy 在 6% 采样达 2.0× 阳性召回增益
- [x] KN-7+:AL 模拟器升级 ✅ 2026-09-23 — 预筛肽池模式 + 批大小消融框架上线;原理性发现:①VDJdb 无法模拟预筛密度增益(阴性=已呈递肽,密度恒 7.0%)②EIG 增益阈值依赖:0.9 严格预筛下增益消失(1.0×)、0.5 宽松预筛下增益仍在(7.4% 处 1.76×)→ 增益是预筛阈值的函数,结构模型 EIG 是必要条件;③批消融:大批(3600×2)24.5% > 小批(300×20)21.6%,与文献"小批多轮"相反,批量存在下界
- [x] KN-4+:VDJdb 标签 ↔ 结构映射 ✅ 2026-09-22 — 281/291 复合物映射成功,2,309 条 CDR3 命中(791 同表位 / 1,517 异表位=交叉反应性证据);`structure_vdjdb_map.tsv`(含 1ao7 A6 TCR 的 8 肽交叉集,与 AF3 排名实验互证)
- [x] KN-8 计算侧:500 万池 EIG 排序 → 前 5 万候选清单 ✅ 2026-09-22(干跑)— 管线全通(4,654,973 唯一肽池 + Kimmtrak/A6 两份清单 + Uniprot 定位);诚实负面:弱模型 EIG 零富集(0/24),真实订购需域适应 EGNN 重排;🔴 实际订购(湿实验)仍由用户决定
- [ ] KN-11:in silico 临床 TCR 深扫描 — 序列空间部分已完成(= KN-8 清单,阴性结果);**结构空间首个可报告结果(2026-09-23)**:144-native 稳定 DA 模型下,两个致死性 off-target(MAGE-A12 模拟 KVAKELVHFL、titin ESDPIVAQY)仅比靶肽低 0.02/0.06(判为"与靶同水平",不 clear),`clinical_scan.py` → `kn11_clinical_scan.tsv`;案例 3(Kimmtrak)提交批已备(`af3_scan_kimmtrak/`),案例 2 无晶体结构;**2026-10-12 复扫**(172-native 5-seed 集成 × 154 候选面板):titin 模拟 ESDPIVAQY **0.948(反超靶肽,家族内 #1)**、靶肽 EVDPIGHLY 0.892(#2)、MAGE-A12 模拟 KVAKELVHFL 0.842(#3)——两个致死脱靶仍"与靶同水平";ILAKFLHWL 0.676 仍部分漏检;Kimmtrak 对照 YLEPGPVTA 0.283 正确低分——**结论在扩大 7 倍的面板上复现**
- [x] **2026-10-12 会话**:①AF3 next60 batch1 全部提交完毕,归档 `af3_next60_batch1/`;复盘发现**按 job_id 判重漏判**:60 个中 24 个的 5 链 FASTA 与更早任务逐字节相同(AF3 只看序列,提交即重复);②判重口径改为**内容哈希**,`already_uploaded.tsv` 升级为内容判重台账(现 391 条 done);③新提交包 `af3_next60/`(batch2,60 个):6 个临床表位家族变体矩阵补全 32 + MAGE/p53/EBV 等临床单点 28,md5 级唯一;④AutoDL 3080 Ti 跑通 KN-4 重训管线(STCRDab 625/634 PDB,1269 graphs,30ep:val displaced AUROC 0.937/AUPRC 0.924,graft 0.731/0.691);⑤**训练/打分归一化对齐修复**(`collate_per_graph`);⑥**重大数据事故与修复**:上一会话把 `af3_predictions/` 里 21 个 job(1ao7/1qrn/1tcr/2vlj/5brz 系)的 CIF 覆盖成同一内容,导致打分全部撞车(5brz 四条肽同分 0.8328)——已从 git 历史/隔离区恢复原文件,全库 md5 级复查干净(af3_predictions 574 / natives 856 / instances 全无重复);⑦**172-native 5-seed DA 重训(修复后,3080 Ti)**:val AUROC 0.949/0.952/0.959/0.947/0.954,val 两两 Spearman 0.943-0.970;**跨实例稳定性 STABLE:mean r=0.609(CI 0.529-0.677,n=154)**,各 TCR 组相关性全部有效(4jff 0.92 / 2gj6 0.81 / 3uts 0.81 / 5brz 0.60);KN-11 临床扫描(修复后):MAGE-A3 靶肽 EVDPIGHLY 0.892、titin 模拟 ESDPIVAQY **0.948(超靶肽,FATAL)**、MAGE-A12 模拟 KVAKELVHFL 0.842(FATAL)、ILAKFLHWL 0.676;Kimmtrak 靶肽 0.283(安全)——与 9 月结论一致;⑧用户又提交 30 个 batch2 job(folds_2026_10_12_07_22)已导入打分,batch2 现有 11 个有结果;⑨**AF3 待传队列收尾**:347 个待传全部入批(batch1-4,batch4=47 个 1E6 变体矩阵,内容唯一 job 清零),判重建包固化为 `scripts/build_af3_next_batch.py`(台账∪磁盘∪已入批∪FASTA md5 四闸),台账扩至 393 条并补登 2 个内容级漏判;benchmarks.md 补齐 2026-10-12 会话全记录(事故修复/collate 归一化/172-native/0.609 判决/临床 FATAL 复扫)
- 项目仓库:`/home/administrator/docs/meta-acdc/`(2026-09-22 起由 opencode 接续推进)

### 文献综述关键发现(2026-09-20)
1. ACDC 原文确认:bioRxiv 2025.07.23.666264(Hong K.-L. 等 14 人)。**缩写展开与计划书不一致**:原文是 "antigen-presenting cell detecting cytokine",计划书写 "Artificial Cell Display of pMHC" → 需与实验室核对
2. NetMHCpan 最新是 **4.2**(2025,Front Immunol),计划书引用的 4.1 已过时 → 预筛管线改用 4.2
3. 直接方法学基线:**MCGLPPI++**(JCIM 2025,几何 GNN 少样本 PPI,含 TCR-pMHC 簇验证)、**ZeroBind**(MAML++)
4. vanilla MAML 有任务冲突/过度记忆陷阱 → 实现 MAML++ / ProtoNet
5. 主动学习文献规则"批数多、每批少、收敛快" → 模拟器需做批大小消融

## 执行计划总表(本会话干实验)

| # | 节点 | 计划时间 | 谁操作 | 通过标准 |
|---|------|---------|--------|---------|
| KN-1 | 🔴 计算环境安装 | 2026-09 | 用户批准安装 | torch GPU 冒烟测试 ✅ 已完成 |
| KN-2 | 文献综述 + 参考文献库 | Month 1-3 | AI | 3 篇综述 + bibtex ✅ 已完成 |
| KN-3 | 数据管线(VDJdb/IEDB/结构集) | Month 1-3 | AI | 标准数据集 + 统计报告 ✅ 已完成 |
| KN-4 | EGNN 原型 + 基准 | Month 4-6 | AI | 诱饵任务 + 对照谱系 ✅ 已完成(v9 侧链特征进行中) |
| KN-5 | 🔴 结构批量预测(AF3 网页版手动提交;TCRmodel2 路线已废弃) | Month 4-6 | 用户协助提交 | 首批 ~100 条 TCR-pMHC 结构入库 |
| KN-6 | 硬负样本挖掘 + 元学习基线 | Month 4-6 | AI | ✅ 已完成(ProtoNet 0.635 / MAML++ 0.625,序列天花板 0.62-0.64) |
| KN-7 | 主动学习模拟器 | Month 6-9 | AI | ✅ v0 完成(6% 采样 2× 增益);升级:结构打分 + 采集函数对比 |
| KN-11 | in silico 临床 TCR 深扫描(Kimmtrak 等,公开数据) | Month 25-27 | AI | off-target 风险排序清单 |
| KN-13 | 风险看板 v1(Web"安全雷达") | Month 28-30 | AI | 部署可用 |
| KN-14 | 论文 1:Meta-TCR-GNN(含 Cα 分辨率诊断) | Month 22-24 | AI 起草/用户投稿 | Nature MI / Bioinformatics |
| KN-15 | 论文 2:主动学习模拟器 + in silico 安全评估 | Month 34-36 | AI 起草/用户投稿 | Bioinformatics / Nat Methods 级 |
| KN-16 | 交付物打包:软件套件 + Docker + 数据集 | Month 37-39 | AI | 开源发布 |

## 技术路线(计划书完整保留,简化表述;湿实验由实验室另行执行)

| 环节 | 内容 | 本会话状态 |
|------|------|-----------|
| ACDC 平台实验 | 5 万肽寡核苷酸池合成 → 慢病毒递呈 → FACS 分选 → NGS 富集分析 | 不执行(技术路线保留) |
| 干-湿闭环 3-5 轮 | AI 挑库 → ACDC 验证 → 模型更新 → 盲区定向采样 | 以模拟器形式实现 ✅ |
| 金标准验证 | 原代 T 细胞杀伤(IncuCyte/ELISPOT)、NSG 小鼠体内 | 不执行(技术路线保留) |
| 安全看板 | Web 交互"安全雷达" | ✅ 在计划内(KN-13) |

## 各阶段执行清单

### Year 1 Q1(2026-09 ~ 2026-11,Months 1-3):整合与文献综述
- [ ] 文献综述:平台(进行中)/ ML 模型(进行中)/ 数据库与主动学习(进行中)
- [ ] 环境搭建(KN-1)
- [ ] 数据管线 v0(KN-3)
- [ ] ACDC 原始数据盘点:向实验室要历史数据清单(🔴 需用户对接实验室)

### Year 1 Q2(Months 4-6):Meta-TCR-GNN 原型
- [x] EGNN 骨架实现 + 接口图构建(10Å 截断)
- [x] AlphaFold 3 结构预测管线(KN-5)— 进行中(用户提交)
- [x] 硬负样本三重筛选(Levenshtein ≤3 / 静电势 Pearson >0.8 / 不确定性 0.4-0.6)
- [x] 元学习框架(MAML / Prototypical)
- [x] 基准评测(VDJdb、IEDB;对照 NetMHCpan、ERGO、NetTCR)— 见
  [`benchmark-report.md`](../benchmark-report.md)(文献数字跨数据集引用)

### Year 1 Q3(Months 7-9):基准测试 + 首批文库下单
- [x] 公开数据基准评测报告 — [`benchmark-report.md`](../benchmark-report.md)
- [x] 首批 50,000 候选的 EIG 排序(全蛋白组 500 万池)— 干跑完成
  (诚实负面:弱代理零富集,待结构重排)
- [ ] 🔴 KN-8:下单 50k 寡核苷酸池(用户)

### Year 1 Q4(Months 10-12):年度答辩 + 文库接收
- [ ] 正式开题报告材料(委员会)
- [ ] 慢病毒质粒构建启动(用户湿实验)
