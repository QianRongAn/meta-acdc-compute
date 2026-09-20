# TCR 特异性预测 ML 与几何深度学习综述

> Meta-ACDC 项目文献综述（撰写日期 2026-09-20）。所有数字经网络检索核对；来源冲突或未能独立核实的条目以 **[未核实]** 标注。完整引用见同目录 `references-ml.bib`。

## 1. TCR–pMHC 结合预测模型（序列模型）

- **ERGO**（Springer 等, *Front. Immunol.* 11:1803, 2020）。自编码器(AE)/LSTM 双变体，输入 TCRβ CDR3 + 肽；训练于 McPAS-TCR 与 VDJdb。三任务基准：TPP-I（已知肽-TCR 配对判别）AUC 0.860/0.840（McPAS/VDJdb），TPP-II 0.810/0.792，TPP-III（肽与 TCR 均未见）仅 0.669。代码 `github.com/louzounlab/ERGO`。**与 Meta-TCR-GNN 的关系**：双编码器+预训练自编码器路线；TPP-III 表明"全新对"场景仅靠序列编码上限很低。
- **ERGO-II**（Springer 等, *Front. Immunol.* 12:664514, 2021）。加入 TCRα CDR3、V/J 基因、MHC 等位基因、T 细胞类型。特征消融：β CDR3 贡献最大、MHC 最小；流感肽 GILGFVFTL 的单肽 AUC 从 0.691（仅 β）升至 0.876（全特征）。代码 `github.com/IdoSpringer/ERGO-II`。**关系**：指导输入特征优先级（成对 αβ + V/J + MHC 结构信息应保留）。
- **NetTCR-1.0**（Jurtz 等, bioRxiv 433706, 2018）。1D CNN，输入肽 + CDR3β（限 HLA-A*02:01）。训练：IEDB 9,015 条 + MIRA 379 条；负样本=错配 + eluted self-peptides + 20 万条健康人 CDR3。>90% 序列同一性的 TCR 同分区防泄漏。代码 `github.com/mniellab/netTCR`。**关系**：负样本构造与按序列相似度分区的先例。
- **NetTCR-2.0**（Montemurro 等, *Commun. Biol.* 4, 2021, DOI 10.1038/s42003-021-02610-3）。成对 CDR3α+β 浅 CNN；配对链显著优于单链（p<0.05 bootstrap）；外部基准 AUC 0.804（McPAS）/0.891（VDJdb）；肽排序任务约 75% 命中（随机 33%）；每肽需 ~150 个独特阳性 TCR 才能 AUC>0.75。**关系**：成对链信息的价值、数据规模对性能的硬约束。
- **NetTCR-2.2**（Jensen & Nielsen, *eLife* 12:RP93934, 2024; bioRxiv 2023.10.12.562001）。pan-specific + peptide-specific 混合训练、dropout 0.6、loss-scaling、TCRbase 相似性集成；训练集 10,239 条/435 肽（IEDB+VDJdb 成对数据）。仅 15 个阳性 TCR 即可达可用性能；IMMREP22 基准上 SOTA；作者指出 IMMREP 基准存在数据冗余与负样本泄漏。代码 `github.com/mnielLab/NetTCR-2.2`。**关系**：少样本任务划分（peptide-specific head）正是 MAML 式元学习的朴素版本。
- **pMTnet**（Lu 等, *Nat. Mach. Intell.* 3:864–875, 2021, DOI 10.1038/s42256-021-00383-2）。三输入（CDR3β、肽、MHC I 等位基因）+ Atchley 因子编码 + 迁移学习。训练 ~30,801 三元组；独立测试 619 对（>30 篇研究）：**AUROC 0.827、AUPRC 0.565**。**关系**：将 AUPRC 引入 TCR 结合预测评估的先例；三要素输入模板。
- **pMTnet-omni**（Han 等, bioRxiv 10.1101/2023.12.01.569599; *Nat. Commun.* 2026, DOI 10.1038/s41467-026-73396-3）。MHC 用 ESM2/ESMFold 结构编码，肽用 Atchley，TCR 用 Vα/Vβ/CDR3α/CDR3β 四编码器（注意力自编码器）；pan-MHC、跨人鼠、I/II 类。训练 111,202 对（69 个数据集）；验证 **AUROC 0.888 / AUPRC 0.741**，比 pMTnet v1 高 ~0.1；掩码 β 链 AUROC 降至 0.685（α 链 0.848）；前瞻患者数据（BEAM-T）AUROC 0.77–0.85；LiL 策略仅用 5 个阳性即实现 PRAME TCR 亲和力成熟预测。**关系**：结构嵌入（ESM2）融入序列模型 + 前瞻外部验证范本。
- **TITAN**（Weber 等, *Bioinformatics* 37(S1):i237–i244, 2021, DOI 10.1093/bioinformatics/btab294; arXiv:2105.03323）。双模态注意力网络（CDR3β + 表位，可选 SMILES 原子级编码）；在 VDJdb 10,599 条（87 表位）+ ImmuneCODE COVID 数据上训练（共 46,290 条、192 表位）；负样本=打乱配对。未见 TCR 10-fold AUC **0.87±0.005**（优于 ImRex 0.61）；严格未见表位 CV 仅 **0.62**，McPAS 来源独立测试 0.78。代码 `github.com/PaccMann/TITAN`。**关系**：未见表位泛化（~0.6）是全部序列模型的天花板——本项目元学习要攻克的目标。
- **ATM-TCR**（Cai 等, *Front. Immunol.* 13:893247, 2022）。双多头自注意力编码器 + MLP 解码；128,142 对（931 表位，VDJdb/McPAS/IEDB）。比 NetTCR/ERGO AUC 高约 2%、recall 高 6.3%；注意力图置信度使 SARS-CoV-2 外推准确率提升 41.11%、recall 25%。**关系**：注意力图作为置信度/主动学习信号。
- **PanPep**（Gao 等, *Nat. Mach. Intell.* 5:236–249, 2023, DOI 10.1038/s42256-023-00619-3）。**MAML 元学习 + 神经图灵机(NTM)**；输入肽 + CDR3β。数据 699 肽/29,467 TCR/32,080 对（IEDB+VDJdb+PIRD+McPAS）。few-shot / zero-shot（PLGS 解耦蒸馏）/ majority 三设置全面优于 pMTnet、ERGO-II、DLpTCR；零样本识别外源肽。代码 `github.com/bm2-lab/PanPep`。**关系**：本项目 Meta-TCR-GNN 最直接对标物——但其仅用 CDR3β 序列、无结构、无主动学习，是我们的差异化空间。
- **TULIP**（Meynard-Piganeau 等, *PNAS* 121(24):e2316401121, 2024）。ProtT5 系无监督 Transformer，**只用正样本训练**。已知表位优于 NetTCR-2.0；对未见表位的泛化随编辑距离衰减（≤4 仍有效）；对 UURA/HRS 负采样方式鲁棒（PanPep 在 UURA 下退化为随机）；与 DMS 实验 EC50 的 Spearman 相关最高 0.47。代码 `github.com/barthelemymp/TULIP-TCR`。**关系**：直接证明"随机负样本偏差"可制造虚假性能——我们的负采样必须按 UURA/EN 双方案并报告两套结果。
- **UniTCR**（Gao 等, *Cell Genomics* 4:100553, 2024）。TCR（Atchley+transformer）与 T 细胞转录组双模态对比学习 + 单模态保持；下游表位结合分类模块（LSTM+交叉注意力）在 majority/few-shot/zero-shot 均优于 KNN、ERGO、TITAN、pMTnet。代码 `github.com/bm2-lab/UniTCR`。**关系**：表示学习 + 少样本评估范式。
- **SCEPTR**（Nagano 等, *Cell Systems* 2025; arXiv:2406.06397）。六条 CDR 环输入 + 自对比学习 + MLM 预训练；下游性能超过 ESM/ProtTrans/TCR-BERT 等更大模型、与 TCRdist 相当；对比学习有效对抗 VDJ 重组分布偏置。代码 `github.com/yutanagano/sceptr`。**关系**：小模型 + 对比学习路线可作为序列编码器。
- **EPACT**（Zhang 等, *Nat. Mach. Intell.* 6(11):1344–1358, 2024, DOI 10.1038/s42256-024-00913-8）。表位锚定对比迁移学习（peptide LM + TCR LM + pMHC 模型）；优于 ERGO-II、NetTCR-2.2、STAPLER；对编辑距离≥5 的表位保持 ~0.7 AUC；用 TCR-pMHC 3D 结构（PDB）微调后预测 CDR–表位残基级距离/接触（与 TEIM-Res 相当）。代码 Zenodo DOI 10.5281/zenodo.10996144。**关系**：证明"结构监督 + 残基级接触预测"可行——可直接改造为 EGNN 辅助任务。
- **TCRoss**（Yang 等, *Brief. Bioinform.* 26(6):bbaf609, 2025）。性质交互 Transformer（氨基酸物化性质 + Linformer 线性注意力），模型不直接读取原始序列以防记忆偏置；配套 **TCRC-200k** 数据集（~218k 条，~110k 高质量阳性；RN/EN 双负采样方案）与 EES 训练策略；湿实验 T 细胞活化验证。**关系**：负采样工程化方案 + 数据集候选。
- **SCORPIO**：多轮检索（含 arXiv/bioRxiv/期刊）**未能定位到该名称的 TCR 预测模型**，疑为误传或未公开方法；相近命名 "SCORPION" 亦无 TCR 相关文献。[未核实]
- **其他 2024–2026**：tcrLM（大语言模型，预训练 2.27B 残基，报告 AUC 0.937/0.933，期刊未核实）[未核实]；PredicTCR（Tan 等, *Nat. Biotechnol.* 42:134–142, 2024）：XGBoost 从单细胞 RNA-seq 抗原无关地预测肿瘤反应性 TCR（几何均值 0.38→0.74）——邻近任务但方法学可参考。

## 2. 几何深度学习（蛋白质）

- **EGNN**（Satorras, Hoogeboom & Welling, ICML 2021, PMLR 139:9323–9332; arXiv:2102.09844）。核心 EGCL 更新：`m_ij = φ_e(h_i, h_j, ‖x_i−x_j‖², a_ij)`；坐标更新 `x_i^(l+1) = x_i^l + C·Σ_{j≠i}(x_i^l−x_j^l)·φ_x(m_ij)`（C=1/(M−1)，实现中常把位移除以 `‖x_i−x_j‖+1` 增强数值稳定）；节点特征 `h_i^(l+1)=φ_h(h_i, Σ m_ij)`。对旋转/平移/反射/置换等变（证明见附录 A）；只需 0/1 阶特征、无需球谐。**实现要点（PyTorch）**：社区实现 `github.com/lucidrains/egnn-pytorch`（torch 原生，依赖少）；批量图可用 `pyg` 或手工 scatter；坐标作为张量前向传播需注意 `.clone().detach()` 防梯度爆炸与质心平移；速度向量变体（式 7）可用于构象采样。**关系**：Meta-TCR-GNN 骨干；评分任务仅需 E(3)-不变输出（图级池化）。
- **GVP-GNN**（Jing 等, ICLR 2021 Spotlight; arXiv:2009.01411）。几何向量感知器：标量+向量双通道、向量 L2 范数并入标量通路、证明 O(3) 不变/等变与通用逼近。蛋白表示为 Cα k=30 近邻图；节点特征=前向/反向 Cα 单位向量 + 推算 Cβ 方向 + 二面角，边特征=方向向量 + 距离 RBF + 序列间隔正弦编码。CATH 4.2 计算蛋白设计：perplexity 5.29、序列回收率 40.2%（优于结构 GNN 6.55/37.3%）。代码 `github.com/drorlab/gvp`。**关系**：其节点/边几何特征构造可作为 EGNN 输入特征模板。
- **GearNet**（Zhang 等, ICLR 2023; arXiv:2203.06125）。关系图消息传递 + 稀疏边缘消息传递（仿 Evoformer 三角注意力）；五种几何自监督预训练任务（multiview contrast、残基类型/距离/角度/二面角预测）；multiview contrast 得 EC Fmax 0.874、GO-BP 0.490；百万级样本预训练达到十亿级序列模型水平。代码 `github.com/DeepGraphLearning/GearNet`（TorchDrug/TorchProtein）。**关系**：多视图对比预训练可作为 EGNN 初始化方案。
- **EquiDock**（Ganea 等, ICLR 2022 Spotlight; arXiv:2111.07786）。刚性蛋白-蛋白对接：成对独立 SE(3) 等变图匹配网络 + 最优传输关键点选择 + 可微 Kabsch；DIPS（~42k 对）预训练 + DB5/DB5.5 微调；比传统对接快 80–500×。代码 `github.com/octavian-ganea/equidock_public`。**关系**：若需"快速生成 TCR-pMHC 对接构象"作为候选构象打分/主动学习批次采样。
- **DiffDock**（Corso 等, ICLR 2023; arXiv:2210.01776）。配体姿态扩散生成模型（平移/旋转/扭转分量上扩散，SE(3) 等变分数）；PDBBind top-1 成功率 38%（RMSD<2Å），vs 传统方法 23%、先前 DL 方法 20%；apo 结构上保持 21.7%。**关系**：构象集成/生成式评分的备选范式；与 EGNN 可组合（EGNN 作等变去噪网络）。
- **综述**：Atz, Grisoni & Schneider, *Nat. Mach. Intell.* 3(12):1023–1032, 2021（"Geometric deep learning on molecular representations"）——对称性、表示层级与药物设计应用的系统梳理。**关系**：为等变/不变设计选择提供理论框架。
- **实践要点**：SE(3)/E(3) 等变网络内存占用与坐标通道数成正比，TCR-pMHC（~400 残基 × 5 链）完全可负担；不显式建模坐标时可退化为不变网络（帧平均或只输出标量）；训练中坐标应相对质心零均值化。

## 3. TCR–pMHC 结构预测

- **AlphaFold 3**（Abramson 等, *Nature* 630:493–500, 2024, DOI 10.1038/s41586-024-07487-w）。扩散式架构，覆盖蛋白、核酸、小分子、离子与修饰残基；抗体-抗原复合物精度显著优于 AF-Multimer v2.3（P=6.5×10⁻⁵，需 1,000 seeds 而非 5）；蛋白-蛋白 DockQ>0.23 成功率显著提升（P=1.8×10⁻¹⁸）；仍依赖 MSA 深度。**AlphaFold Server**：免费但仅限非商业使用；每日作业限额约 10–20 个（来源冲突：TAMU 资料称 20/天，第三方综述称 10/天，2024-07 曾扩容）[未核实精确值]；**输出条款禁止用于训练生物分子结构预测模型或自动化结合预测系统**——本项目若批量调用需走条款允许的科研公开用途边界，结构训练数据建议改用本地 TCRmodel2/ColabFold。
- **AlphaFold-Multimer**（Evans 等, bioRxiv 10.1101/2021.10.04.463034）。多聚体特化 AF2：4,433 复合物上 67%/69% 异源/同源界面预测成功；训练按 384 残基裁剪（内存约束）；对无模板的 TCR-pMHC/抗体对接近失败（作者自述抗体差；MSA 中 CDR 区噪声大）。
- **TCR-pMHC 结构预测现状（关键证据）**：Bradley (*eLife* 12:e82813, 2023) 证明 AF-Multimer 在 130 个 TCR-pMHC 上质量高度可变、常见错误对接；其 AF_TCR 混合模板管线（12 种对接几何 × 每链独立模板 + model_2_ptm + 93 结构两轮微调）decoy AUC **0.82**（YLQ/ELA 系统 0.97）。第三方基准（bioRxiv 2025.01.12.632367，STCRDab-22）：**AF3 平均 DockQ 0.499**、TCRmodel2 0.566、AF-Multimer 0.490、TCRdock 0.143、HADDOCK 0.098；70 复合物基准（*Brief. Bioinform.* 2025）MSA 类方法（AF3）优于 PLM 类；CDR3 pLDDT 可用于重排序；ipTM>0.85 一般可靠、≤0.65 常伴随大误差；评审共识 "AF3 cannot reliably predict TCR docking"。**结论**：现成结构预测不足以直接作交叉反应判据——需要学习型（几何深度学习）评分/精修，这正是本项目定位。
- **ColabFold**（Mirdita 等, *Nat. Methods* 19:679–682, 2022, DOI 10.1038/s41592-022-01488-1）。MMseqs2 服务器式 MSA（UniRef100/PDB70/环境库）+ AF2/RoseTTAFold；单 GPU 日均可产近 1,000 个结构；CASP14 自由建模 TM-score 0.826。**关系**：批量生成 TCR-pMHC 结构数据的主力工具。
- **TCRmodel2**（Yin 等, *NAR* 51(W1):W569–W576, 2023, DOI 10.1093/nar/gkad356）。AF 2.2/2.3 改造：TCR/MHC 限定 MSA 库（52,096/137,991 条）、pMHC 复合物模板（884 I 类 + 44 II 类）、模板直接用输入序列；~15 分钟/复合物；置信度（ipTM+pTM 组合）与 DockQ 相关 r=0.75（≥0.85 可靠，≤0.49 不可靠）。服务器 `tcrmodel.ibbr.umd.edu`。**关系**：本地可部署的开源结构生成器（优于受条款限制的 AF Server）。
- **ImmuneBuilder**（Abanades 等, *Commun. Biol.* 6:575, 2023, DOI 10.1038/s42003-023-04927-7）。ABodyBuilder2 / NanoBodyBuilder2 / TCRBuilder2：抗体 CDR-H3 基准 RMSD 2.81Å（AF-Multimer 2.90Å），纳米抗体 2.89Å，GPU 上 ~5 s；4 结构集成给出逐残基误差估计（~1Å 阈值可把 CDR-H3 RMSD 4.46→2.53Å）。代码 `github.com/oxpig/ImmuneBuilder`。**关系**：快速 TCR 结构先验 + 不确定度估计（供主动学习）。

## 4. 临床 TCR 交叉反应案例（安全设计约束）

- **Morgan 2013**（*J. Immunother.* 36(2):133–151, DOI 10.1097/CJI.0b013e3182829903）：HLA-A*02 限制性 MAGE-A3（KVAELVHFL，A2 转基因鼠 TCR + A118T 亲和力增强）；9 例患者（7 黑色素瘤、1 滑膜肉瘤、1 食管癌），5 例缓解；3 例神经毒性、2 例死亡（坏死性白质脑病）；死因为 **MAGE-A12 在脑内低水平表达**——MAGE 家族内序列同源交叉反应（亦识别 MAGE-A3/A9/A12，较弱 A2/A6）。
- **Cameron 2013**（*Sci. Transl. Med.* 5(197):197ra103, DOI 10.1126/scitranslmed.3006034）：HLA-A*01 限制性 MAGE-A3（EVDPIGHLY），噬菌体展示亲和成熟 a3a TCR；2 例患者在输注后 5 天内死于急性心衰；事后丙氨酸扫描锁定 **titin 肽 ESDPIVAQY**（与 EVDPIGHLY 9 位中 4 位不同）；同一 TCR 还识别 MAGE-A6、MAGE-B18 肽。教训：体外常规筛选未能预见。
- **Linette 2013**（*Blood* 122(6):863–871, DOI 10.1182/blood-2013-03-490565）：亲和力增强 A3A TCR（HLA-A*01/MAGE-A3）在骨髓瘤与黑色素瘤两项试验中各 2 例死亡（心源性休克，心肌酶大幅升高）；心脏组织无 MAGE-A3 表达；**iPSC 分化跳动心肌细胞杀伤实验**证实 titin 交叉识别。教训：功能筛选需组织特异细胞模型。
- **Raman 2016**（*Sci. Rep.* 6:18851, DOI 10.1038/srep18851）：解析 MAG-IC3 TCR 与 A1-MAGE-A3（2.62Å）及 A1-titin（2.40Å）复合物晶体结构；**两条肽骨架 Cα RMSD 仅 0.285Å、整体复合物 RMSD 0.418Å、TCR 交叉角同为 57°、CDR 环构象一致**——"直接分子模拟（direct molecular mimicry）"机制；据此理性设计 TCR 突变（如 R31β）提高靶/脱靶区分度。**对本项目的含义**：交叉反应本质是 pMHC 表面形状与化学互补性的重叠，天然适合用等变几何网络建模。
- **Tebentafusp/Kimmtrak**（2022-01-26 FDA 批准；Nathan 等, *NEJM* 385:1196–1206, 2021）：高亲和可溶 TCR（gp100:280–288 YLEPGPVTA / HLA-A*02:01）+ 抗 CD3 scFv 双特异分子。Ph3 IMCgp100-202（n=378, 2:1）：中位 OS **21.7 vs 16.0 个月（HR 0.51, 95%CI 0.37–0.71, P<0.0001）**；3 年 HR 0.68、5 年 HR 0.67（2026 AACR）；ORR 9% vs 5%。安全性：皮疹/发热/低血压为主（gp100 在黑素细胞正常表达所致的 on-target/off-tumor 毒性，皮肤自限），≥3 级 CRS <1%，黑框警告 CRS；无治疗相关死亡。**对比教训**：商业化 TCR 疗法靠**剂量递增 + 脱靶组织毒性可管理**获得批准，而 MAGE-A3 亲和力增强 TCR 的致死性来自**未预见的表位交叉反应**——即本项目要预测的风险类型。
- **机制小结**：(i) 家族内同源交叉（MAGE-A3/A6/A9/A12/B18，序列相似）；(ii) 非同源分子的结构模拟（titin，序列不同但骨架/表面化学兼容）；(iii) on-target/off-tumor（gp100–黑素细胞）。第 (ii) 类无法用序列比对捕获，只能靠结构/几何学习。

## 5. 对模型设计的启示（Meta-TCR-GNN 具体决策）

1. **输入特征**：成对 CDR3α/β（证据：NetTCR-2.0、pMTnet-omni β 链主导）+ V/J 基因 + 肽 + MHC；结构侧采用 GVP 风格节点特征（Cα 前向/反向单位向量、Cβ 方向、二面角）与边特征（距离 RBF + 方向 + 序列间隔）喂 EGNN；序列侧可选 SCEPTR/ESM-2 嵌入做初始化（注意内存与延迟）。
2. **负采样**：弃用纯随机错配；采用 UURA + EN 混合（TULIP 教训、TCRC-200k 方案），并按编辑距离/结构相似度构造"难负样本"（相近 pMHC 的错配）；**两套负采样下分别报告结果**，防止虚假高性能。
3. **评价指标**：以 **AUPRC 为主**（1:10–1:100 类不平衡下 AUROC 失真），辅以 AUROC、AUROC@0.1（低假阳性区，NetTCR-2.2 的教训）；按表位分层报告 seen/unseen（编辑距离分桶）性能；保留前瞻外部验证（pMTnet-omni BEAM-T 范式）。
4. **基准**：IMMREP22 数据集（17 表位、VDJdb 四聚体/十聚体分选、5:1 负正比）+ IEDB/VDJdb 标准 TCR split 与 epitope split；少样本设定对齐 NetTCR-2.2（≥15 阳性/肽）与 PanPep（<10 阳性 few-shot）；注意 IMMREP 已知泄漏问题，需自建去冗余验证集。
5. **元学习设计**：以"表位（或结构聚类后的 pMHC 构象簇）"为任务（PanPep 式），MAML 内环适配 + 外环跨任务初始化；与 PanPep 的差异点=加入 3D 结构模态与主动学习回路；NTM 记忆模块可选做 zero-shot 扩展。
6. **结构数据来源**：训练结构用本地 TCRmodel2/ColabFold/ImmuneBuilder 批量生成（AF Server 条款限制）；残基级距离/接触矩阵（TEIM-Res、EPACT 微调范式）作为 EGNN 辅助监督，直接呼应 Raman 2016 的接触位点机制。
7. **等变性与实现**：EGNN 输出图级不变标量做二分类；坐标更新层数 3–6 即可；质心零均值化、位移归一化（除以 ‖x_i−x_j‖+1）、梯度裁剪；必要时用速度向量变体做构象扰动集成（测试时增强）。
8. **主动学习闭环**：候选表位按模型不确定性（注意力图/预测方差/结构置信度 ipTM≤0.65 过滤）+ 结构可及性排序，优先湿实验验证高风险交叉反应对（对照案例：titin ESDPIVAQY 类）——合成免疫细胞（synthetic immune cells）即为此提供可控负对照与验证平台。
