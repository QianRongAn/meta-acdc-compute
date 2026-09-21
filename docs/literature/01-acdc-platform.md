# ACDC 平台与高通量 TCR-pMHC 筛选技术综述

> Year 1 Q1 文献综述 · 2026-09-20 · 与计划书 Section I/III 对应。
> ⚠️ 标记 [未验证] 的条目未逐一核对原始出处。

## 1. ACDC 平台原文(Hong et al., 2025, bioRxiv)

**出处**: Hong K.-L., Gao B., Stalder L., Ehling R.A., Horvath M., Wehrle S., Forster J.L., Junet V., Kucharczyk J., Lalevee S., Diringer M.-C., Yang Q., Vazquez-Lombardi R., R.S.T. *Predicting TCR antigen specificity at proteome-scale with synthetic immune cells and machine learning.* bioRxiv 2025.07.23.666264 v1(2025-07-27 发布)。

- **平台构造**(HEK293 逐步 CRISPR 工程):
  - STAT5 驱动的 mRuby2 荧光报告(感应共培养 T 细胞旁分泌 IL-2 信号)——以激活信号而非结合信号为读出
  - 泛 HLA-I 敲除(保守外显子 4 gRNA,消融 HLA-A/B/C 本底)
  - 单等位整合着陆垫(CCR5 位点,Cas9-GFP → BFP 双表型标记)
  - HLA-A\*02:01 敲入 + **18 个不同 HLA-I 等位基因**的系列细胞系
  - B2M 敲除后用 B2M-抗原融合构建体(T2A 自切、CMV 启动子、ER 信号肽)回补,实现**基因组编码的抗原递呈**
- **筛选流程**: 3NNK 肽库(MAGE-A3、gp100)基因组整合 → 可溶性 TCR-Fc(sTCR-Fc)结合 → FACS 分选结合/非结合群体 → 深测序富集分析,**3 轮选择**
- **ML 部分**: 深测序数据训练 MLP + "RH" 打分方法;筛了 2 个亲和力增强的 sTCR(含已获批药物 **Kimmtrak**);在全人类蛋白组尺度预测 off-target,实验验证了**与靶表位编辑距离大**的新 off-target 肽
- **注意**: 计划书把 ACDC 展开为 "Artificial Cell Display of pMHC",但 bioRxiv 原文缩写展开为 **"antigen-presenting cell detecting cytokine"**。🔴 **建议与实验室核对**——两者含义不同(前者强调递呈、后者强调细胞因子检测),答辩时被问到会有风险。

## 2. 所在实验室相关工作(平台/ML 语境)

- **TCR-Engine**(Immunity, 2022): CRISPR 工程 T 细胞系上高通量改造 TCR 功能与特异性;30 个 TCR、437 个单点变异、~260,000 组合变异;发现 TCR 结合与激活的**不一致性**(binding vs activation 分离)——这直接支持我们"激活信号为标签"的设计
- **Mason et al., Cell Systems 2024(视角文章)**: "用 ML 预测适应性免疫受体特异性本质上是**数据生成问题**";全库仅 ~702 个 TCR-pMHC 共晶结构(vs 蛋白结构 >20 万)——结构数据的稀缺是核心瓶颈,证明 AF3 结构预测 + 干湿闭环的必要性
- **TouCAN / CALM**: 对比学习 TCR 聚类与特异性预测(ESM 编码),免疫特异性基础模型方向
- **Engimmune Therapeutics**(2021 年衍生公司,2022 年种子轮 CHF 15.5M): 可溶性 TCR 疗法,基因组编辑+功能筛选+深测序+ML 管线

## 3. 竞争/互补平台(高通量 TCR-pMHC 互作图谱)

| 平台 | 出处 | 原理 | 吞吐 | 局限 |
|---|---|---|---|---|
| **ACDC** | Hong 2025 bioRxiv | HEK293 合成 APC + IL-2 报告 + 基因组编码肽库 | 蛋白组尺度预测 + 定向验证 | 需逐个 TCR 建系;噪声仍存 |
| **SABR / SABR-II** | Joglekar lab, Nat Methods 2019/2024 | 肽-MHC 与 TCR 信号域融合的双功能受体,NFAТ-GFP 读出 | SABR-II: 4,075 表位/库,~20 min/复孔×3/TCR | 文库规模受限于 Jurkat 转导 |
| **T-Scan** | Kula 2019, Cell | 慢病毒递送抗原库 + Granzyme B 切割荧光报告 | 全基因组尺度 | 只能检测杀伤(非全部激活)信号 |
| **MCR** | Akama-Garren, Sci Immunol | 肽-MHC-II 胞外域 + TCR 信号域融合,报告细胞富集 | 转录组尺度肽库 | 以 MHC-II 为主 |
| **Trogocytosis** | Li/Joglekar, Nat Methods 2019 | APC 从相互作用 T 细胞提取膜蛋白 | 单对识别 | 亲和力依赖 |
| **酵母展示(MHC-II)** | Birnbaum lab, JBC 2023 / eLife | 肽与 MHC-II 解耦表达,单库筛任意等位基因 | 11,040 SARS-CoV-2 肽 × 3 等位基因 | 需可溶性 TCR;MHC-I 展示难 |
| **VelociRAPTR** | bioRxiv 2025 | 酵母展示 + pMHC 假型慢病毒 VLP | **47M TCR 变异 × 92 pMHC** 同时筛 | 结合信号,非激活信号 |

**对本项目的含义**:
- ACDC 是唯一"基因组编码抗原 + 细胞因子激活读出 + 可工程化等位基因"的组合,信号最接近真实 T 细胞功能,这正是计划书噪声优势论的基础
- 竞争平台大多只测**结合**(SABR 类除外),激活读出是 ACDC 差异化优势
- 计划书成本断言(蛮力 ACDC ~$8M、液滴 scRNA-seq $1.5–3M/10^8 组合)来自实验室内部估算,文献中无直接可比数据 [未验证]——综述中建议改述为"实验室内部测算"

## 4. 对模型设计的启示

1. ACDC 标签是**激活富集信号**(FACS 分箱 + NGS 富集),不是亲和力常数——回归/排序目标优于二分类;信号去卷积(文库内串扰)是预处理重点
2. 数据分布:每轮库设计由 ML 决定 → **训练集分布随时间漂移**,必须用校准锚点(小规模重复验证)对抗 covariate shift
3. 3 轮 FACS 选择在结合筛选阶段引入选择偏差,富集值需做轮次间归一化
4. Kimmtrak 案例是现成的金标准测试案例(与我们 Module 3 一致)
