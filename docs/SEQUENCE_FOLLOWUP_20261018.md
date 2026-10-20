# 完整工程DNA、D1/D6与SPR条件续查（2026-10-18）

**本轮补到1D3HMCys的公开完整工程DNA参考；没有补成临床实际批次证明。严格SPR池仍31/51条，原分区16/12/3，测试仅2个独立块，继续阻断训练。** 原10条测量条件、2条精确表值、2条温度未取得足够新证据，未因检索次数增加而放行。独立验收核心仍7条/5个受体谱系/4篇论文/3个来源组，全部未评分。新证据与失败响应分别登记，旧冻结材料和历史结果保留。

## 补到的完整工程DNA

追回[Lisa Borkner 2009博士论文](https://epub.uni-regensburg.de/13413/1/Doktorarbeit_v7_Bibliotheksexpemplar.pdf)，PDF130–132页/印刷123–125页附录6.2.1。该实验研究直接公开 **pMP71-1D3hmcys** 碱基序列，非反向翻译或按常数区模板拼接。

- 公开参考质粒按原文连续碱基提取得 **6338 nt**；长度是本次计算值。
- 原文β1103–2023、P2A2024–2089、α2090–2899，对应 **1797 nt编码插入（包含末端终止密码子）→598 aa蛋白**。
- β307 aa + P2A22 aa + α269 aa；α段另有末端终止密码子。β–P2A–α顺序、CDR3α `CALGGGADGLTF`、CDR3β `CASSLSGTTADEQYF`、α180C和β191C均核对一致。
- 逐行保留PDF页码、提取行、原始DNA及翻译；已看原PDF图面核对有色基因区、边界与工程位点。

[参考DNA与来源](sequence_followup_20261018/1D3_published_reference_DNA.json)、[DNA FASTA](sequence_followup_20261018/1D3_reference_DNA.fasta)、[蛋白FASTA](sequence_followup_20261018/1D3_reference_protein.fasta)。这是2009年研究质粒参考；2014[制造论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4208561/)及2015[病例原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4817886/)描述同名受体/MP71/鼠常数区/半胱氨酸/P2A，但不能据此认证临床批次逐碱基一致。病例原文明确报告MCB及最终输注产品做过序列完整性检查；尚缺的是公开批次序列与此参考DNA的比对证据。未把“作者做过QC”误写成“没有做过QC”，也未新增临床风险标签。

## D1/D6仍未补到完整可变区的原因

取得[原研究TableS1工作簿](https://ars.els-cdn.com/content/image/1-s2.0-S2211124724008842-mmc2.xlsx)：只有一个214×22的表，提供V/J、CDR3及供者计数，没有完整可变区或工程DNA。第36行对应D6；第38行对应D1配对CDR3；第37行同时出现两个β链及“or”注释，保留歧义，不擅自挑一条。详见[工作簿核查](sequence_followup_20261018/TW10_workbook_audit.json)。

[研究全文](https://ndownloader.figshare.com/files/48843463)PDF16页数据声明只列8F5A（KS1复合物）与8F7M（pMHC突变体），额外资料需向作者索取；本轮公开来源未取到D1/D6全变量链。另一个PDB **8V4Z也叫D1**，却识别流感肽LPFEKSTIM，V基因为TRAV1-1/TRBV7-8，与TW10受体不同，不能替代。D1/D6保持配对CDR3/VJ及真实25°C SPR证据的保留候选，未补造胚系序列，未加入完整身份验收核心。

## 不同来源SPR的实际处理

新取回[2026 Science TCR156研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC13004167/)全文XML及32页原始补充材料。原文为人PBMC来源的TCR156、TLMSAMTNL/HLA-A*02:01，Garcia/Witte等来源组；结构库V基因注释为TRAV12-2/TRBV9，但该注释尚未独立复核到V等位基因。此组合**已有于保留核心PIK3CA TCR4**，因此即使区别于上一轮TRBV19受体，也不能算新增V基因家族；只能算新的靶肽/研究来源候选。

|受体变体|原始平衡Kd，µM|相邻动力学计算Kd，µM|配对序列来源|
|---|---:|---:|---|
|WT|30|24|9NMU|
|S32Hα|61|98|9NMV|
|S32Mα|5.3|2.7|9NMX|
|S32Qα|9.8|6.1|9NMY|
|S30E32Qα|11.5|7.8|9NMW|

补充PDF9页明确25°C、BiacoreT100、链霉亲和素芯片约200RU、TCR分析物经3C去拉链；PDF18页FigS6A表值已按图面核查。**原文固定相写的是pMHC tetramer**，未默改为单体，也没有补填未报告的缓冲液、流速和重复数。平衡Kd和koff/kon计算值分开保存；S6B的S30Q标注与主文/表/PDB的S32Q不一致，单列为原文标注问题。发表位点32对应该沉积可变区序列第33位，未盲目按硬编码位置造突变。

这5条作为[待澄清测量格式的SPR保留集](sequence_followup_20261018/PAP_SPR_pending.json)，**不是5个独立受体，也未增加严格核心计数**。原91条亲和力和结构v1–v7没有命中；原保留的IMMREP25元数据有3条保守同链/域家族匹配，逐条公开，不以新PDB号声称从未暴露。WT及4变体均加入强制保留库，禁止训练或模型选择；未运行模型评分。S32E/S30E另有原表数值，但未独立取得对应完整沉积变量链，单列数值证据，未按WT替换位点补链。

其他候选：[KRAS-G12V/HLA-A3论文](https://doi.org/10.1002/eji.202451079)明确是HLA转基因小鼠产生的**鼠TCR**，人MHC不代表人受体，排除人源核心。[2025交联研究](https://doi.org/10.1016/j.jbc.2025.108445)原始TableS1有5个25°C SPR nM数值，但均为亲和力成熟受体，尚缺完整测量方法/谱系独立性，不因字母A–E不同就计为新家族。HD14原方法仍未取到；数个新PDB的主引用为ToBePublished，不能据结构当作实测SPR。见[全部候选判定](sequence_followup_20261018/external_candidate_dispositions.json)。

## 仍缺的10条条件、2条表值、2条温度

四条旧SI已能证实精确数值（C12C两条、CF34/LC13两条），缺的是当前实验完整条件，不能称“什么原文也没有”，也不能称“已解决”。

|记录|受体/来源|仍缺与本轮处理|
|---|---|---|
|ATLAS:1|LC13；PMID15640805|原论文精确值与当前SPR条件未完整核对；Nature返回HTML，JCU公开索引附件此前401|
|ATLAS:1058|5934；PMID42190883|Rac1研究公开索引/订阅入口不足以核对原始实验条件|
|ATLAS:508|DD31；PMID25355921|公开仓储与出版社尚未取到足够的原始SPR方法与表值|
|ATLAS:514/515|C12C；PMID23521884|SI精确4.58/1.10µM已核实；当前SPR温度未明确；Cell主文403，不能借用热变性25–95°C或旧协议温度|
|ATLAS:532/537|CF34/LC13；PMID19167249|SI精确8.89/15.07µM已核实；SI只有表，无完整SPR方法；Cell主文403|
|ATLAS:542|SB27；PMID15849183|当前SPR原文条件未取得；JBC403，旧SI只有结构图|
|ATLAS:634|B7；PMID15670602|原始热力学/SPR方法未取得；检索命中的出版入口仍为访问错误，不能把响应内数字当表值|
|ATLAS:91|RA14；PMID19542454|AAI/OUP入口403；HAL索引没有可获取主文件|
|ATLAS:149/154|LC13/SB27；PMID20483993|同一个TableS2的10.4/11.2µM尚未核到；已有25°C方法不替代表值；PNAS403，EuropePMC附件请求是296字节错误XML|
|ATLAS:714/716|LS01/LS10；PMID28250417|已查主文SPR段及SI没有明确当前SPR温度；不能用结晶室温、仪器默认值或相关论文推断|

[逐条缺口](sequence_followup_20261018/remaining_SPR_gaps.json)与[51条原始账本](sequence_followup_20261018/review_ledger.json)保留相同fold/block。上述14项需要对应原文Methods/附件、原始数值表，或可审计且明确指向本次实验的作者说明；访问失败不等于资料永久不存在。本轮没有发送索取邮件或申请受限访问。

## 临床完整构建体尚缺什么

|对象|已经有|确实尚未补到|
|---|---|---|
|1D3HMCys|本轮6338nt研究质粒/1797nt完整参考插入|与临床MCB、输注批次逐碱基对应的公开序列证明|
|WT1 C4|既有1809nt/602aa完整公开参考插入|完整临床载体及批次对应；不能换成C4-DLT|
|HPV16 E7|既有1800nt参考插入，临床协议写7310bp MSGV1/PG13|7310bp完整载体DNA及实际临床构建体逐碱基对应|
|KRAS10mer/9mer|10mer599aa完整工程蛋白参考、组合治疗原文|完整临床工程DNA；9mer A/B/C中实际使用哪一个；不能按蛋白反向翻译或以相同V基因选构建体|
|p53 4196-AV6/BV11|完整配对变量片段及4196来源→4349治疗患者对应|最终鼠常数区/P2A/MSGV1完整临床编码DNA|
|Foy37|36对一致变量片段，TCR919冲突隔离|完整工程DNA与批次对应；论文“full length”片段不补造成整链|
|DMF4原始试验|主文α-IRES-β/MSGV1AIB|完整实际DNA；不能替换成后来P2A格式|
|afami-cel|临床主文和同一工程谱系|完整实际工程DNA及临床序列证明|

[临床缺口账本](sequence_followup_20261018/clinical_remaining_gaps.json)记录当前与历史检索来源。索引病例、同一受体多名患者、重复发表的同一病例不增加独立受体数量；参考DNA新增不改变既有临床结局证据等级。原临床18例冻结集原样保留。

## 验证与发布

新增14项边界检查，完整223项：**220通过、3项因旧结构夹具缺失跳过，无失败**。覆盖DNA坐标/读框/终止/P2A/工程Cys/两条CDR3、禁止冒充临床批次、编号TCR FASTA的明确链选择、禁止用相邻动力学Kd替换平衡Kd、pMHC四聚体格式和缺失条件保留、强制保留库不可通过空自定义库绕过。

本地检查发布v8继续`training_allowed=false`：31条、原分区16/12/3，测试仍2个独立块。未安装软件、未重排分区、未拟合优化器、未评分验收集。全部旧冻结文件、原工作树与6项关键历史输入做SHA/状态校验；原始资料与失败响应本地保留，GitHub提交派生证据、来源URL/SHA、代码和测试记录。见[完整性证明](sequence_followup_20261018/integrity.json)、[测试日志](sequence_followup_20261018/unit-tests.log)、[封存清单](sequence_followup_20261018/review_manifest.json)。
