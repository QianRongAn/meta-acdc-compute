# SPR续查与独立受体验收资料（2026-10-18，第二轮）

本轮继续核查剩余18条SPR原文、4条温度和临床工程序列。**补全原51条中的2条，另2条追回正文但精确数值仍缺；新增保留3条真实SPR测量，来自2个受体谱系。没有重训、拟合或评分外部验收集。没有重新洗牌。**

当前严格池 **25/51条**，原训练/验证/测试分割为 **16/8/1**，连通块为10/5/1。测试块仍不足3，训练闸门继续关闭；新增外部面板也尚不足以支持稳健的泛化结论。版本v5只供检查，旧发布和原始数据均保留。

## 原51条的新核查结果

|记录|原文核查与处置|
|---|---|
|ATLAS:724 F50 / GILGFVFTL / A*02:01|[PMID28931605](https://pmc.ncbi.nlm.nih.gov/articles/PMC5683187/) Results、Figure1B、SPR Methods明确76±4µM、BiacoreT100、25°C，生物素化pMHC固定SA芯片1000RU、TCR分析物；1:1平衡拟合。进入原训练fold；误差按原文保存，不转CI|
|ATLAS:736 TCR589 / IPLTEEAEL / B*35:01|[PMID30053426](https://pmc.ncbi.nlm.nih.gov/articles/PMC6140336/) Figure2B和Results明确4µM；STAR Methods明确BiacoreT100、25°C、固定pMHC100–300RU、SA芯片、HBS-P+加0.1%BSA、平衡拟合。回到原测试fold。原文还显示有结合却无刺激，不能赋予功能阳性|
|ATLAS:149 LC13、ATLAS:154 SB27|[PMID20483993](https://pmc.ncbi.nlm.nih.gov/articles/PMC2890827/) 正文及SPR Methods已取得，明确25°C/Biacore3000。精确10.4/11.2µM数值所在TableS2仍未取得；Figure1曲线已人工查看，不能据曲线臆造精确数值。以`quarantine_primary_value_unverified`继续隔离|

[51条账本](spr_followup_20261018/review_ledger.json)保存每行原始字段、继承fold、原文定位和SHA；[可读清单](spr_followup_20261018/review_51.tsv)、[26条未放行记录](spr_followup_20261018/unresolved_records.json)分开保存。统计：25条通过、14条仍缺足够原文、2条正文已得但数值表缺、4条温度未证实、2条肽数值错配、1条完整TCR-CD3格式分存、1条最初测量来源待追溯、2条方法/引用冲突。

## 新增独立验收候选与实际接纳

|谱系/来源|真实测量|配对序列、条件与判定|
|---|---|---|
|T4H2，[2021PNAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC7848747/)，PMID33468649|IMDQVPFSV/HLA-A*02:01，1.8µM|25°C、BiacoreT200，TCR固定CM5、pMHC分析物，1:1全局平衡拟合。具名TCR双链来自PDB6VM9；保留外部，禁止训练|
|SILv44，同一2021研究|IMDQVPFSV/HLA-A*02:01，210µM|同上实验条件；具名配对TCR序列来自6VM7。6VM7晶体肽是WT；这里只借其配对受体身份，未把WT晶体当成Met肽结构|
|T4H2，[2024PNAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC11295079/)，PMID39042676|IMDQVPFSV/HLA-A*02:01，6.8±0.5µM|25°C、BiacoreT200，pMHC固定CM5、TCR分析物；单结合位点平衡拟合，误差是拟合SE。与2021同受体同肽，不是第二个受体，不跨研究取平均|

**新增为3条观察、2个按现有保守序列规则划分的受体谱系、2篇论文、1个实验室来源组。T4H2与SILv44都使用TRBV19，不能称为新增两个不同TRBV基因家族。** 它们与旧91条KD、MIRA/临床保留库、前轮新增保留库及历史结构v1–v7未命中现有对应链/CDR3近缘排除规则；也与之前外部核心互不命中。历史v1–v4已退役结构仍参与检查，避免只与幸存结构比较。[独立性证据](spr_followup_20261018/external_independence_audit.json)。这不证明未知预训练语料中绝无暴露。

外部常规核心从原4条/3谱系/2论文扩为 **7条/5谱系/4论文，按来源合并仅3个研究组**。原4条逐字段保留；新增3条未评分、未训练、未选择模型。之前6条HLA-E/ImmTAV挑战仍单独保存，不混入常规核心。T4H2的1.8与6.8µM差异必须作为条件敏感性报告；不能把测量差异藏在均值里。未来评估须按受体及研究来源聚类，当前样本量不能支持精确泛化CI或临床风险解释。[冻结协议](spr_followup_20261018/external_test_protocol.json)、[7条定量核心](spr_followup_20261018/external_classical_panel.json)、[仅本轮新增](spr_followup_20261018/external_classical_additions.json)。

### 找到了但未充入独立测试的数据

- [P1-15，JCIInsight202235](https://insight.jci.org/articles/view/202235)：有真实25°CSPR；PDB28IL确实存在。其beta近缘命中历史训练NYN-I/8YE4及AF3版本，不能当新独立家族。原文给出的betaCDR3 motif与沉积序列还存在差异（CASSSGGGYEQYF对CASSESGGYEQYF），未自动替换。
- [NYN-I/NYN-II，2024JBC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11342781/)：NYN-I对应8YE4早已用过；NYN-II也未获得通过近缘排除的证据，未凭换克隆名宣称独立。正文13.5µM不被二次引用13.6覆盖。
- [TP53R248Q，JCI196613](https://www.jci.org/articles/view/196613)：真实SPR1.86µM；SPR段未报温度，其他段25°C是ELISA。已找到关联专利CN118772260B，但被测克隆、精确配对序列及完整表达构建的对应仍未核实，继续隔离；与临床p53R175H受体不混同。
- 2024T4H2文中8±2µM对应**非天然norleucine肽**；不编码成L或M，不当作IMDQVPFSV重复测量。整条spike蛋白不与TCR结合也不能生成肽–MHC不结合标签。

完整[候选判定清单](spr_followup_20261018/external_candidate_dispositions.json)。

## 临床构建核查结果与真实缺口

补取[1D3HMCys制造原文PMID25143008](https://pmc.ncbi.nlm.nih.gov/articles/PMC4208561/)和[致死事件原文PMID25896248](https://pmc.ncbi.nlm.nih.gov/articles/PMC4817886/)。确认pMP71、鼠恒定区、额外半胱氨酸、密码子优化、P2A等设计；论文明确报告实施过MasterCellBank/输注产品的序列完整性校验。因此缺口的准确表述是“**本次未取得公开完整编码序列和批次特异的测序报告**”，不能写成研究者没有做过序列验证。2015病例仍是2022队列同一患者，未重复计数。

[KRAS临床全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9531755/)确认9mer GADGVGKSA与10mer GADGVGKSAL分别转导、合并输注，一名患者的组合疗效不能拆成两个独立受体成功案例。PMC索引有补充附件，但实际返回HTML检查页；NEJM附件403。未取得可核对的完整临床插入片段，未据构型描述拼接参考DNA。

此前C4完整公开参考插入1809nt/602aa、E7前临床插入1800nt和Foy36对一致可变区片段继续保留；本轮**未新增完整临床工程编码盒**。仍缺Foy完整表达盒且TCR919冲突、p53-4196最终鼠恒定区/P2A/MSGV1构建、1D3完整编码插入、KRAS9mer/10mer精确临床链配对和插入片段。C4/E7参考盒也不是临床批次测序证明。逐项来源、已具备资料和缺项见[临床缺口](spr_followup_20261018/clinical_remaining_gaps.json)。旧18候选临床验收冻结完全保留，未评分。

## 4条实验温度仍不能补填

LS01/LS10（PMID28250417）及NYBR1/a7b2（PMID32310221）重新读取全文SPR Methods；LS正文有结晶室温，NYBR1引用旧实验协议，但没有足够依据把当前SPR温度确认为25°C。先前补充材料也未核实到温度。不借用MD298K、结晶温度或仪器默认值。四条按当前严格温度政策继续隔离；这不是宣称原测量无效。[温度核查](spr_followup_20261018/temperature_followup.json)。

## 剩余14条原文条件未满足

以下为记录数；重复论文仍分别列对应观察。公开索引给出付费/订阅、出版社403或资料库401/403时保留访问记录；不能把HTTP200的摘要页或验证页当全文，也不能把这些访问失败说成文献永久不可得。

|记录|受体|PMID|当前可补性|
|---|---|---|---|
|ATLAS:1|LC13|[15640805](https://pubmed.ncbi.nlm.nih.gov/15640805/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:1058|5934|[42190883](https://pubmed.ncbi.nlm.nih.gov/42190883/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:145|LC13|[20064448](https://pubmed.ncbi.nlm.nih.gov/20064448/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:148|LC13|[20064448](https://pubmed.ncbi.nlm.nih.gov/20064448/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:508|DD31|[25355921](https://pubmed.ncbi.nlm.nih.gov/25355921/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:514|C12C|[23521884](https://pubmed.ncbi.nlm.nih.gov/23521884/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:515|C12C|[23521884](https://pubmed.ncbi.nlm.nih.gov/23521884/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:532|CF34|[19167249](https://pubmed.ncbi.nlm.nih.gov/19167249/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:537|LC13|[19167249](https://pubmed.ncbi.nlm.nih.gov/19167249/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:542|SB27|[15849183](https://pubmed.ncbi.nlm.nih.gov/15849183/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:634|B7|[15670602](https://pubmed.ncbi.nlm.nih.gov/15670602/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|
|ATLAS:727|D30|[28945243](https://pubmed.ncbi.nlm.nih.gov/28945243/)|已得部分补充表；主文Methods/温度未取到，资料库返回403|
|ATLAS:728|D30|[28945243](https://pubmed.ncbi.nlm.nih.gov/28945243/)|已得部分补充表；主文Methods/温度未取到，资料库返回403|
|ATLAS:91|RA14|[19542454](https://pubmed.ncbi.nlm.nih.gov/19542454/)|出版社/已查公开资料库未返回足够正文与SPR条件；元数据或预览不能放行|

## 校验与保存

新增具名TCR FASTA选择器，必须同时匹配克隆名、受体标记与alpha/beta角色，拒绝误取MHC alpha链。新外部资料校验KD单位、肽字母、配对CDR3、温度、文件SHA和沉积链身份；缺项及任何训练/风险/二分类标签提升都会拒绝。新增受体进入强制保留库，即使传空自定义库也不能绕过训练排除。

本轮全套193项检查：190通过、3项因旧结构测试夹具缺失跳过；无失败。原始源文件保存在本地版本化快照；GitHub提交证据账本、来源URL/SHA、判定与代码，不把错误PDF后缀或原始数据库重新包装为有效数据。旧61份冻结文件、原工作树及6项关键历史输入均另做完整性校验。详见[测试日志](spr_followup_20261018/unit-tests.log)、[完整性结果](spr_followup_20261018/integrity.json)和[审核封存](spr_followup_20261018/review_manifest.json)。
