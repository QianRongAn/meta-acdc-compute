# SPR原文、温度、不同V基因候选与KRAS工程参考（2026-10-18）

**本轮新增6条达到严格池条件：原14条待补原文中补全4条，原4条温度中解决2条。当前31/51条，原分区16训练、12验证、3测试；测试仅2个独立块，训练继续阻断。** 没有重排fold、拟合模型、评分外部验收集或改动历史实验结果。最终本地检查发布v7；旧发布保留。

## 新增达标证据

|记录|原文与条件|判定|
|---|---|---|
|ATLAS:727/728，D30，GTSGSPIVNR/GTSGSPIINR|[帝国理工接受稿](https://spiral.imperial.ac.uk/entities/publication/c55987f2-6336-410a-9d2c-e37c11088237)，PDF30页Table1：136.0±7.0/142.5±2.5µM；PDF20页SPR：25°C、Biacore3000/T200，抗体捕获TCR于CM5，pMHC分析物1.56–400µM，1:1模型|保持原验证分区；误差是3次独立实验各重复测量的SEM，不改成SD或独立受体数|
|ATLAS:145/148，LC13，EEYLKAWTF/EEYLQAFTY|[原研究补充材料](https://ars.els-cdn.com/content/image/1-s2.0-S107476130900510X-mmc1.pdf)，PDF17页TableS3：1.54±0.04/49±0.20µM；PDF4页Methods：25°C、Biacore3000，W6/32捕获HLA-B*44:05于CM5，LC13分析物。FigureS3对应mimotope/allotope肽身份|保持原验证分区。保留直接Kd列；不换成相邻的动力学计算Kd 1.90/21µM；±类型未核实，不转置信区间。所需原始资料在补充材料内齐备，不冒称已取得整篇主文|
|ATLAS:747，NYBR1，47pM；ATLAS:751，WT1-a7b2，70nM|[共同作者RoryCrean博士论文](https://purehost.bath.ac.uk/ws/portalfiles/portal/202462532/Rory_Crean_Thesis.pdf)，同题Chapter4；PDF97/108页对应相同受体与数值，PDF115页明确SPR全在25°C、PBS、30µL/min、Biacore3000、链霉亲和素偶联CM5、single-cycle|补足温度，保留原期刊数值来源和原测试分区；明确温度来自一手作者论文，不是假称期刊Methods已写出温度。两条仍属同一研究/同一split block|

[51条证据账本](diversity_followup_20261018/review_ledger.json)、[可读TSV](diversity_followup_20261018/review_51.tsv)、[20条仍未放行清单](diversity_followup_20261018/unresolved_records.json)保留原始记录、原fold/block、新来源定位和SHA。

## 不同V基因、不同团队的真实SPR

追回Gras团队[2024 CellReports研究](https://doi.org/10.1016/j.celrep.2024.114555)的[大学仓储全文](https://ndownloader.figshare.com/files/48843463)及[补充材料](https://ars.els-cdn.com/content/image/1-s2.0-S2211124724008842-mmc1.pdf)。正文PDF18页明确25°C、BiacoreT100/T200、固定HLA-B*57:01–TW10于SA约2000RU，TCR作分析物；补充TableS2/S4提供配对CDR3/VJ与Kd。已核对表格和Methods图面，PDF文本抽取误把µM读成mM，派生数据按图面µM保存。

|受体|V基因|TW10真实SPR Kd|当前用途|
|---|---|---|---|
|D1|TRAV4/TRBV20-1|47.7±3.3µM|保留为独立验收候选；完整双链可变区未取得，尚不能进入现有模型的严格核心|
|D6|TRAV4/TRBV20-1|127.1±7.7µM|同上；同一研究、同一V基因组合的第二克隆，不是第二团队|
|KS1|TRAV21/TRBV5-1|226.9±8.1µM|命中旧91条KD及历史8F5A结构/AF3，禁止充作新独立验收受体|

表中误差是至少2次独立实验、各做重复测量的SEM。D1/D6与上轮TRBV19不同，来自Gras/LaTrobe–Monash/UWA来源组；但该团队也参与过开发集研究，不宣称它是完全未见研究团队的挑战。D1/D6配对CDR3在旧91条、MIRA/临床保留库、之前追加库和所有历史结构v1–v7下无现有近缘规则命中。这仅是已有序列规则下的排除结果，不能证明未知完整链绝无近缘关系。**未从V/J和CDR3拼接参考胚系补造完整可变区，也未降低严格核心的身份条件。**

两个新候选已进入强制保留库，防止之后混入训练。当前可用于既定模型的外部核心仍是**7条/5谱系/4论文/3来源组，未评分**；不能把D1/D6算进去宣称已扩充就绪核心。[候选与来源](diversity_followup_20261018/external_differentV_candidates.json)、[逐库重叠核查](diversity_followup_20261018/TW10_identity_audit.json)、[验收协议](diversity_followup_20261018/external_test_protocol.json)。

同文的T242N“NB”结果保留原有限定：不写无限大KD，不自动赋二分类阴性。WT实验TCR最高浓度500/430µM已核实，但未擅自转成突变体实验的检测限。

另找到[Ragon/Walker/Birnbaum团队IW9研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC12889509/)，原文有25°CSPR和1.53/1.92/6.53/7.98µM四项数值。该版本是**未同行评审预印本**；具体被测配对序列/V基因仍未取得，附件返回HTML验证页、bioRxiv返回429。未凭四个TCR名字凑四个独立家族。G9是MHC-II反向、肽无关结合且SPR温度未报，另存；2025PIK3CA新论文仍对应既有TCR3/TCR4，不重复计新谱系。[候选隔离原因](diversity_followup_20261018/external_candidate_dispositions.json)。

## KRAS完整工程蛋白参考，仍缺精确临床DNA

从[EP3494133B1](https://patentimages.storage.googleapis.com/a6/63/45/76d90e3a55d9a0/EP3494133B1.pdf)SEQ56/57/58核出KRAS G12D 10mer受体的**alpha269aa、beta303aa、连接肽27aa**。全部三字母残基行逐一对照PDF60/61/63/64页；现有WindowsOCR仅辅助定位，未直接信任其混乱的列顺序及识别字符。专利同族[Example2/Table7](https://patents.google.com/patent/US20240209059A1/en)明确TRAV12-2/TRBV10-2配对及alpha–连接肽–beta顺序；据此保存**599aa公开工程蛋白参考**，含鼠恒定区及原文工程改造，而非只有CDR3。

[参考序列与逐行转录](diversity_followup_20261018/KRAS10mer_published_protein_reference.json)、[FASTA](diversity_followup_20261018/KRAS10mer_reference.fasta)明确区分：公开参考蛋白、完整编码DNA、完整载体、2022实际输注批次。这次补到前者，**未取得后面三者**，未反向翻译制造所谓原始密码子优化DNA。

9mer还有TRAV4/TRBV5-6的A/B/C三个公开构建，不能只凭同一靶点或V基因挑一个对应临床9mer。2022病例是一名患者的9mer/10mer组合治疗，不把这份专利参考自动绑定临床疗效或拆成两个独立病例。

其余临床缺口仍是：Foy公开表是可变区片段，未包含完整工程盒且TCR919冲突；p53-4196缺最终鼠恒定区/P2A/MSGV1精确插入；1D3制造论文有设计和序列校验说明，但公开完整DNA/批次序列报告未取得。C4完整公开编码参考和E7前临床编码插入仍保存；均不冒充临床批次序列证明。[临床逐项缺口](diversity_followup_20261018/clinical_remaining_gaps.json)。旧临床18候选验收冻结和结局归属未改动、未评分。

## 尚未补齐的项目及原因

原4条温度缺口只剩**LS01/LS10**。主文SPR段及已查补充材料均没有当前SPR的明确温度；结晶室温、培养温度和仪器惯例不能补填。新补足的NYBR1/a7b2有同题一手作者论文的直接25°C声明，证据等级与猜默认温度不同。[温度账本](diversity_followup_20261018/temperature_followup.json)。

**ATLAS:149/154**仍缺PMID20483993同一个TableS2中的两条精确数值；全文已有25°C方法，Figure1曲线不能替代10.4/11.2µM精确表值。出版社附件403，PMC附件返回HTML验证页，作者仓储附件401；不把这些响应当PDF。

剩余**10条原文条件不足**中，本轮又追回4条数值：C12C两条由补充FigureS3核实4.58±0.11/1.10±0.05µM；CF34/LC13两条由TableS1核实8.89±1.31/15.07±1.67µM。但C12C补充Methods没有当前SPR温度，CF34材料只有数值表；仍不放行。C12C的25–95°C是热稳定性实验，未冒充SPR条件。

|记录|受体|来源|缺项或访问结果|
|---|---|---|---|
|ATLAS:1|LC13|[PMID15640805](https://pubmed.ncbi.nlm.nih.gov/15640805/)|JCU公开索引有PDF链接但实际401；出版社预览/403；不能当全文完成|
|ATLAS:1058|5934|[PMID42190883](https://pubmed.ncbi.nlm.nih.gov/42190883/)|出版社/公开索引只有订阅入口，未找到可核对实验条件的开放原文|
|ATLAS:508|DD31|[PMID25355921](https://pubmed.ncbi.nlm.nih.gov/25355921/)|JCU公开索引有PDF链接但实际401；出版社预览/403；不能当全文完成|
|ATLAS:514|C12C|[PMID23521884](https://pubmed.ncbi.nlm.nih.gov/23521884/)|FigureS3精确数值已核实；Methods只引用旧SPR协议，当前温度未明确；ORCA主文403|
|ATLAS:515|C12C|[PMID23521884](https://pubmed.ncbi.nlm.nih.gov/23521884/)|FigureS3精确数值已核实；Methods只引用旧SPR协议，当前温度未明确；ORCA主文403|
|ATLAS:532|CF34|[PMID19167249](https://pubmed.ncbi.nlm.nih.gov/19167249/)|原始补充TableS1数值已核实；完整SPR方法和温度未取得|
|ATLAS:537|LC13|[PMID19167249](https://pubmed.ncbi.nlm.nih.gov/19167249/)|原始补充TableS1数值已核实；完整SPR方法和温度未取得|
|ATLAS:542|SB27|[PMID15849183](https://pubmed.ncbi.nlm.nih.gov/15849183/)|JCU公开索引有PDF链接但实际401；出版社预览/403；不能当全文完成|
|ATLAS:634|B7|[PMID15670602](https://pubmed.ncbi.nlm.nih.gov/15670602/)|出版社/公开索引只有订阅入口，未找到可核对实验条件的开放原文|
|ATLAS:91|RA14|[PMID19542454](https://pubmed.ncbi.nlm.nih.gov/19542454/)|期刊/OUP403；HAL记录只有文献信息，没有可获取的主文件|

公开入口访问失败并不表示文献永久不可得。新一轮[作者仓储路线核查](diversity_followup_20261018/remaining_primary_routes.json)、全部下载尝试与内容分类有记录；成功的SI与失败的HTML分别登记。

## 验证和发布边界

完整检查209项：206通过、3项因旧结构夹具缺失跳过，无失败。新增检查覆盖蛋白转录/链顺序、禁止反向翻译冒充DNA、候选不能绕过保留库、KS1历史暴露、NB不能变成绝对阴性、LC13直接Kd与计算Kd分列、同研究两测试记录不算两个独立块。31条池的原fold和block全部继承，测试3条仅2个block，故训练闸门仍关闭。

原工作树、6项关键历史输入与全部旧冻结文件校验不变；原始PDF/序列下载和历史结果在本地版本化保存。GitHub提交来源URL/SHA、完整派生证据、判定及校验代码。见[完整性证明](diversity_followup_20261018/integrity.json)、[测试日志](diversity_followup_20261018/unit-tests.log)、[封存清单](diversity_followup_20261018/review_manifest.json)。
