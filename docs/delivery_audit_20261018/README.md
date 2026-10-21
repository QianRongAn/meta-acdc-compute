# 交付审查复核与处置（2026-10-18）

结论：历史结构评分及其临床解读不能引用。报告指出的交付脱节是真实问题，尤其MAG-IC3/a3a混同、闸前临床结果、错误引用和硬编码图表；本轮已修正现行产物与消费入口。此前“28项已修”不能等同数据、模型、结果和交付都已闭环。

但报告混合了旧归档与当前分支，也含几条过强推断。不能把整份报告直接当已验证事实或执行指令。外部原文保留于`external_review.md`，输入哈希见`review_receipt.json`；39项的证据、当前状态和未完成事项见`adjudication.json`。

## 两张表究竟是什么

- 临床7行全部由本地找回的5张per-seed表取均值精确重建；`clinical_table_reconciliation.json`记录每个种子的原值和均值。例如靶肽0.89164→0.8916、titin0.94832→0.9483。
- 21行relabeled集成表的均值和SD，全部可由同包105行逐AF3模型分数精确重建（每job5个模型，SD为样本SD）；见`relabeled_table_reconstruction.json`。
- 因此两表不是“同一五种子输出的两份副本”。排序反转是真实的跨评分来源/汇总口径差异，来源与版本混用仍是交付缺陷。数值聚合找回不证明旧checkpoint当时如何产生分数，不能择一称为当前正确结论。

原目录找回了6张评分表、5个checkpoint及旧dataset。评分表以原始字节交付在`recovered/`，权重/cache与原数据在本地保留并单独打包。`recovered_inventory.json`记录哈希与体积，`legacy_checkpoint_metadata.json`显示5个seed存在但缺release指纹。没有为旧运行补造receipt、重训或重新给临床肽打分。

## 身份、引用与QC

官方RCSB5BRZ FASTA与实际CIF的五条链逐字相同，长度197/241/9/275/100。由[Raman2016原文](https://www.nature.com/articles/srep18851)明确：5BRZ/5BS0是MAG-IC3；临床a3a是另一工程构建体。同源不等于同一给药受体。未把MHC判成A02，也未否定原文0.285Å引用。

旧表按af3_predictions联合键有4/7临床行带QC flag且仍有分数。现行`kn11_clinical_scan.tsv`逐行WITHDRAWN，分数、靶分数、差值和hazard全部清空；旧内容原样归档。临床参考表与其生成器明确区分a3a、MAG-IC3、Morgan-A02、tebentafusp；不把病例结局扩展为模型标签。

新核查还纠正两条报告沿用的旧PMID：23470321实际为共刺激综述（Morgan目标论文为[23377668](https://pubmed.ncbi.nlm.nih.gov/23377668/)），34525232实际为儿科COVID研究（tebentafusp生存论文为[34551229](https://pubmed.ncbi.nlm.nih.gov/34551229/)）。Cameron为23926201，MAG-IC3结构论文为26758806。九肽KVAELVHFL不等于十肽KVAKELVHFL，后者继续未验证，不转挂到A02病例。

## 数值离散度能说明什么

中位SD0.2512、13/21>0.20、5/21>0.30以及5brz四个差值小于对应SD都已复现。但“差值<SD”不等于“没有任何信息”，也不是等效性或非显著性检验。需要配对差值、协方差、预定义效应/等效界限和明确实验单位。两种SD是不同汇总轴，不能直接相加或把它们称作校准后的测量误差。

本轮撤回的是没有证据支持的临床解读，没有从旧汇总表推一个新的“最小可辨差”来恢复结论。

## 版本差异与不应照做的建议

当前DA的8种复合物命名均归一到5brz，训练用release冻结的家族索引；AA3TO1覆盖诱饵比较、math导入、长度不等拒绝、缺骨架抛错、小界面拒绝等在本轮之前已存在。当前data包完整，torch2.9.1支持torch.func。

旧data.zip的.git/shallow明确记录a3da597边界，属于浅导出；当前仓库不是浅库，复核时237个可达提交，报告列出的6个提交全部可定位。交付浅克隆确实妨碍审计，但不能据此推断曾squash或提交不存在。新交付包含可达历史bundle。报告使用的原code.tar.gz未在可用路径找到，不能声称已逐字验证该归档内全部8处代码。

`NA`不是ASN糖基化的标准一字母占位表示；`GS`也不能自动当标签。[UniProt P04439](https://www.uniprot.org/uniprotkb/P04439/entry)的信号肽为1–24、成熟链25位即从GSHSM开头。不能为了迎合报告，把合法序列统一剪掉。构建体标签/变异必须有逐条原始依据。

188候选、72个PDB前缀和59个精确A/B链对已复现；全体角色不能仅凭链字母定为受体。347/73所依赖的确切pending面板与聚类规则未取得，不伪造相同计数。0.920/0.741/0.692的最大值不是0.933，后者与这三个种子的具体run关系也不能靠叙述补全。

## 本轮实际改变

1. 临床生成器、扫描入口、仪表盘及现行表不再产生或展示临床分数/危险标签；旧表原样保留并由`data/processed/artifact_status.json`登记撤回，改名复制也会因哈希被拒绝。
2. 退休三个旧DA shell入口，避免导入旧清单后再遇到缺release而失败；不启用监控、重训或自动派发。
3. 评分输出新增结构/QC/checkpoint/release哈希receipt；集成要求当前release、预声明种子、receipt及QC，明确between_seed_means和AF3_models_of_checkpoint_means，保留退化seed并使质量状态不通过。
4. 旧硬编码绘图脚本归档；当前入口仅从哈希绑定的开发比较报告作图，生成`derived_figures/`及provenance，不画旧0.933/0.644故事。图像已检查。旧图文件保留历史身份，不作为现行结果。
5. 增加必须存在的真实CIF/官方FASTA fixture、两表聚合重建、QC/撤回消费、来源receipt、临床构建体与回归测试；旧测试中“高分应产生fatal标签”的不支持规则也被替换为必须拒绝。

## 没有完成、也不能假装完成的事

- 缺旧运行时模型/数据/分数的完整绑定日志；恢复现有文件及重建均值不等于完整重演所有历史结果。
- 未取得报告原code.tar.gz及确切347面板，未重现其全部环境/统计。
- 没有新增真实阴性、独立临床结局验证、湿实验、已发表工具同任务比较或物理验证的诱饵；这些不是修代码能补出的证据。
- 保留原有严格新SPR验收门槛。新亲和力开发结果、已消费测试及未评分预留池原样不变；本轮没有训练、调参或再评分。旧结构release因入口代码变化需新版本指纹，继续拒绝，未重写旧manifest以伪装可用。

代码入口复核与结果可追溯性改进，不构成任何结构结合预测或临床模型的放行。

## 验证与交付

279项测试，276通过、3项既有可选fixture测试跳过；新增14项事故测试无跳过，另增加集成旧CLI拒绝测试。完整日志见`verification-tests-release.log`。首次回归发现的旧仪表盘schema错误及2项旧CLI契约测试失败保留在前序日志中，已修正后重跑通过。

原工作树清洁，6个关键历史输入、旧亲和力模型及其两份manifest、新旧冻结验收集均未变，见`integrity.json`。复用现有venv，在本独立工作树用`PYTHONPATH=src`运行，避免共享venv的旧editable路径指向原工作树；未重复安装。

代码和评分表随Git提交；旧权重/100MB缺陷cache作为本地追溯补充包交付，不伪装成新release，哈希与12项成员验证见`historical_package.json`。新代码快照与完整Git历史bundle另附本地输出，原314MB数据zip不覆盖。

## 第二轮修复（2026-10-21）：39 项中残余的工程缺口

第一轮（2026-10-18）完成撤回、身份、QC、追溯与测试主修复；以下为第二轮补上的
可执行工程项。没有恢复任何已撤回数字，也没有补造旧运行的证据。

| 复核项 | 第二轮处置 |
|---|---|
| #32 Fisher 非独立 | `rank_stats.py` 新增 `--families`：先按受体家族池化、再跨家族合并；原“按独立结构”口径保留但明确为上界。纯函数 + 单测锁定（`tests/test_review_round2.py`）。 |
| #14 固定 +z 位移 | `make_displaced_peptide_decoy` 改为沿“肽质心→MHC 沟槽质心”方向平移；缺 MHC 或质心退化即报错；docstring 仍明确这只是构造规则。方向依赖性与报错路径有单测。 |
| #17 链分类无测试 | 角色标注路径（`classify_chains` 只信 depositor/reference 注释）补边界与冲突单测；长度窗口启发式已从当前代码移除，不做序列裁剪。 |
| #20 无随机肽基线 | 新增 `meta_acdc.analysis.random_peptide_baseline`：构造匹配的随机 9 肽嫁接对照打分。仅在“已验收 release + 指纹匹配 checkpoint”下可运行；无有效模型时硬拒绝，绝不产出假数字。 |
| #2 delta 与 SD | 新增 `meta_acdc.analysis.resolution_floor`：种子轴与 AF3 模型轴分开报告 per-job SD 与对参考的 delta，只给描述性统计并显式声明“delta<SD 不是显著性/等效性检验”。 |
| #16/#33 面板口径 | 新增 `scripts/build_panel_registry.py` + `data/processed/candidate_panel_registry.tsv(+summary.json)`：已打分面板复现独立复核计数 188/72/59/28/12（逐项一致）。待传面板随后由 [`docs/forensic_followup_20261021`](forensic_followup_20261021/README.md) 定位并逐字核验为仓库 `af3_pending_upload/`（347 FASTA 与首加入提交 `4f730d2` 字节一致）：347/111/91、精确 αβ 对 98、MHC 45、B2M 11、NA 前缀 53、GS 前缀 279——与外部报告一致。本注册表待传面板改用 FASTA 角色头（alpha/beta/peptide/mhc/b2m），现输出 98/45/11 与核验一致。 |
| #1/#21 五种子来源 | [`docs/forensic_followup_20261021`](forensic_followup_20261021/README.md) 在隔离环境用原始源码 `5c70164` + 5 个旧 checkpoint + 884 份 CIF 重放 940 个 job×seed 均值：938 个与旧表四位小数一致（2 个末位差 0.0001）；另找回 5 份逐结构 per-seed 原表（4420 个结构×seed，4417 个一致）。分数保持 FORENSIC_WITHDRAWN；21 行表的原 checkpoint 仍未闭环（0/3 探测匹配），两套结果继续撤回。 |
| #34 清单反列+重复 | `load_manifest` 单点加护栏：重复行直接拒绝；旧 `kn5_submission_list.tsv`（tcr_a/tcr_b 反列）消费时强制 stderr 声明。`da_rescore_gpu_20261012.sh` 退役（会覆盖已找回的历史种子表）。 |
| #39 bib 占位 | ZeroBind 条目逐字核实（DOI `10.1038/s41467-023-43597-1`，PMID 38030641，标题/作者/卷页修正，注明其为 DTI 方法）；ADR 2026 综述 DOI 核实（`10.1016/j.addr.2025.115762`，PMID 41421504）；`xxxxx` 占位清除。 |
| #22 提交号审计 | AUDIT 引用的 6 个提交号在当前仓库全部可定位：`98773bd`（28 项主修复）、`f59dd45`（P1 修复）、`a3b8715`（跟进对齐）、`21926a9`（QC 表重建）、`60fb071`（二线修复）、`58df6ab`（复核回执）。 |
| #21/#23 | 5 张 per-seed 表 + 21 行集成表继续以原字节交付于 `recovered/`；第二轮新增回归测试 21 项（全部不 skip），覆盖 Fisher 分层、位移方向、清单护栏、面板计数复现、resolution_floor 与基线拒绝路径。 |
| #19 同任务外部工具 | 由 [`docs/six_blockers_20261021`](six_blockers_20261021/README.md)（`15eb5c6`）落地：官方 ERGO-I（固定 McPAS LSTM checkpoint，源码/权重 Git blob 与固定提交一致）在同一已消费 MIRA 功能观测面板运行 10,000 行并完整保留，逐肽宏平均 AUROC 0.5112（原 SPR 基线 0.5001）、宏平均 AP 0.1166（0.1197）。明确为共同功能终点的关联对照：不是同一 SPR 回归任务、不做显著性主张、不恢复“结构优于序列”结论；TCRMatch 因原始 VDJdb 快照无该面板人源 TRB 行而判不适用并如实记录。 |
| #21 21 行权重 | `15eb5c6` 预先限定 8 个不同哈希的兼容旧 checkpoint 对全部 105 份结构逐一重放：节点/边计数 105/105 一致，但分数仅 0 或 1/105 四位小数匹配——无完整匹配，不从偶然吻合认领原模型；旧表继续撤回。原审查 tar 仍未找到，用评分提交 `6797863` 重建源码包并记录哈希，明确重建包≠原档。 |

第二轮验证：`python -m unittest discover -s tests` 共 307 项全部通过（9 项既有可选 fixture 测试按
仓库惯例在本地原始数据缺失时跳过；第二轮 21 项 + 取证 7 项无跳过）。
另修复 3 项与第二轮无关的既有测试/入口缺陷：`test_cif_parse_and_identity` 改为当前身份分配路径
（`annotate_prediction`，审查 #17 后不再按链长/字母猜角色）；`attribute.py` 旧 CLI 拒绝无 release
指纹的旧 checkpoint 并要求显式 `--reference-json`/`--peptide`（测试改为断言拒绝）；
`test_fresh_spr_reserve` 在本地原始快照缺失时按既有模式跳过而非报错。
