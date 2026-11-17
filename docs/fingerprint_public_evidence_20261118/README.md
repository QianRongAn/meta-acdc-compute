# KRAS 构建体原件、3995 读段与扩大结构检索（2026-11-18）

本轮定位了 3995 对应的公开样本和外显子运行，扩大到 1,125 项结构元数据。Sim 2020 原始 SI、能绑定 6ULN／6ULR 的订购插入序列或实物质粒测序、9d 原始配对读段仍待取得。新增合格独立跨研究正对照为 0，评估尚未冻结。

## Sim 2020 补件与晶体表达原件

[Sim 2020 正文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7293613/)列出了 pnas.1921964117.sapp.pdf，文件约 1.8 MB。Europe PMC 的文章记录为 hasSuppl=Y、isOpenAccess=N，其 supplementaryFiles 接口明确返回文章非开放获取的错误。PMC、旧 NCBI、Europe PMC 和 PNAS 附件入口本轮返回 HTML 挑战页、访问拒绝或下载准备页。实际浏览器打开 PMC 附件也停在 reCAPTCHA 页面。HTTP 200 的 HTML 没有按成功 PDF 计入。[文章记录](provenance/sim_europepmc_record.json) · [接口错误原文](sources/sim_europepmc_supp.zip) · [全部请求状态和哈希](retrieval_receipts.json)

OpenAlex 记录的来源仍为期刊、PubMed 和 PMC。DataCite 按本论文关联标识符检索的 5 项结果均为该研究的 PDB 记录。两项目录检索都没有给出表达质粒、订单文件或测序原件。GenScript 公开站点的本轮检索返回自动访问拒绝。现有论文方法只说明供应商、胞外段合成及 pET30a 的 NdeI/XhoI 克隆体系，缺少能绑定到具体晶体样品的插入片段和质粒记录。尚不能判断 SI 是否包含相关序列，也不能断言作者或供应商没有保留原件。

本轮从正确的 [wwPDB 版本化仓库](https://www.wwpdb.org/ftp/pdb-versioned-ftp-site)取得了 6ULN 和 6ULR 目录及 v1.4 坐标。两份修订记录从 2020-05-27 的 v1.0 开始，后续均为主版本 1，目录目前各仅列出 v1.4。wwPDB 规则将聚合物序列变更计入主版本变化，因此公开记录未显示一次后续序列大版本修订。原始 v1.0 文件本轮仍未取得，这项版本推断保留为间接证据。[6ULN 目录](history/6uln_versions.html) · [6ULR 目录](history/6ulr_versions.html) · [修订记录解析](audit.json)

α M10V 的差异位置仍对应两份结构中有原子的 D13，β D1M 对应只在 SEQRES 声明的 E1。公开克隆核酸分别支持 α ATG 和 β GAC。PDBe 当前实体序列与此前比对一致，TCR 的 mutation_details 仍为空。结合上一轮逐密码子核验，可把差异定位在公开克隆与沉积构建体之间。具体发生在设计、合成、克隆、蛋白末端处理还是沉积声明，现有原件不足以区分。β 起始 Met 与 NdeI 构建解释相容，未把这一解释标为已证实。[6ULN 实体记录](provenance/6uln_pdbe.json) · [6ULR 实体记录](provenance/6ulr_pdbe.json) · [上一轮完整序列核验](../fingerprint_kras_crystal_constructs_20261117/README.md)

## 3995 的原始数据追踪

[Tran 2015 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7445892/)给出 PRJNA298330。本轮读取全部运行并核对每个样本别名，确认下列映射。

| 来源 | 公开标识与核查结果 |
| --- | --- |
| 3995 样本 | SAMN04158402，ENA XML 的 alias 为 3995，NCBI 对应 SRS1112231 |
| 肿瘤外显子 | SRX8130908／SRR2637256，标题为 3995 exome tumor |
| 正常外显子 | SRX1331481／SRR2637257，标题为 3995 exome normal |
| 本项目唯一 RNA-Seq | SRX2148238／SRR4203093，原始记录的样本别名为 3948 |
| 项目全清单 | 17 次运行，16 次 WXS、1 次 RNA-Seq |
| 另外两个相关项目 | PRJNA298310 有 4 次运行，PRJNA298376 有 1,161 次运行，保留公开清单，尚未建立到 3995 同一 T 细胞克隆的对应关系 |

[3995 样本 XML](samples/SAMN04158402.xml) · [原始实验记录](samples/3995_sra_records.xml) · [RNA 样本实验记录](samples/prjna298330_sra.xml) · [完整项目清单](retry/prjna298330.tsv)

NCBI 的 3995 AND (Tran OR Rosenberg OR KRAS) 全字段检索有 3 条实验记录，其中 2 条是上述 3995 外显子，另一条是无关小鼠肠道数据。样本别名核查避免了仅凭数字命中认领数据。ENA 的样本标题模糊检索也包含大量无关数字，原始响应留在本地，公开审计采用能绑定患者和研究的实验记录。

目前 9d 的公开配对编码序列仍来自 AU2016323017B2 SEQ 28/29，β nt385 保留 N。缺少带原始克隆或同细胞对应的双链读段、RACE 克隆测序或色谱，无法补足该碱基或把这份专利编码序列绑定到晶体订购物料。3995 的外显子记录确实存在，但缺少完成已重排 αβ 克隆配对所需的证据。

## 扩大的独立结构筛查

本轮做了两条检索，并对成功响应核对哈希和元数据完整性。

| 检索 | 范围和结果 |
| --- | --- |
| RCSB 全库名称与标题 | 对实体描述或标题中 T cell receptor、T-cell receptor、TCR 做并集查询，返回 838 项沉积。补取此前库存外的 533 项元数据，与旧库存合并共 1,125 项 |
| STCRDab 全量快照 | 984 行、634 个 PDB，其中 class-I αβ 配对 354 行。用直接编号记录形成 155 组观察域查询，3 行因缺链编号或沉积对应不符排除 |
| 观察域发现检索 | 两条序列须分别出现在不同沉积实体中。得到 148 次跨研究且作者名无交集的命中，按 PDB 对去重为 88 对。不同长度片段可重复命中同一对，这些数字不代表独立受体数 |
| 完整变量区检索 | 沿用 44 个边界已核实的锚点，在全部 1,125 项中精确查双链。作者名无交集的跨研究命中为 7N1F–7RTR、4EI6–4ELK、4EI6–4ELM |
| 新候选官方组装 | 对 13 项候选下载官方生物组装及原始 PDB 的 COMPND 角色记录，以既有链角色和接触规则提取单元。13 项都未产生完整 class-I αβ TCR–pMHC 单元 |

[原始 RCSB 查询](rcsb_broad_query.json) · [RCSB 结果](broad/rcsb_broad_hits.json) · [STCRDab 快照](additional/stcrdab_summary.tsv) · [逐候选审计](audit.json)

13 项为 2BNU、2VLM、2XN9、2XNA、3QH3、3UTP、4GRM、7AMP、7AMQ、7AMR、7AMS、9NRJ、9NDS。游离 TCR、TCR–抗体、超抗原复合物均缺少既定目标所需的完整单元。9NRJ 是 TCR–T-cell engager 复合物，9NDS 是 mini-HLA–肽–TCR，后者缺 B2M，也被既定完整单元规则排除。没有补原子或构造补全链。

完整变量区命中的 4EI6 是游离 NKT TCR，4ELK 是另一研究的游离 NKT TCR，4ELM 为 CD1d–脂质体系。这两对保留作来源线索，未纳入当前 class-I pMHC 评估。7N1F–7RTR 仍受 7RTR β E28 缺界面重原子阻断，前轮官方组装复验记录继续适用。其余观察域命中需核实完整 N 端和构建体身份，片段相等没有触发准入。

这一轮检索使用现有完整锚点和 STCRDab class-I 观察域，覆盖范围仍受查询序列、数据库实体描述、已存编号和作者名标准化限制。还存在 3 个排除的编号行，未用推断补齐。查询结果没有证明公开数据库中不存在其他可用正对照。

## 未完成项及所需原件

| 未完成项 | 当前缺口 |
| --- | --- |
| Sim 2020 原始 SI | 文件被公开文章引用，下载路径尚未返回可核读 PDF，接口非 OA 限制和浏览器挑战已留档 |
| 两处差异具体引入步骤 | 缺带插入边界的 GenScript 设计或订单文件、具体晶体表达质粒的测序记录，或能联系到样品的原始 N 端鉴定 |
| 3995 原始 9d 配对读段 | 已定位的明确患者读段为外显子，缺原始双链克隆或同细胞关系，专利核酸仍有一个未定碱基 |
| 合格独立跨研究正对照 | 完整变量区命中仍受目标体系或原子完整性限制，新候选也未通过完整单元要求 |

训练 release 和评估协议尚未生成，独立验收集未访问。

复验命令如下，输出路径须尚不存在。

    PYTHONPATH=src:. python -m scripts.audit_public_construct_evidence --out /tmp/public-construct-evidence-replay.json

[保留源文件和哈希](retained_source_manifest.json)记录离线复验所用的本轮来源，旧库存沿用各自归档收据。[全部请求收据](retrieval_receipts.json)同时保留下载失败、非 PDF 响应及本地留存状态。
