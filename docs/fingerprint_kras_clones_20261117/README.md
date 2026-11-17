# KRAS 9a／9d 原始克隆与晶体构建体审计（2026-11-17）

两个原始配对克隆在 α 变量区第 10 位均为 Met，β 变量区首位均为 Asp。6ULR／6ULN 沉积对应位置分别为 Val 和 Met。差异已经定位在原始克隆与晶体构建体之间，具体引入步骤仍缺直接记录。[9a 专利，EP3494133B1 表 1、表 5 和序列表](https://patentimages.storage.googleapis.com/a6/63/45/76d90e3a55d9a0/EP3494133B1.pdf) · [9d 专利，AU2016323017B2 实例 2、3 及序列表](https://patentimages.storage.googleapis.com/db/74/cd/810734a40ea5c1/AU2016323017B2.pdf)

| 原始克隆 | α 原始／沉积成熟段起始 | β 原始／沉积成熟段起始 | 构建体对应 |
| --- | --- | --- | --- |
| 9a，4095 | LAKTTQPISM / LAKTTQPISV | DAGVTQSPTH / MAGVTQSPTH | EP3494133B1 SEQ 63、64，DNA 388、400 nt。翻译后对应 6ULR D/E，两条 CDR3 均精确匹配 |
| 9d，3995 | LAKTTQPISM / LAKTTQPISV | DAGVTQSPTH / MAGVTQSPTH | AU2016323017B2 SEQ 9、10，蛋白 129、132 aa，含各自前导区。对应 6ULN D/E，两条 CDR3 均精确匹配 |

9a 的 EP 序列逐碱基与同族 AU2023233125B2 的机器可读清单交叉确认。9d 的 AU 公告附原始 ST.25 序列表，EP3901169A1 实例 2、3 把这对变量区追溯至 3995 的 KRAS 反应性 T 细胞，并说明随后接小鼠恒定区表达。两份专利与 Sim 等论文 Table 1 的 CDR3 和克隆对应吻合。源页、序列表转录和 SHA256 见 [patent_listing_receipts.json](patent_listing_receipts.json)。[9a 同族清单](https://patents.google.com/patent/AU2023233125B2/en) · [9d 同族专利](https://patents.google.com/patent/EP3901169A1/en) · [论文勘误版 Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC7959486/)

两条 D 链的 SEQRES 在成熟段之前带 MAG，差异 Val 对应 PDB 的 D13，两条结构都有其原子。STCRDab 的 α 编号比专利成熟段短一位，末端 Ser 实际紧邻并存在于 SEQRES。延至相同长度后两个 α 都只剩 M10V。两个 β 的编号段只剩 D1M，首位 Met 在 SEQRES 里，晶体 ATOM 未建模该位。Sim 等方法记录胞外段合成、NdeI/XhoI 克隆至 pET30a，在 BL21(DE3) 表达。β 首位 Met 与使用翻译起始位点的构建解释相容。实际订购序列、晶体质粒或首端蛋白测序仍缺，因此 α M10V 的具体环节、β 是有意替换还是构建起点处理，仍未证实。原论文 SI 附件在 PMC、出版社及备选入口均返回访问错误或 HTML，相关请求留档。[Sim 等原文方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC7293613/)

## 9a 单细胞读段核查

从公开 GSE136394 的 4095 肿瘤样本 GSM4047960／SRR10039337 读取前 480,000 个 spot（总计 28,317,417）。使用 NCBI SRA Normalized VDB 的原始 READ 和 QUALITY 列保留技术 R1。原提交 FASTQ 的匿名 S3 访问返回 Requester Pays 拒绝，ENA 单文件 FASTQ 前缀只有 R2。原始读段和导出留本地，公开归档保留固定 accession、行范围、哈希、参数和复验代码。[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136394) · [SRA](https://www.ncbi.nlm.nih.gov/sra/SRR10039337)

343,246 个 spot 通过条形码和 UMI Q20。以专利核酸精确匹配两链 CDR3，20 个原始 16 碱基条形码同时出现两链。13 个条形码还各自读到了 α ATG 和 β GAC，待查密码子逐碱基 Q30、两侧锚序列 Q20。它们支持 9a 原始克隆的 Met/Asp。在这 13 个条形码中，未见 PDB 型 α Val 或 β Met。整个筛查窗口也存在零星其他密码子，已按实际观察保留。两个 CDR3 都各有至少两个 UMI 的条形码为 0。这里是同条形码分段读证据，未完成细胞纠错、双细胞识别或全长配对组装，亦未验证整份 SRA 文件的 MD5。

9d 找到了来自原始分选克隆的已公开配对**蛋白**序列。当前未找到 3995 的原始配对核苷酸读段或克隆色谱。GSE136394 和早期论文 PRJNA342632 的 ENA 公开 run 清单分别有 128 和 49 项，标题均无 3995。这个检索范围不覆盖未标样本编号的其他项目或未开放材料。

## 独立跨研究正对照

11 月 17 日 RCSB 新发布和修订查询均返回空结果。对上轮 317 项 class-I 目录作已沉积 αβ 全构建体筛查，得到 227 个唯一链对组，其中不同论文且作者无交集的完全相同沉积链对为 0。这项筛查依赖实体描述和全构建体相等，恒定区或标签变化的情形需用变量区搜索。另以本轮四条专利的成熟区查 RCSB，两个 α 查询均无命中，β 的局部命中仅为已知 6ULR、6ULN、7SU9。β 局部相似不抵消已证实的 N 端差异。6ULR／7SU9 仍有四位共同作者，7SU9 官方组装仍缺固定规则可分配的完整接触单元。[前轮结构和原子 QC](../fingerprint_weekly_controls_20261116/README.md) · [本轮机器审计](audit_final.json)

6ULN／6ULR 的官方生物组装原子 QC 结论继续成立。原始克隆证据解释了差异出现在构建体中，现行完整变量区规则仍拒收这两份结构。新增合格跨研究正对照与同 V 天然对照均为 0，未冻结评估，独立验收集未访问。

离线复验两链的专利—坐标对齐、公共项目清单和独立结构筛查可运行以下命令。如本地有 [read_receipt.json](read_receipt.json) 中的 SRA 导出，可追加 `--reads PATH` 精确重放读段证据。

    PYTHONPATH=src:. python -m scripts.audit_kras_clone_provenance --out /tmp/kras-clone-replay.json
