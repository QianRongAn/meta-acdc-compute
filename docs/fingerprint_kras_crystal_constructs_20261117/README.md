# KRAS 9a／9d 晶体构建体与独立正对照续查（2026-11-17）

找到 9d 的公开配对核苷酸序列，确认两个 N 端位置在专利序列中分别为 Met 和 Asp。9d 原始核苷酸读段和 6ULN／6ULR 晶体表达质粒的订购或测序记录仍缺；现有证据不能把两个差异精确归到某次设计、合成、克隆、翻译处理或沉积。新增合格独立结构正对照为 0。

## 9d 配对核苷酸

[AU2016323017B2 公告](https://patentimages.storage.googleapis.com/db/74/cd/810734a40ea5c1/AU2016323017B2.pdf) 的 [0054]（PDF 第 23 页）列出配对的 α SEQ ID 28、β SEQ ID 29；PDF 第 82 页附完整 ST.25 记录。以同一公告的蛋白 SEQ ID 9、10 逐密码子复核。专利实例 2、3 追溯受体至 3995 的 KRAS 反应性 T 细胞，并说明接小鼠恒定区做功能重表达；[Tran 等人](https://pmc.ncbi.nlm.nih.gov/articles/PMC5178827/)和[结构论文勘误中的克隆表](https://pmc.ncbi.nlm.nih.gov/articles/PMC7959486/)提供克隆与识别背景。

| 链 | 专利公开核苷酸 | 翻译核验 |
| --- | --- | --- |
| α | SEQ 28，387 nt，含前导区；与 9a SEQ 63 的前 387 nt 完全相同 | 129 aa 与 SEQ 9 完全一致；成熟段第 10 位为 ATG → Met |
| β | SEQ 29，396 nt，含前导区；第 385 位为 N，落在第 129 个密码子 NTA 的首碱基 | C→CTA、T→TTA 均译为 Leu，与 132 aa 的 SEQ 10 完全一致；A→ATA 为 Ile、G→GTA 为 Val，均与蛋白清单不符。成熟段第 1 位为 GAC → Asp，不受此歧义影响 |

两条记录均标为 synthetic。它们公开了与 9d 相符的配对编码序列，缺少能把每个碱基直接追溯到 3995 原始克隆的读段或色谱。FASTA 保留原始 N，未择取某个同义密码子。这对序列不是 6ULN 的 pET30a 质粒记录。[FASTA](9d_disclosed_nucleotides.fasta) · [源页转录](au_page82.txt) · [比对结果](audit.json)

## 6ULN／6ULR：能定位到哪里

- **α M10V：** 9a 与 9d 专利变量区均为 Met；两份 PDB 的 D13 都有建模的 Val，成熟段其余位置相同。首端还声明了额外的 MAG。M10V 在成熟段内部，不能由翻译起始 Met 的加入或去除解释。原始克隆到沉积结构之间的具体环节待核实。
- **β D1M：** 两份专利均为 Asp，PDB E1 声明为 Met，但两份结构都没有 E1 原子。2020 原文记载胞外区由 GenScript 合成，经 NdeI/XhoI 放入 pET30a、在 BL21(DE3) 表达。NdeI 位点可提供起始 ATG，所以表达起点的处理是一个解释；没有插入片段边界、原始设计或蛋白 N 端测序，不能区分设计时替换、使用翻译起点、蛋白末端加工与仅在沉积序列中声明 Met。[Sim 2020 方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC7293613/) · [6ULN](https://www.rcsb.org/structure/6ULN) · [6ULR](https://www.rcsb.org/structure/6ULR)
- **沉积和质粒库：** 两份 mmCIF 的 TCR 实体没有填突变说明；唯一列出的序列差异条目是 KRAS 肽 G12D。空字段不能推断 TCR 非工程构建。Addgene 以 TCR9a 和 TCR9d 原样检索各显示 0 项，只能说明该目录及查询没有命中。[沉积字段摘录](deposition_annotations.json) · [访问状态与哈希](retrieval_receipts.json)
- **7SU9 和补件：** [eLife 2022](https://elifesciences.org/articles/75670) 的试剂表把 9a 两条可溶链的来源都指向 Sim 2020，没有列出质粒编号。文中记录胞外 α 1–208、β 1–244 的合成与克隆体系，并明确用 6ULR 作 7SU9 的分子置换搜索模型。7SU9 与 6ULR 两条沉积 TCR 全序列一致，仍有共同作者，也没有通过既定完整接触单元门。7SU9 的 [补充文件 1](https://cdn.elifesciences.org/articles/75670/elife-75670-supp1-v2.docx) 是晶体数据收集和精修表，不含表达质粒或订购序列。试剂来源指向同一研究不等于同一实体质粒已获证实。
- **尚未取得的直接材料：** 需能绑定 6ULN／6ULR 样品的 GenScript 原始订购文件、带插入片段和接头的实物质粒测序、或蛋白末端测序记录。Sim 2020 SI 和 Tran 2016 附件的 PMC／期刊链接本轮仍返回挑战页或拒绝，虽有 HTTP 200 也不是 PDF；PMC 公共对象目录未提供这两份附件。不能据此声称两个附件一定没有构建信息。尚未找到 3995 原始双链读段、克隆色谱或能绑定原始克隆的原始组装记录；此前 4095 的单细胞读段不代表 3995。

## 独立结构候选和放行

对既有 44 个边界已核实的完整 αβ 变量区锚点，在 592 项存档 RCSB 元数据中直接做双链精确搜索：允许恒定区或标签不同，不依赖实体名称；α 和 β 必须在不同实体中各出现一次。得到 17 个非自身锚点—沉积命中，8 个同研究、6 个跨研究且作者重合、2 个来源未确定、1 个跨研究且作者名无交集。这是有界检索，实体匹配不能自行证明两条链配对或原子 QC 合格。同日公开发布／修订的前轮检索收据见[上一轮审计](../fingerprint_kras_clones_20261117/audit_final.json)。

唯一作者名无交集的命中是此前已确认身份的 YLQ7 **7N1F／7RTR**。本轮用两份 [7N1F](https://www.rcsb.org/structure/7N1F)／[7RTR](https://www.rcsb.org/structure/7RTR) 官方生物组装复验，逐一核对沉积原子身份并沿用固定原子和图构建规则。7N1F 通过，7RTR 的 β 链 E28 仍缺 CD、CG、OE1、OE2；没有补原子。6ULN／6ULR 的官方组装原子 QC 继续成立，原始克隆与沉积变量区的两个差异仍阻止按原始身份纳入。未更新数据准入或划分，未冻结评估，未训练，未访问独立验收集。

离线重放：

    PYTHONPATH=src:. python -m scripts.audit_kras_crystal_constructs --out /tmp/kras-crystal-constructs-replay.json
