# 数据集卡片:Meta-ACDC 策展数据集(Deliverable 5)

> 本文件按 "Datasheets for Datasets" 风格记录本项目策展/派生的数据集:
> **来源、处理流程、统计、用途、局限**。原始第三方数据遵循其各自许可;
> 本仓库只分发**派生索引与小体积策展文件**(大文件由脚本复现,见
> [`SOP.md`](SOP.md))。

## 0. 总览

| 名称 | 文件 | 规模 | 是否入库 |
|---|---|---|---|
| VDJdb 清洗表 | `data/processed/vdjdb.clean.tsv` | 199,289 对 | 否(脚本复现) |
| IEDB 清洗表 | `data/processed/iedb.clean.tsv` | 2,302,095 条 | 否(脚本复现) |
| **结构↔特异性图谱** | `data/processed/structure_vdjdb_map.tsv` | 2,308 命中 | **是** |
| KN-5 提交清单 | `data/processed/kn5_submission_list.tsv` | 2,251 对 | **是** |
| 临床金标准 | `data/processed/clinical_gold_standard.tsv` | 3 案例 | **是** |
| 天然 AF3 正样本清单 | `data/processed/af3_native/manifest.tsv` | 144 | **是** |
| AF3 提交批(天然/候选) | `data/processed/af3_native*/`、`af3_batch1|2/` | 数百 FASTA | **是** |
| 交叉反应图谱 | `data/processed/structure_vdjdb_map.tsv` | 179/281 TCR | **是** |

---

## 1. 来源(provenance)

| 数据 | 版本/日期 | 来源 | 许可 |
|---|---|---|---|
| VDJdb | 2026-06-03 | vdjdb.cor.ac.ba | 见上游 |
| IEDB `mhc_ligand_full` | 2026 | iedb.org | 见上游 |
| STCRDab 结构 | 2026 | opig.stats.ox.ac.uk | 见上游 |
| UniProt Swiss-Prot(人) | 2026 | uniprot.org | CC-BY 4.0 |
| AlphaFold3 预测 | 2026-09 | alphafoldserver.com(网页版,条款合规,手动提交) | AF3 输出条款 |

## 2. 处理流程

1. **VDJdb**:去重 → 标准化 → `ingest_vdjdb.py` → `vdjdb.clean.tsv`
   (199,289 对;标签 0/1;score ∈ [0,3])。
2. **IEDB**:过滤 → 标准化 → `ingest_iedb.py` → `iedb.clean.tsv`
   (2.3M;544 等位基因)。
3. **结构集**:STCRDab 抓取 + **内容型链角色分类**(`classify_chains`,
   解决 PDB 链字母不一致);多复合物拆分 291 → 423 训练复合物。
4. **结构↔特异性映射**(`vdjdb_structure_map.py`):把每个复合物的
   TCR CDR3 与 VDJdb 记录匹配 → `structure_vdjdb_map.tsv`,标注
   `match = same`(同表位)或 `different`(**异表位 = 交叉反应证据**)。
5. **KN-5 清单**(`kn5_candidates.py`):从映射生成 (TCR α/β, 肽, MHC)
   提交对,优先交叉反应。
6. **AF3 导入**(`native_manifest.py` / `import_af3.py`):**链指纹匹配**
   天然提交 FASTA,交叉反应候选路由到独立目录;`relabel_instances_by_chains.py`
   按晶体链指纹回贴真实 pdb(修复 tcr_a/tcr_b 列写反导致的 48/73 错标)。

## 3. 统计

### 3.1 结构↔特异性图谱(`structure_vdjdb_map.tsv`)

- 281 个复合物映射成功(291 中)
- 2,308 条 CDR3 命中:**791 同表位 / 1,517 异表位(交叉反应)**
- **179/281 TCR 有交叉反应记录**;每 TCR 异表位中位 2、最大 8
- 列:`pdb, chain, cdr3, vdjdb_epitope, pdb_peptide, match`

### 3.2 VDJdb 清洗表

- 178,974 独特 CDR3-β;2,068 独特肽;170 MHC 等位基因
- 标签:0(92.9%)/ 1(7.1%);CDR3 长度众数 14–15

### 3.3 临床金标准(`clinical_gold_standard.tsv`)

| TCR | 靶肽 | 脱靶 | 致死 | 证据 |
|---|---|---|---|---|
| A3A(MAGE-A3,HLA-A*01) | EVDPIGHLY | ESDPIVAQY;ILAKFLHWL;KVAKELVHFL | 是 | PMID:23999400/05 |
| A2-restricted MAGE-A3 | KVAELVHFL | KMVELVHFL;KMAELVHFL | 是 | PMID:23470321 |
| Kimmtrak(gp100) | YLEPGPVTA | YLEPGPVTV;YLEPGPVTL | 否 | PMID:34525232 |

## 4. 用途(intended use)

- 训练/评估 TCR-pMHC 化学相容性模型(诱饵任务、排名、安全扫描)
- 提供"结构 ↔ 免疫学验证特异性"的连接(交叉反应图谱)
- 生成 AF3 提交批(天然正样本 / 交叉反应候选)

## 5. 局限与伦理(limitations)

- **标签噪声**:VDJdb score=0 多为"未观察到应答"而非系统验证阴性;
  score=1 为弱阳性 → 报告须谨慎。
- **覆盖偏倚**:STCRDab 偏向研究充分的 TCR(如 A6/MAGE 家族),
  交叉反应图谱不能代表全 TCR 组库。
- **pdb 标签**:提交清单 `tcr_a/tcr_b` 列曾写反(已按链指纹修复);
  下游用链对匹配,勿依赖单一列。
- **AF3 质量**:TCR-pMHC DockQ ~0.5–0.57;低 ipTM(<0.75)结构不可信
  (见 `af3_qc.tsv` 可靠性闸门)。
- **湿实验数据缺失**:本数据集**不含** ACDC 湿实验验证(计划书
  Deliverable 5 的"百万级 ACDC 验证条目"需实验室产出);当前为
  **公开数据策展 + 结构派生**版本。
- **二次分发**:原始 VDJdb/IEDB 数据请从其官方获取;本仓库只分发
  派生索引。

## 6. 复现

```bash
python src/meta_acdc/data/download_vdjdb.py && python src/meta_acdc/data/ingest_vdjdb.py
python src/meta_acdc/data/download_iedb.py  && python src/meta_acdc/data/ingest_iedb.py
python src/meta_acdc/structure/vdjdb_structure_map.py
python src/meta_acdc/structure/kn5_candidates.py
python src/meta_acdc/data/clinical_cases.py
```

统计报告:`python src/meta_acdc/data/stats_report.py` → `docs/dataset-report.md`。

## 7. 引用

使用本数据集请引用计划书与相应上游数据库(VDJdb / IEDB / STCRDab /
UniProt)。见 [`../README.md`](../README.md) 与
[`literature/references-*.bib`](literature/)。
