# SOP — Meta-ACDC 计算管线复现手册(KN-16 交付物)

> 目标:从零复现本文档引用的全部计算结果(数字见 `docs/benchmarks.md`)。
> 参考主机:Ubuntu + GTX 1050 Ti 4GB / 8GB RAM;CPU-only 亦可(训练慢 ~3-5×)。
> 回归测试:`.venv/bin/python -m unittest discover -s tests -v`(无数据的用例自动跳过)。
> **一键健康检查:`scripts/verify.sh`**(导入 + 测试 + manifest/CIF 一致性 +
> 稳定性判决 + 临床扫描;数据缺失自动 SKIP)。

## 0. 环境两条路径

**A. 裸机(推荐,支持 GPU)**
```bash
bash scripts/setup_env.sh          # uv venv + torch 2.9.1+cu126 + 科学栈
.venv/bin/pip install -e . --no-deps
.venv/bin/pip install mhcflurry    # 预筛基准(2.2.1,torch 后端)
```
注意:Python ≥ 3.13 时 mhcflurry 下载器会因 `pipes` 移除而崩,用 shlex shim:
```bash
.venv/bin/python -c "
import sys, types, shlex
m = types.ModuleType('pipes'); m.quote = shlex.quote; sys.modules['pipes'] = m
sys.argv = ['x', 'fetch', 'models_class1_presentation']
from mhcflurry.downloads_command import run; run()"
```

**B. Docker(CPU)**
```bash
docker build -t meta-acdc .
docker run --rm -v $PWD:/work meta-acdc python -c "import meta_acdc; print('OK')"
```

## 1. 数据管线(KN-3)

```bash
.venv/bin/python src/meta_acdc/data/download_vdjdb.py     # -> data/raw/vdjdb/
.venv/bin/python src/meta_acdc/data/download_iedb.py      # -> data/raw/iedb/
.venv/bin/python src/meta_acdc/data/ingest_vdjdb.py       # -> data/processed/vdjdb.clean.tsv (199,289)
.venv/bin/python src/meta_acdc/data/ingest_iedb.py        # -> data/processed/iedb.clean.tsv (2.3M)
.venv/bin/python src/meta_acdc/data/stats_report.py       # -> docs/dataset-report.md
.venv/bin/python src/meta_acdc/data/clinical_cases.py     # -> clinical_gold_standard.tsv
```

结构集:STCRDab 抓取 + 自动链分类(见 `structure/fetch_structures.py`;
`data/raw/structures/complexes.tsv` 为 291 个单复合物清单)。

## 2. 结构 ↔ 特异性资源(KN-4+)

```bash
.venv/bin/python src/meta_acdc/structure/vdjdb_structure_map.py
# 281/291 复合物映射 VDJdb;2,309 CDR3 命中(791 同表位 / 1,517 异表位)
```

## 3. EGNN 训练与诱饵协议(KN-4,论文 1 核心)

```bash
.venv/bin/python src/meta_acdc/models/dataset.py     # 界面图 + 位移/嫁接诱饵
.venv/bin/python src/meta_acdc/models/train_egnn.py  # v9.1:嫁接 0.933 / 位移 0.980
```
诱饵协议(防泄漏核心):位移 = 肽链刚性 +6Å;嫁接 = Kabsch 侧链移植
(骨架/节点集/肽标志不变);节点集 pin;供体序列去重;PDB 分组切分。
消融与对照谱系见 `docs/benchmarks.md`。

## 4. 元学习基线(KN-6)

```bash
.venv/bin/python src/meta_acdc/models/train_protonet.py   # 0.635(留出表位)
.venv/bin/python src/meta_acdc/models/train_maml.py       # MAML++ 0.625
```

## 5. 主动学习模拟(KN-7 / KN-7+)

```bash
.venv/bin/python src/meta_acdc/active_learning/simulate_al.py --mode eig --seeds 2
.venv/bin/python src/meta_acdc/active_learning/simulate_al.py --mode random --seeds 2
# KN-7+:ACDC 式预筛肽库 + 批大小消融
.venv/bin/python src/meta_acdc/active_learning/simulate_al.py \
    --pool-mode prefiltered --rounds 6 --batch 1500 --seeds 3
.venv/bin/python src/meta_acdc/active_learning/simulate_al.py --ablation
```

## 6. MHC 预筛基准(计划书 Module 0 假设)

```bash
# 本地(MHCflurry,无排队):top-2% 阈值 0.56 召回 / ~50× 压缩
.venv/bin/python src/meta_acdc/data/mhcflurry_bench.py --n 50
# NetMHCpan-4.2(在线,受 DTU 队列影响):
.venv/bin/python src/meta_acdc/data/netmhcpan_bench.py --n 50
```

## 7. 蛋白组 5 万候选清单(KN-8 计算侧交付物)

```bash
curl -sL 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(organism_id:9606)%20AND%20(reviewed:true)' \
    -o data/raw/human_sprot.fasta
.venv/bin/python src/meta_acdc/active_learning/select_candidates.py --stage prefilter
.venv/bin/python src/meta_acdc/active_learning/select_candidates.py --stage rank \
    --tcr CASSLGRYNEQFF   # Kimmtrak/gp100 天然前体;A6 TCR: CAVTTDSWGKLQF
# -> data/processed/acdc_pool5m.tsv.gz(500 万预筛池)
# -> data/processed/kn8_top50k_*.tsv(EIG 前 5 万,含蛋白来源定位)
```

## 8. AF3 结构打分与排名(KN-5,湿实验对接侧)

```bash
# AF3 预测经官方网页版提交(条款合规,人工操作);CIF 放入 data/raw/af3_predictions/
.venv/bin/python src/meta_acdc/structure/score_predictions.py \
    --model data/processed/egnn_dataset.model.pt \
    --out data/processed/prediction_scores.tsv
# 域适应训练(AF3 天然正样本,提交后;--af3 直接吃 CIF,全诱饵协议):
# ⚠️ 天然提交务必带 --unique-fallback:否则无精确链对匹配的 job 会退化为
#    candidates[0],造成不同 job 塌缩到同一 pdb + 把天然误标成交叉反应
.venv/bin/python src/meta_acdc/structure/import_af3.py \
    --src "$FOLDS_DIR" --list data/processed/kn5_submission_list.tsv \
    --out data/raw/af3_native_predictions --unique-fallback \
    --report data/processed/af3_native/import_report.tsv
# 权威做法:用 af3_native*/ 的天然 FASTA 链指纹回贴真实 job_id(见
# docs/benchmarks.md §2026-09-23 真实 native DA),再重建 manifest.tsv
# 一键:导入天然 + 把交叉反应候选路由到 af3_predictions/(不丢弃)
.venv/bin/python src/meta_acdc/structure/native_manifest.py \
    --src "$FOLDS_DIR" --out data/raw/af3_native_predictions \
    --manifest data/processed/af3_native/manifest.tsv \
    --report data/processed/af3_native/import_report.tsv \
    --candidate-out data/raw/af3_predictions \
    --list data/processed/kn5_submission_list.tsv
.venv/bin/python src/meta_acdc/models/train_domain_adapt.py \
    --af3 data/raw/af3_native_predictions \
    --af3-manifest data/processed/af3_native/manifest.tsv --epochs 150
# ⚠️ 规模要求:候选跨实例判决过 0.5 门需 ~144 个天然正样本(63→塌缩,
#    89→0.479,144→0.729 STABLE);正样本越多越稳,务必先攒够再报告排名
# 下游:临床金标准结构扫描(KN-11,稳定模型上)
.venv/bin/python src/meta_acdc/structure/clinical_scan.py \
    --scores data/processed/prediction_scores_ensemble.tsv \
    --out data/processed/kn11_clinical_scan.tsv
# KN-8 结构重排:给定治疗 TCR 结构 + 候选肽表 -> AF3 5 链提交 FASTA
.venv/bin/python src/meta_acdc/structure/af3_scan_batch.py \
    --pdb 5brz --peptides data/processed/kn8_top50k_CAVTTDSW.tsv \
    --n 200 --out data/processed/af3_scan_5brz
# 提交后按上面流程导入/打分,用稳定 DA 模型对候选肽重排(替代弱序列代理)
# 域适应模型(pLDDT 置零训练+打分一致):
.venv/bin/python src/meta_acdc/structure/score_predictions.py \
    --model data/processed/da_seed0.model.pt --mask-plddt --out ...
.venv/bin/python src/meta_acdc/structure/da_stability.py --scores \
    data/processed/da_scores_seed{0,1,2}_ensemble.tsv   # 跨实例稳定性判决
# 重提感知导入 + AF3 vs EGNN 方差分解(2026-09-23):
.venv/bin/python src/meta_acdc/structure/import_af3_instances.py \
    --src "/mnt/c/.../folds_YYYY_MM_DD_HH_MM"      # -> data/raw/af3_instances/{A,B1,B2,B3}
.venv/bin/python src/meta_acdc/structure/af3_variance.py \
    --out data/processed/af3_resubmission_variance.tsv
.venv/bin/python src/meta_acdc/structure/af3_rmsd.py \
    --out data/processed/af3_resubmission_rmsd.tsv   # 模型无关肽 RMSD
.venv/bin/python src/meta_acdc/structure/af3_qc.py \
    --out data/processed/af3_qc.tsv                  # ipTM/链完整性可靠性闸门
```

# 域适应正样本(AF3 预测风格的天然复合物),待用户网页提交:
.venv/bin/python src/meta_acdc/structure/af3_batch.py \
    --evidence same --n 1000 --out data/processed/af3_native
# -> 174 个天然 5 链 FASTA,提交后按 §8 导入并混入 train_domain_adapt.py
```

## 9. 风险看板(KN-13)

```bash
.venv/bin/python src/meta_acdc/dashboard/server.py --port 8000
# http://localhost:8000 — 热度图 / 排名 / 金标准案例 / 交叉反应证据 / 诚信横幅
```

## 10. 关键陷阱记录(踩过的坑)

- Python ≥ 3.13:`pipes` 模块移除,mhcflurry 下载器需 shlex shim(见 §0)
- `SimpleHTTPRequestHandler`(Py ≥ 3.12):类属性 `directory` 失效,必须
  构造器传参(`partial(Handler, directory=...)`)
- mhcflurry `Class1PresentationPredictor.predict` 的 alleles 是基因型语义
  (≤6),单等位基因传一元素列表
- v9 首训 0.988 系诱饵 atoms 字段泄漏(已回退重做)——"看似突破先自查泄漏"
- AF3 网页版五链任务必须用多链 FASTA;高置信任务重跑很稳(ipTM Δ≤0.07),
  但低 ipTM(~0.5,JM22 类)任务 AF3 自身会飘(Δ0.21、肽 RMSD 达 11 Å)——
  **用 ipTM 做可靠性闸门**;旧的"重提方差 0.94"系打分模型侧假象
  (2026-09-23 修正);**不可信的是单实例 EGNN 打分,不是单次 AF3 提交**
- 任何 OOD 打分结论必须过跨实例稳定性检验(Rashomon 效应)
- **pdb 标签不可信**:`kn5_submission_list.tsv` 的 `tcr_a/tcr_b` 两列**写反**
  (tcr_a=β、tcr_b=α)。导入务必按**链对**匹配(`import_af3.tcr_pair` +
  `matches_chain_pair`),并在导入后跑
  `relabel_instances_by_chains.py`(按 CIF 链指纹回贴真 pdb)核对;
  2026-09-23 曾因此把 48/73 条结构标签弄错(1ao7↔1qrn、1g6r/1mwa→1tcr、
  1oga→2vlj、2bnr→2f53)

## 湿实验接口(不在计算范围,技术路线保留)

ACDC 平台对接点:(a) 50k 候选清单 = §7 输出(KN-8,用户订购);
(b) FACS/NGS 回读数据 → 替换 `simulate_al.py` 的 oracle 即成真闭环;
(c) 模型侧重训零改动(见 `docs/milestones/KEY-NODES.md` 技术路线表)。
