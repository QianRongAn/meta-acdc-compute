# Meta-TCR API 文档(Deliverable 3)

> 本包(`meta_acdc`)是 Meta-ACDC 计算套件的可复用接口:**AlphaFold3 结构
> 接口 + EGNN 编码器 + 主动学习采样**。既可命令行调用(见文末 CLI),
> 也可作为库导入。所有函数均有类型注解与 docstring。

安装:

```bash
pip install -e .              # 仅核心(无 torch)
pip install -e ".[ml]"        # 含 torch 训练栈
```

---

## 1. 结构接口(`meta_acdc.structure`)

### 1.1 `graph` — 界面图构建(核心)

```python
from meta_acdc.structure.graph import (
    parse_pdb, classify_chains, build_interface_graph_from_residues,
    TCR, MHC, B2M, PEPTIDE, InterfaceGraph, Residue,
)
```

- `parse_pdb(path) -> list[Residue]`:解析 PDB,每个残基含全部重原子坐标。
- `classify_chains(residues) -> dict[chain, role]`:**基于内容**判定链角色
  (TCR / PEPTIDE / MHC / B2M),解决 PDB 链字母不一致。
- `build_interface_graph_from_residues(residues, interface_radius=10.0,
  edge_cutoff=8.0, forced_nodes=None) -> InterfaceGraph`:
  构建界面图。节点 = 残基;边 = 质心距离 < `edge_cutoff`。
  `forced_nodes` 固定节点集(诱饵生成防泄漏用)。
- `InterfaceGraph`:`.n_nodes`, `.node_features`(26 维),
  `.node_coords`, `.edge_index`, `.edge_features`(16 维 = 12 RBF + 4 接触),
  `.residue_info` `(chain, resname, resid)`, `.node_keys`。

### 1.2 `cif` — AlphaFold3 mmCIF 解析

```python
from meta_acdc.structure.cif import parse_cif, parse_all_models
residues = parse_cif("job_model_0.cif", model="1")   # -> list[Residue]
```

pLDDT 落在 `Residue.plddt`(来自 B-factor 列)。与 PDB 管线共用同一
`Residue` 结构,下游无需改动。

### 1.3 `score_predictions` — EGNN 打分

```python
from meta_acdc.structure.score_predictions import score_graph
s = score_graph(net, graph, device="cuda", mask_plddt=True, temperature=1.0)
```

- `mask_plddt`:域适应模型训练时 pLDDT 置零,打分须一致。
- `temperature`:来自 `models.calibrate` 的温度,用于标定概率。
- CLI:`meta-acdc-score --cifs <dir> --model <ckpt> --mask-plddt --out <tsv>`
  (可重复 `--model` 做集成;同时输出 `*_ensemble.tsv`)。

### 1.4 `import_af3` — 提交清单匹配

```python
from meta_acdc.structure.import_af3 import (
    read_job_sequences, tcr_pair, matches_chain_pair, load_manifest,
)
```

- `tcr_pair(seqs, peptide) -> (alpha, beta)`:按保守基序(如 `YFCAS`)
  区分 α/β,顺序无关。
- `matches_chain_pair(alpha, beta, cand)`:按**链对集合**匹配(兼容提交
  清单 tcr_a/tcr_b 列写反)。
- CLI:`meta-acdc-native-manifest --src <folds> --out <native_dir>
  --manifest <tsv> --candidate-out <cand_dir>`(链指纹匹配天然正样本,
  并把交叉反应候选路由到独立目录)。

### 1.5 `da_stability` — 跨实例稳定性判决门

```python
# CLI
meta-acdc-stability --scores seed0_ensemble.tsv seed1_ensemble.tsv seed2_ensemble.tsv
```

对多个模型实例的候选分求两两 Spearman,自动**剔除退化实例**(常数输出),
给出 mean r + bootstrap 95% CI + 判决(STABLE 当 mean r ≥ 0.5)。

### 1.6 `ensemble_mean` — 固定多实例集成

```python
from meta_acdc.structure.ensemble_mean import load_seed   # 或 CLI
```

对多个 `*_ensemble.tsv` 求均值(自动排除退化 seed),输出报告用集成分。

### 1.7 `clinical_scan` — 临床金标准扫描(KN-11)

```python
meta-acdc-clinical --scores prediction_scores_ensemble.tsv \
    --clinical clinical_gold_standard.tsv --out kn11_clinical_scan.tsv
```

对每个临床案例,按家族内排名判读靶肽与脱靶;`--margin 0.1` 定义
"与靶同水平"的致命模拟判定。

### 1.8 `attribute` / `alanine_scan` — 可解释性

```python
meta-acdc-attribute --cif X.cif --model da_seed0.model.pt --mask-plddt
meta-acdc-mutate    --cif X.cif --model da_seed0.model.pt --mask-plddt \
                    --mutant TRP --cumulative
```

- `attribute`:梯度归因到边的接触通道,输出 top 残基对(含肽相关接触)。
- `alanine_scan`:逐位突变(ALA=保守对照,TRP/ARG=冲突探针);
  `--cumulative` 给出贪心累积突变曲线(分布式 vs 热点)。

### 1.9 `af3_batch` / `af3_scan_batch` — 提交批生成

```python
from meta_acdc.structure.af3_batch import build_job, fasta_block
```

- `af3_batch`:从 KN-5 提交清单生成 5 链 FASTA(每 job)。
- `af3_scan_batch`:给定治疗 TCR 结构 + 候选肽表 → 5 链 FASTA
  (KN-8 结构重排用)。

---

## 2. 模型(`meta_acdc.models`)

### 2.1 `egnn` — 等变图神经网络

```python
from meta_acdc.models.egnn import EGNN, EGCL
net = EGNN(node_dim=26, edge_dim=16, depth=6, hidden=128)
logit = net(h, x, edge_index, edge_attr, batch)          # (n_graphs, 1)
feat  = net(h, x, edge_index, edge_attr, batch, embed_only=True)
```

- E(3) 等变;读出 = 全局均值 + 全局最大 + **肽节点均值** + 坐标范数。

### 2.2 `dataset` — 诱饵构造(可复现性核心)

```python
from meta_acdc.models.dataset import (
    make_decoy, make_displaced_peptide_decoy, build_dataset,
)
```

- `make_displaced_peptide_decoy(residues, shift=6.0)`:肽刚性平移,几何破坏。
- `make_decoy(residues, foreign_peptide)`:Kabsch 骨干对齐移植外来侧链
  ——骨架/节点集/肽标志不变,仅化学改变。
- 防泄漏:节点集 pin(`forced_nodes`)、供体序列去重、PDB 分组切分。

### 2.3 `train_egnn` / `train_domain_adapt`

```bash
meta-acdc-train-egnn            # 晶体诱饵任务(基线)
meta-acdc-train-da --af3 <native_dir> \
    --af3-manifest <manifest.tsv> --seeds 0 1 2 3 4 --epochs 150
```

`train_domain_adapt` 保存 `da_seed{k}.model.pt`(含 `val_auroc` 元数据)。

### 2.4 `calibrate` / `eval_da_crystal`

```bash
python -m meta_acdc.models.calibrate --model da_seed0.model.pt --out T.json
python -m meta_acdc.models.eval_da_crystal --da da_seed0.model.pt ...
```

- `calibrate`:温度标定(最小化 BCE),报告 ECE。
- `eval_da_crystal`:域适应是否损害晶体诱饵任务的回归检查。

---

## 3. 主动学习(`meta_acdc.active_learning`)

### 3.1 `acquisition`

```python
from meta_acdc.active_learning.acquisition import expected_information_gain
eig = expected_information_gain(probs)   # probs: (n_members, n_candidates)
```

### 3.2 `simulate_al`

```bash
python -m meta_acdc.active_learning.simulate_al --mode eig --seeds 2
python -m meta_acdc.active_learning.simulate_al --ablation
```

EIG/熵/方差/随机采集函数 + ε-greedy;支持预筛池模式与批大小消融。

### 3.3 `select_candidates`

```bash
meta-acdc-select --stage prefilter          # 蛋白组 -> 5M 池
meta-acdc-select --stage rank --tcr <CDR3>  # EIG 排序 -> top-50k
```

`--max-aa-fraction 0.5` 过滤低复杂度伪影。

---

## 4. 看板(`meta_acdc.dashboard`)

```bash
meta-acdc-dashboard --port 8000
# http://127.0.0.1:8000 — 热度图 / 每 TCR 排名 / 金标准 / 交叉反应证据 /
# 可靠性闸门 / 域适应状态 / 提交进度 / KN-11 临床扫描
```

零依赖(标准库 `http.server`);`build_data(...)` 返回可 JSON 化的 dict。

---

## 5. CLI 入口一览

| 命令 | 作用 |
|---|---|
| `meta-acdc-dashboard` | 启动 TCR-Safety-Radar 看板 |
| `meta-acdc-score` | EGNN 打分(支持集成、温度、pLDDT 置零) |
| `meta-acdc-stability` | 跨实例稳定性判决门 |
| `meta-acdc-ensemble` | 多实例集成均值(排除退化 seed) |
| `meta-acdc-clinical` | 临床金标准结构扫描 |
| `meta-acdc-attribute` | 梯度归因可解释性 |
| `meta-acdc-mutate` | 逐位/累积突变扫描 |
| `meta-acdc-native-manifest` | 天然导入 + 候选路由 |
| `meta-acdc-select` | 蛋白组候选预筛与 EIG 排序 |
| `meta-acdc-train-egnn` | 晶体诱饵任务训练 |
| `meta-acdc-train-da` | 域适应训练 |

---

## 6. 典型配方

**从 AF3 结果到稳定排名**:

```bash
# 1) 导入(天然 + 候选路由)
meta-acdc-native-manifest --src folds_YYYY --out data/raw/af3_native_predictions \
    --manifest data/processed/af3_native/manifest.tsv \
    --candidate-out data/raw/af3_predictions
# 2) 域适应(5 seed)
meta-acdc-train-da --af3 data/raw/af3_native_predictions \
    --af3-manifest data/processed/af3_native/manifest.tsv --seeds 0 1 2 3 4
# 3) 打分 + 判决
for s in 0 1 2 3 4; do
  meta-acdc-score --cifs data/raw/af3_predictions \
    --model data/processed/da_seed${s}.model.pt --mask-plddt \
    --out data/processed/da_scores_seed${s}.tsv
done
meta-acdc-stability --scores data/processed/da_scores_seed{0,1,2,3,4}_ensemble.tsv
# 4) 集成 + 临床扫描
meta-acdc-ensemble --scores data/processed/da_scores_seed{0,1,2,3,4}_ensemble.tsv \
    --out data/processed/prediction_scores_ensemble.tsv
meta-acdc-clinical
```

或直接:`scripts/run_da_af3.sh <folds_dir>`(一键全流程)。
