# Meta-ACDC 项目科学性审查报告

**审查对象**：`meta-acdc-code.tar.gz`（代码）+ `meta-acdc-compute-data.zip`（数据与结果）
**审查日期**：2026-10-18
**审查方式**：源码逐行阅读 + 实际执行验证脚本 + 用 RCSB 权威序列核对结构链身份 + 文献原文核对
**总体判定**：**当前状态下，本项目不存在任何可引用的定量结论。**

---

## 0. 一句话结论

项目自带一份 28 条的外部审计（`docs/AUDIT.md`）和一份 README 顶部撤回声明，这个诚实性记录值得肯定。但我在审计之外发现了 **9 个更严重的问题**，其中 3 个属于"结果本身不可信"级别：

1. **交付的两张结果表对同一批结构给出互相矛盾的分数，且定性结论反转**（第 1 节）
2. **被解释的分数差（0.02–0.06）远小于分数自身的测量离散度（中位数 0.25）**（第 2 节）
3. **临床"金标准"案例存在受体身份替换：标签写 A3A，结构里是 MAG-IC3；且 4 行结果的结构已被自家 QC 判为不可信，却仍带 hazard 判定交付**（第 3 节）

同时，审计声称"已修"的 28 条中，**至少 8 条在交付的代码里并未真正修复**（第 6 节）。

---

## 1. 【致命】两张结果表自相矛盾，且结论反转

这是本次审查最重要的发现。

数据包里唯一交付的集成打分表是
`data/processed/prediction_scores_ensemble.relabeled.tsv`（21 行），
临床扫描结果是 `data/processed/kn11_clinical_scan.tsv`（7 行）。
**两者都自称是"5-seed 域适应集成"的输出，但对同一 job 给出完全不同的分数。**

| job_id | kn11 临床扫描 | 交付的集成表 | 差值 | 集成表 ±SD |
|---|---|---|---|---|
| `5brz_ESDPIVAQY` | **0.9483** | 0.7062 | +0.2421 | ±0.1472 |
| `5brz_EVDPIGHLY` | **0.8916** | 0.6828 | +0.2088 | ±0.2566 |
| `5brz_KVAKELVHFL` | 0.8417 | 0.6499 | +0.1918 | ±0.3065 |
| `5brz_ILAKFLHWL` | **0.6761** | **0.8095** | −0.1334 | ±0.2512 |

**后果：定性结论直接反转。**

- 在 kn11 表里，5brz 家族排名为 ESDPIVAQY(0.948) > EVDPIGHLY(0.892) > KVAKELVHFL(0.842) > ILAKFLHWL(0.676)。项目据此写下核心结论："2/3 致死脱靶与靶同水平，ILAKFLHWL 部分漏检"（README:378、README:345）。
- 在交付的集成表里，排名变成 **ILAKFLHWL(0.8095) > ELQHGLYAL(0.7920) > ESDPIVAQY(0.7062) > EVDPIGHLY(0.6828) > KVAKELVHFL(0.6499)**。那个被称为"漏检"的 ILAKFLHWL 变成了**家族第一，比靶肽高 0.127**；而被称为"反超靶肽、FATAL"的 ESDPIVAQY 只比靶肽高 0.023。

这不是数值抖动，是**同一个实验的两套结论**。README、paper1、paper2 全部引用的是 kn11 那一套（0.948/0.892/0.842/0.676），而这套数字在交付的数据里**没有任何来源表可以复现**。

**追溯失败**：`clinical_scan.py:38` 默认读 `data/processed/prediction_scores_ensemble.tsv`。我实测该文件、以及生成它所需的 `da_scores_seed{0..4}_ensemble.tsv`（5 个）、`egnn_dataset.pt`、`da_seed*.model.pt`（全部）**均不在数据包中**。`kn11_clinical_scan.tsv` 是一个孤儿文件——输入缺失，而唯一可比的表与它矛盾。

> 严重程度：**致命**。项目最响亮的结论（"模型不放过致死模拟肽"）建立在一张无法追溯、且与同包另一张表矛盾的数字上。

---

## 2. 【致命】被解释的效应量小于测量噪声本身

即使不追究上面哪张表是对的，两张表都记录了一个决定性事实：**单个 job 在 5 个 AF3 模型间的标准差，中位数是 0.2512。**

交付集成表 21 个 job 的 SD 分布：

- min 0.0105，中位数 **0.2512**，均值 0.2194，max 0.4260
- SD > 0.20 的：**13/21**
- SD > 0.30 的：**5/21**

把 5brz 家族里"与靶肽的差值"和"该 job 自己的 SD"放在一起看：

| job | \|delta vs target\| | 自身 SD | 可分辨？ |
|---|---|---|---|
| 5brz_ILAKFLHWL | 0.1267 | 0.2512 | **否，delta < SD** |
| 5brz_ELQHGLYAL | 0.1092 | 0.2646 | **否** |
| 5brz_KVAKELVHFL | 0.0329 | 0.3065 | **否** |
| 5brz_ESDPIVAQY | 0.0234 | 0.1472 | **否** |

**全部 4 个 delta 都小于对应 job 自身的标准差。**

而 README:346 和 paper1 恰恰是把这些 0.02–0.06 的差值当作科学信号来解读的：

> "两个**致死**脱靶（MAGE-A12 模拟肽、titin）分别只低 0.02 / 0.06——模型无法把它们与靶肽分开"

这句话在方法上是站不住的。**在一个自身重复离散度就有 ±0.15~0.31 的测量上，0.02 的差值不携带任何信息**——它既不能证明"无法区分"，也不能证明"与靶同水平"。要主张"分数接近"，必须先给出分数分辨率的下界，而这个下界比被比较的差值大一个数量级。

注意这个 SD 还**只是** AF3 结构采样分量。种子间分量（就是那个 mean pairwise r = 0.644 的稳定性分析所度量的东西）是**额外的**，未包含在内。真实离散度更大。

> 严重程度：**致命**。临床安全判定的全部推理，建立在一个小于测量噪声的差值上。

---

## 3. 【致命】临床"金标准"案例的受体身份替换 + QC 闸门被绕过

### 3.1 打分用的结构确实是忠实的 5BRZ，但它不是临床致死的那个受体

我用 RCSB 官方 FASTA（`https://www.rcsb.org/fasta/entry/5BRZ`）逐字比对了
`data/raw/af3_predictions/5brz_EVDPIGHLY_model_0.cif` 的五条链：

| 项目 CIF | 角色 | 长度 | 对应 RCSB 5BRZ | 比对 |
|---|---|---|---|---|
| chain A | tcr | 197 | Chain D（TRAV21, TCR α） | **逐字相同** |
| chain B | tcr | 241 | Chain E（TRBV5-1, TCR β） | **逐字相同** |
| chain C | peptide | 9 | Chain C（EVDPIGHLY） | **逐字相同** |
| chain D | mhc | 275 | Chain A（HLA-A1 α 链） | **逐字相同** |
| chain E | b2m | 100 | Chain B（β2-微球蛋白） | **逐字相同** |

结构本身没问题，是 5BRZ 的忠实重预测。**但受体身份有问题。**

Raman 2016 原文（Sci Rep 6:18851，PMID 26758806，我已核对全文）明确区分了两个受体：

- **a3a**：临床给药、导致 2 名患者心脏死亡的亲和力增强 TCR
- **MAG-IC3**：5BRZ/5BS0 晶体结构里的受体。原文措辞是"a soluble higher affinity variant of the same HLA-A\*01-restricted MAGE-A3 wild type TCR from which a3a was derived was engineered for use as an **ImmTAC reagent**"

即 **5BRZ 里结晶的是 ImmTAC 试剂 MAG-IC3，不是致死的临床受体 a3a**。二者同源但不同。

而 `kn11_clinical_scan.tsv` 的 `tcr` 列写的是：

```
A3A affinity-enhanced MAGE-A3 TCR
```

`clinical_gold_standard.tsv` 第一行写的是 `A3A (affinity-enhanced MAGE-A3 TCR, HLA-A*01)`。

**项目把 MAG-IC3 的结构打分，标注为 A3A 受体的结果，并以此声称验证了 A3A 的致死交叉反应。** 这是受体替换，AUDIT.md #6 修了 PMID 和 5BS0 归属，但没有触及这一点。

### 3.2 KVAKELVHFL 是跨案例拼接

`clinical_gold_standard.tsv` 自己第 2 行写明：

```
A2-restricted MAGE-A3 TCR (A118T affinity-enhanced) | KVAELVHFL | ... | PMID:23470321
```

即 KVAELVHFL/KVAKELVHFL 属于 **A\*02:01 限制性**的另一个案例（Morgan 2013，MAGE-A12 神经毒性）。

但 `kn11_clinical_scan.tsv` 把 KVAKELVHFL 挂在 `A3A-MAGEA3 / pdb=5brz / restriction=HLA-A*01` 下，当作 A3A 案例的 off-target 打分。

**一个 A\*02:01 案例的肽，被放到 A\*01 限制性受体（MAG-IC3 + HLA-A1）的结构上打分，然后作为 A3A 致死案例的脱靶结果报告。** 项目自己的两个数据文件互相矛盾。

### 3.3 QC 闸门：声称强制，交付数据仍是闸前口径

`af3_qc.tsv` 共 261 行，标记分布：`''`=147、`no_confidence`=96、`low_iptm`=18。
我按 `score_predictions.py` 的 (instance, job) 联合键口径统计 188 个候选 CIF：

**pass = 93，flagged = 95** —— 与 AUDIT.md #26 记录的"93 通过 / 95 带标记被闸"完全一致。

现在看 `kn11_clinical_scan.tsv` 的 7 行：

| case | peptide | score | QC 状态 |
|---|---|---|---|
| A3A-MAGEA3 | EVDPIGHLY | 0.8916 | pass |
| A3A-MAGEA3 | ESDPIVAQY | 0.9483 | pass |
| A3A-MAGEA3 | ILAKFLHWL | **0.6761** | **low_iptm (ipTM=0.540)** |
| A3A-MAGEA3 | KVAKELVHFL | 0.8417 | pass |
| A3A-cross-structure | EVDPIGHLY | **0.8829** | **no_confidence** |
| KIMMTRAK | YLEPGPVTA | **0.5476** | **no_confidence** |
| KIMMTRAK | YLEPGPVTL | **0.9106** | **no_confidence** |

**7 行里 4 行的结构已被项目自己的 QC 判为不可信，却仍带着具体分数留在交付数据中。**

而 `score_predictions.py:99-139` 的实现是：QC 缺失或有 flag 时**写空 score、不给数字**。所以 kn11 这张表**不可能**是闸门生效后的产物——它是闸前口径。

最严重的是 ESDPIVAQY 那行同时带着 `hazard_flag = "lethal-mimic scores at/above the intended target"`。**一个已被自家 QC 部分否定的数据集，仍在输出一条"致命模拟肽达到或超过靶点"的危险判定。**

AUDIT.md 第三轮自己承认"0.644/188 一律标注闸前口径"，但 kn11 这张**数据文件**没有任何此类标注，且被 README/paper1/paper2 直接引用。

> 严重程度：**致命**。三重问题叠加：受体标错、跨案例拼接、QC 闸门未落实。

---

## 4. 【严重】审计声称"已修"的问题，在交付代码里并未修复

我逐条实测了 AUDIT.md 的修复声明。**至少 8 条与代码实际状态不符。**

### 4.1 #2 残基编码：graph.py 修了，诱饵生成没修

AUDIT.md:250 声称"改走 AA3TO1"。`graph.py:266` 确实改了。但同一条坏映射在**诱饵构造路径**里原封不动：

```
src/meta_acdc/models/dataset.py:260           native_seq = "".join(r.resname[:1] ...)
src/meta_acdc/models/dataset.py:270           if "".join(r.resname[:1] for r in candidate) != native_seq
src/meta_acdc/models/train_domain_adapt.py:132 native_seq = "".join(r.resname[:1] ...)
src/meta_acdc/models/train_domain_adapt.py:145 if "".join(r.resname[:1] for r in cand) != native_seq
```

这是**嫁接诱饵供体去重的守卫条件**——用来避免"供体肽 == 天然肽"导致诱饵等于正样本。它仍在把 ARG/ASN/ASP→A、GLU→G、LYS→L、PHE→P、TRP/TYR→T 折叠。

我实测了项目自己面板里的真实肽对：

```
native AAGIGILTV   donor NAGIGILTV
  坏映射 → AAGIGILTV vs AAGIGILTV  判定"相同" → 跳过诱饵
  真序列 → AAGIGILTV vs NAGIGILTV  实际不同
```

**后果：诱饵集的构成是错的，且没有任何测试覆盖这条路径。** `AAGIGILTV`/`NAGIGILTV`、`LLFGYPVYV`/`LLFGKPVYV` 这类只差一个残基的肽对项目大量存在（面板 91 个唯一肽里同家族变体极多），去重判断在这些对上系统性失效。

### 4.2 #4 训练/验证切分：train_egnn 修了，产出 0.644 的 DA 路径没修

AUDIT.md:252 声称"组键归一化到 PDB/复合物 + 受体家族前缀"，`train_egnn.py:145-153` 确实实现了。

**但报告 0.644、0.609、临床扫描结果的模型是 `train_domain_adapt.py` 训练的，它用的是自己那份未修的 `pdb_group()`（L194-195）：**

```python
def pdb_group(name: str) -> str:
    return name.split("-disp")[0].split("<-")[0]
```

我用项目自己的代码实测同一个生物学复合物的各种命名：

```
5brz                     -> group '5brz'
5brz_0                   -> group '5brz_0'
5brz_1                   -> group '5brz_1'
5brz-disp                -> group '5brz'
5brz<-3uts               -> group '5brz'
af3-5brz_EVDPIGHLY       -> group 'af3-5brz_EVDPIGHLY'
af3-5brz_EVDPIGHLY-disp  -> group 'af3-5brz_EVDPIGHLY'
tcr-5brz_EVDPIGHLY       -> group 'tcr-5brz_EVDPIGHLY'

=> 一个复合物被切成 5 个"独立"组
```

这正是 AUDIT #4 描述的缺陷本身：**晶体版可以进训练集、它的 AF3 重预测版进验证集**。而域适应的全部前提就是"把 AF3 预测的天然复合物混入训练"，所以这条泄漏路径在 DA 模型里**必然被触发**。

进一步确认修复未接线：

```
train_domain_adapt.py 导入 family_map？   False
train_domain_adapt.py 有 --group-by？     False
train_domain_adapt.py 有 --family-map？   False
train_domain_adapt.py 调用 train_egnn.pdb_group？ False

run_da_af3.sh            : 用 train_domain_adapt.py=True, 传 --family-map=False, --group-by=False
da_flow_gpu_20261012.sh  : 同上
da_retrain_collatefix.sh : 同上
```

**修好的分组逻辑存在于一个不产出报告数字的脚本里；产出报告数字的脚本用的是坏的那个。**

> 这意味着 `val AUROC 0.947–0.959` 和 0.644 都在存在验证集污染的条件下取得，"域适应提升"不能与"泄漏提升"区分开。

### 4.3 #20 固定预算消融：声称已修并验证，实际该路径根本无法运行

AUDIT.md:268 声称"末批补齐到同一 budget；各 batch size 花费相同"。

实测：`simulate_al.py:290` 是

```python
rounds_b = max(1, math.ceil(budget / b))
```

而该文件的 import 列表是 `argparse, csv, sys, numpy`（+ 3 个项目内模块）——**从未 `import math`**。

pyflakes 全仓扫描，唯一的 undefined name 就是这一处：

```
src\meta_acdc\active_learning\simulate_al.py:290:31 undefined name 'math'
```

直接执行同一语句确认：`NameError: name 'math' is not defined`。

**`python simulate_al.py --ablation` 在打印任何东西之前就崩溃。** README:244 报告的"大批(3600×2)24.5% > 小批(300×20)21.6%"无法从此代码复现，AUDIT 声称的验证也不可能是对这条路径做的。

### 4.4 #1 / #7 / #13 / #22：归属的模块在两个归档中都不存在

AUDIT 把这四条的修复分别归给
`data/ingest_vdjdb.py`、`data/ingest_iedb.py`、`data/splits.py`。

实测两个归档：

```
code tarball  : src/meta_acdc=True   src/meta_acdc/data=False
data zip      : src/meta_acdc=False  src/meta_acdc/data=False
```

**整个 `meta_acdc.data` 包不存在。** 三个模块因此在导入期就失败：

```
可导入: 48 / 51
FAIL meta_acdc.active_learning.simulate_al   ModuleNotFoundError: No module named 'meta_acdc.data'
FAIL meta_acdc.models.meta                   ModuleNotFoundError: No module named 'torch.func'
FAIL meta_acdc.models.train_maml             ModuleNotFoundError: No module named 'torch.func'
```

依赖 `meta_acdc.data.splits` 的还有 `select_candidates.py:209`、`benchmark_baseline.py:148`。

**后果**：这 4 条修复（VDJdb 标签语义、TRA/TRB 分链、切分聚类、IEDB 终点语义）**无法被审阅、无法被运行、无法被验证**。README:502 给出的复现命令
`python src/meta_acdc/data/download_vdjdb.py && python src/meta_acdc/data/ingest_vdjdb.py`
指向不存在的文件。同时 `train_maml.py`/`meta.py` 依赖 `torch.func`（PyTorch 2.0+），与项目声明的 torch 2.9.1 兼容但与当前可运行的最低环境不符——MAML++ 那条 0.625 的元学习结果同样无法复现。

### 4.5 #6 临床 PMID：代码改了，交付的数据文件没改

AUDIT.md:254 声称"PMID 改为已核实的 23926201"。`kn11_clinical_scan.tsv` 里确实是 `PMID:23926201`。

但 `clinical_gold_standard.tsv`（同一目录、同为交付数据）仍是：

```
PMID:23999400; PMID:23999405
```

AUDIT #6 原本就指出 **PMID:23999405 对应的是稀土材料电子结构研究**，Cameron titin 原文是 23926201。**这个已知错误的 PMID 仍留在交付数据里，且与同目录另一张表的引用不一致。**

### 4.6 #26 QC 闸门：代码实现了，交付数据是闸前的

见第 3.3 节。`af3_qc.tsv` 的 93/95 划分与 AUDIT 记录一致，说明闸门代码是对的；但 `kn11_clinical_scan.tsv` 里 4 行带 flag 却有分数，说明**交付的结果表不是闸门生效后的产物**。

### 4.7 #9 bootstrap 分层：改成按 PDB 聚类，仍不是按受体聚类

`da_stability.py:88-89` 的聚类键是 `job.split("_", 1)[0]`，即 **PDB id**。

我对实际打分的 188 个候选 CIF 解析 TCR 链序列后统计：

```
候选数                     : 188
bootstrap 聚类单位 (PDB id) : 72      <- da_stability.py 重采样的单位
不同 TCR 链序列            : 59      <- 真实生物学独立单位
有 >1 个候选的受体          : 28
单个受体最多候选数          : 12
```

**72 vs 59**：同一个受体序列被多个 PDB id 携带（晶体 + 重预测 + 拆分复合物），被当作独立聚类重采样。注释里写的是"resample RECEPTOR groups"，实现的是 PDB 组。**报告的 CI [0.572, 0.706] 偏窄。**

### 4.8 测试入口路径不存在

AUDIT.md:230 给出的验证入口是 `python -m pytest src/meta_acdc/tests -q`。
实际测试在仓库根的 `tests/`，`src/meta_acdc/tests` 不存在。

---

## 5. 【严重】诱饵任务的三个实现缺陷

0.933 这个"侧链接触突破"是全文的核心论点，它完全建立在嫁接诱饵的有效性上。我直接调用项目的诱饵构造函数实测，发现三个未被审计记录的问题。

### 5.1 长度不匹配时按索引截断，没有锚点重排

`dataset.py:119` 用 `zip(native_pep, peptide)` 配对，然后 `swap = {id(r): f ...}`。

我实测 9 肽受体 + 13 肽供体：

```
grafted residues: [(1,'PRO'),(2,'LYS'),(3,'TYR'),(4,'VAL'),(5,'LYS'),
                   (6,'GLN'),(7,'ASN'),(8,'THR'),(9,'LEU')]
=> 13 个供体位置只用了 9 个；对齐方式是索引，不是 MHC 锚点位置
```

面板里的肽长 8–13 不等（`RLPAKAPL`=8、`KAPYDYAPI`=9、`RLPAKAPLLGCG`=12、`PKYVKQNTLKLAT`=13），**长度不匹配的嫁接是常规情形而非例外**。按索引对齐意味着 P2/P9 锚点残基被错位移植——供体的锚点落到受体的非锚点位置上，产生的界面在物理上无意义，与"保持骨架、只换化学"的设计意图不符。

### 5.2 供体骨架不完整时静默回退，但残基名仍被改写

`dataset.py:79-97` `_transplanted_atoms`：

```python
src = _backbone_coords(foreign); dst = _backbone_coords(native)
if src is None or dst is None:
    return list(native.atoms)      # 静默回退：保留天然原子
```

而调用方 `make_decoy`（L125）**无条件写入供体的 resname**：

```python
out.append(Residue(r.chain, f.resname, r.resid, ...atoms))
```

我实测（供体 ARG 缺 C 原子，低分辨率结构中常见）：

```
graft result : ['C','CA','CB','CD1','CD2','CG','N','O']   <- 全是天然 LEU 的原子
通过 make_decoy() 后: chain P: resname=LEU(写入的是供体名) atoms=['C','CA','CB','CG','N','O']
```

**结果是一个自相矛盾的残基：节点身份 one-hot 说它是供体氨基酸，几何与侧链接触直方图说的是天然氨基酸。** 而侧链接触直方图正是 0.933 所依赖的唯一决定性特征（消融：去掉它 → 0.501）。这类样本既不是有效正样本也不是有效负样本，且**回退是静默的，不计数、不告警、不记录**。

### 5.3 位移诱饵沿固定全局 +z 轴平移

`dataset.py:146-148`：

```python
shifted_atoms = [(n, x, y, z + shift) for n, x, y, z in r.atoms]
```

沿**沉积文件坐标系的全局 +z** 平移 6 Å，没有任何检查确认 +z 对该条目指向 MHC 沟槽外侧。对某些复合物，肽会侧移甚至移向 TCR。

而 `build_interface_graph_from_residues:360-362` 有一个静默兜底：

```python
if len(keep) < 5:
    keep = set(range(len(residues)))   # 界面太小 → 退回全复合物
```

docstring 声称位移诱饵"STAYS WITHIN the 10A interface selection"，这个保证是**取向依赖的**，且失败会被上面的兜底吞掉。位移诱饵 AUROC 0.946–0.980 因此无法排除部分样本是"肽移入 TCR 界面"或"节点集退化为全复合物"的伪影。

### 5.4 根本性问题：诱饵从未被证明不结合

AUDIT.md:287 自己承认："嫁接诱饵没有实验依据证明不结合；在缺少真实负样本前，0.933 只能读作'区分天然结构与按该规则构造的诱饵'"。

这个自我限制是对的，但**它实际上否定了 0.933 与本研究目标之间的全部联系**。研究目标是预测真实脱靶风险；任务实际是识别"哪些结构是脚本按规则造出来的"。这两者之间没有任何已验证的桥梁。结合 5.1–5.3 的实现缺陷，模型很可能学的是**构造规则的可检测痕迹**（锚点错位、身份/几何不一致、固定轴向位移），而非分子相容性。

---

## 6. 【严重】面板冗余：名义样本量远大于有效样本量

### 6.1 待上传 FASTA 面板（347 个文件）

```
总 FASTA 文件              : 347
唯一 PDB 前缀              : 111
唯一肽序列                 : 91
精确唯一 (α,β) 受体对       : 98
98% 一致性下的受体簇        : 73      <- 真实生物学独立单位
唯一 MHC 序列              : 45
唯一 B2M 序列              : 11
```

最大的冗余簇：`3uts / 5c07 / 5c09 / 5hyj` 四个 PDB 共享同一受体，共 32 个文件；`3h9s / 3pwp / 3qfj` 三个 PDB 24 个文件；`5c0a / 5c0b` 16 个文件；`3d39 / 3d3v` 16 个文件。

我逐字比对了同族文件的链内容，例如 `1qsf_LLFGFPVYV` / `2gj6_LLFGYPVYV` / `3d39_LLFGYPVYV` / `3h9s_LLFGYPVYV`：**α 链、MHC、B2M 完全相同，β 链仅个别残基差异，肽相同**。这些是同一受体的近重复结构，不是独立观测。

### 6.2 已打分面板（188 个 CIF）

```
候选数                     : 188
bootstrap 聚类单位 (PDB id) : 72
不同 TCR 链序列            : 59
单受体最多候选              : 12
唯一 MHC 重链序列          : 25（最大的两个各占 52、46 个）
```

**188 个候选实际上只覆盖 59 个受体和 25 个 MHC 重链，且约 52% 的候选集中在 2 个 MHC 序列上。**

由此产生的具体问题：

1. **稳定性 CI 偏窄**（见 4.7）：按 72 个 PDB 组重采样，而非 59 个受体。
2. **组内排名不可报告**：README:320 自己记录了 5/17 组的组内 r < 0.5，但结论仍写"全面板级排名不受影响"。当面板本身由近重复结构构成时，"全面板级 Spearman"度量的很大程度是**模型能否对同一受体的重复结构给出一致分数**，而不是能否在生物学上区分候选。
3. **面板规模口径从未冻结**：README:137 说"172 天然 + 154 候选"，README:302 说"188 候选"，README:317 说"154 候选面板"，`af3_pending_upload` 有 347 个待上传 FASTA，实测 `af3_predictions` 有 188 个 model_0 CIF。四个数字并存，无法确定任一报告结果对应哪个集合。

### 6.3 序列伪影被当作天然序列提交给 AF3

```
beta 链以 "NA" 开头 : 53 / 347 个文件
MHC 链以 "GS" 开头  : 279 / 347 个文件
MHC 链以 "MG" 开头  : 16 个文件
```

β 链 N 端变体统计：`NAGVTQTPKFQV`(53)、`GVTQTPKFQVLK`(40)、`DAGVIQSPRHEV`(37)、`AVTQSPRSKVAV`(30)…

这些是从 PDB SEQRES 带出的 **N-糖基化位点占位（ASN→NA）与表达标签（GS/MG）残留**，不是天然序列。它们被原样写进 FASTA 提交给 AF3，导致：

- 同一受体的不同 PDB 条目产生**长度不同的伪变体**（这正是 6.1 里"98 个精确受体对 vs 73 个 98% 簇"差异的来源之一），人为膨胀了名义独立性；
- `classify_chains`（`graph.py:126-137`）用**长度窗口**判角色（TCR 靠 YFC/YLC 基序、MHC >250、B2M 90–105、肽 5–20），标签残留直接改变落入哪个窗口。实测面板链长：α 112–243、β 194–250、MHC 212–277、B2M 96–101。**MHC 下界 212 已低于 `>250` 的判据**，靠的是 `150<=n<=250 → MHC(class II)` 这条兜底——也就是说面板里的 I 类 MHC 有大量条目是被**当作 II 类**分对的。这条分类逻辑没有测试覆盖，一旦某条链落到 105–150 或 >250 的间隙就会被判成 `OTHER`，然后静默并入 pMHC 侧。

---

## 7. 【严重】图表不是派生产物

`scripts/make_figures.py`（745 行，无 argparse 参数）逐函数检查数据来源：

| 图函数 | 是否读数据文件 | 硬编码的三位小数 |
|---|---|---|
| `fig1_ablation` | **否** | 0.933, 0.924, 0.899, 0.606, 0.501, 0.500, 0.495 |
| `fig2_heatmap` | 是 | — |
| `fig3_resubmission` | 是 | — |
| `fig5_instability` | 是 | — |
| `fig6_domain_adapt` | 是 | **0.644** |
| `fig7_crossreactivity` | 是 | — |
| `figp2_recall` | **否** | — |
| `figp2_prefilter` | **否** | — |
| `fig8_cross_method` | 是 | 0.008 |

`fig1_ablation`（L101-111）把整条消融序列写成字面量：

```python
configs = [
    ("Full (v9.1)", 0.933, True),
    ("− side-chain extent", 0.924, False),
    ("− residue identity", 0.899, False),
    ("− physicochemical", 0.606, False),
    ("− contact histogram", 0.501, False),
    ("Cα-only (v5-v8)", 0.500, False),
    ("Sequence baseline", 0.495, False),
]
```

我不认为这是数据造假——这些数字应当来自训练日志。但问题是**方法学上的**：

1. **图与数据完全解耦**。重跑实验不会改变图，`docs/figures/` 里已提交的 PNG/SVG 无法与任何原始输出对账。
2. **结合第 4.4 节**（模型 checkpoint、`egnn_dataset.pt`、per-seed 打分表全部缺失），**没有任何途径验证这些字面量曾经由某次运行产生过**。
3. `fig6_domain_adapt` 一边读数据一边硬编码 0.644，意味着即使数据变了，图上那个决定性数字也不会变。

一个可复现研究的图，应该是从结果表机械生成的。

---

## 8. 【严重】缺失的对照与基准

### 8.1 无任何同数据集的已发表工具对比

`docs/benchmarks.md:1240-1246` 记录：ERGO-II、NetTCR-2.2、MCGLPPI++、ZeroBind **均因环境兼容性问题未运行**，仅引用文献数字作"区间参照"。

于是论文的核心论证链——"序列方法卡在 0.62–0.67（TITAN 0.62、ERGO 0.669），所以必须转向结构"——是**跨数据集、跨切分、跨评估协议的数字并列**。这些工具在自己的基准上的表现，不能推出它们在本项目数据上的表现。"序列天花板"这个前提从未在本项目数据上被任何外部工具验证过。

### 8.2 无随机肽基线

模型对**完全不相关的随机肽**打多少分？整个项目里没有这个基线。

这是解释第 1、2 节那些数字的必要参照。如果模型对任意肽都给出 0.6–0.8，那么 0.676（ILAKFLHWL）的"部分漏检"和 0.948（ESDPIVAQY）的"FATAL 反超"都只是同一分布内的波动。没有这个基线，**任何绝对分数都无法解释**。

### 8.3 MHC 预筛基准样本量不足且未报不确定性

50 个已知结合肽 + 50 个组成打乱诱饵，得出"top-2% 保留 56%、压缩 ~50×"。

n=50 的召回率 95% CI 约为 **[42%, 70%]**（±7% 标准误），论文未报告。更关键的是 AUDIT #10 已指出**"top-2% 阈值外推为蛋白组压缩 50 倍"需要收缩**——从 100 个肽的组成打乱对照，外推到 465 万肽的真实蛋白组，两者难度不在一个量级。而且"组成打乱"是最容易的负样本类型（破坏锚点基序），远低于真实筛选场景的难度。

### 8.4 计划书承诺的 80–90% 灵敏度无任何支撑

`committee_proposal.md` 声称主动学习"expected to achieve 80–90% Sensitivity in identifying known-category high-risk off-targets"。项目中最接近的结果是 AL 模拟的 2.0× 召回增益——**该结果已因 VDJdb score 被当作结合标签而全面撤回**（AUDIT #1），且其消融路径本身无法运行（第 4.3 节），且在严格预筛下增益消失（README:246 记录 0.9 阈值下为 1.0×）。

### 8.5 湿实验、文库合成、动物实验全部未做

计划书 Module 2（寡核苷酸池合成→慢病毒→FACS→NGS 闭环）、Module 3（原代 T 细胞杀伤 + NSG 小鼠）、Deliverable 5（百万级 ACDC 验证数据集）、Deliverable 6（实物高危抗原文库）——`docs/milestones/KEY-NODES.md` 明确标注"不执行（技术路线保留）"。

这本身是合理的范围界定，README:24 也说明了。**但 paper2 的标题仍是 "An AI-driven 'Dry-Wet' Loop"，摘要仍以"digital twin of the dry-wet loop"立论**，正文未清楚交代没有任何湿实验数据。这会误导审稿人对证据强度的判断。

---

## 9. 【中等】循环论证与选择性报告

### 9.1 自设门限当作客观标准

"多种子 Spearman ≥ 0.5 才可报告排名"这个门限是项目自己定的。AUDIT #9 已承认"0.5 门槛是自设的，没有真实结合标签验证其意义"，代码也加了正确的免责声明（`da_stability.py:110-122`）。

**但叙事框架没有改**：README:301-302 的表格用"判决"列写"稳定""判决不变"，README:377 写"Rashomon 危机解除"，`run_da_af3.sh` 末行输出"STABLE iff mean pairwise Spearman r >= 0.5"。读者接收到的信息是"通过了客观检验所以结果可信"，而实际是"通过了一个无外部效度校准的自设阈值"。

### 9.2 CI 下界压线仍判 STABLE

`benchmarks.md:641-643` 记录 144-native 时 bootstrap 95% CI 为 **[0.495, 0.854]**，下界 0.495 **低于** 0.5 门限，且 n=21 候选。按项目自己的规则应判 provisional，文中写"STABLE(初步)"并在后续简称为 STABLE。

### 9.3 best-of-N 种子选择性报告

paper1 L165-167 承认嫁接任务 3 个种子给出 0.920 / 0.741 / 0.692，**跨度 0.23 AUROC**，然后写"all reported numbers use the best of 3 seeds, and we recommend multi-seed reporting"。

推荐多种子报告，但自己全文所有 headline 数字都是单种子最优值。0.933 这个被反复引用的数字，是三个种子中最好的那个，而最差的是 0.692——**接近 0.606（去理化标量消融）**。

### 9.4 种子质量过滤规则是事后制定的

paper1 §2.6 L278-279："train several seeds and report a *seed-filtered* ensemble (non-degenerate, val-AUROC-qualified)"。

"哪个种子塌缩"是**观察之后**才知道的。没有预注册的纳入标准，这是事后选择。README:312 记录"5 个种子里有 1 个（val AUROC 0.804）会塌缩成常数输出"——注意那个被剔除的种子 val AUROC 是 0.804，**并不低**。剔除标准到底是什么，未见明确定义。

### 9.5 "单调改善"的叙事跳过了退化点

paper1 L462 称不稳定性"vanishes monotonically (0.447 → 0.479 → 0.729)"。但 README:298 的表格自己记录了 **63 个天然正样本时 Spearman = NaN（种子塌缩）**。序列是 0.447 → **NaN** → 0.479 → 0.729，中间有一个未定义点。AUDIT #10 已把"扩样单调解决不稳定"列入需撤回的表述，但 paper1 Discussion 仍在用。

README:305-308 后来也承认"改善并非单调，正样本风格（预测风 vs 晶体风）比数量更关键"——这与 paper1 L462 直接冲突，两处并存。

### 9.6 Fisher 合并把非独立结构当独立

`rank_stats.py:82-90` 的代码注释自己承认："this treats the listed structures as independent; receptors from the same family are not, so treat this as an upper bound on significance until a family-stratified version is run."

这个自我限制是正确的，但**"family-stratified version"从未实现**。`benchmarks.md` 报告的 pooled Fisher P=0.66（作者诚实标注 underpowered）用的仍是这个非独立版本。

---

## 10. 【中等】可复现性与工程完整性

### 10.1 复现所需的关键产物全部缺失

| 需要的文件 | 状态 |
|---|---|
| `da_scores_seed{0..4}_ensemble.tsv` | **全部缺失** |
| `prediction_scores_ensemble.tsv`（clinical_scan 的输入） | **缺失** |
| `egnn_dataset.pt` | **缺失** |
| `da_seed*.model.pt` | **全部缺失** |
| `af3_qc.tsv` | 存在（261 行） |
| `prediction_scores_ensemble.relabeled.tsv` | 存在（21 行） |
| `colabfold_scores_ensemble.tsv` | 存在（172 行） |
| `prediction_scores.relabeled.tsv` | 存在（105 行） |

**0.644、CI [0.572, 0.706]、val AUROC 0.947–0.959、全部临床分数，都无法从交付物复现。** 只能靠重跑 AF3 + 重训（而 AF3 网页提交本身不保证可重复）。

### 10.2 版本控制审计链断裂

数据仓 git 历史只有 **1 个 commit**：

```
a3da597 Recover published 1D3 engineered DNA and quarantine source-checked SPR candidates
```

AUDIT.md 里引用了 `98773bd`、`f59dd45`、`a3b8715`、`21926a9`、`60fb071`、`58df6ab` 等提交号。**这些提交在交付的仓库中都不存在。** 无法确认哪个数字产生于哪次运行、哪次修复之前还是之后。这对一个以"逐条审计修复"为核心叙事的项目是根本性缺陷。

### 10.3 测试覆盖的盲区恰好是所有关键路径

实测：`23 tests OK, 6 skipped`。跳过的 6 个正是：

```
test_cif_parse_and_identity        fixture missing: af3_instances/A/1ao7_LLFGYPVYV_model_0.cif
test_qc_table                      fixture missing: af3_qc.tsv
test_relabel_map_present_and_sane  fixture missing: af3_relabel_map.tsv
test_features_appended_and_respond structure fixtures missing
test_attribution_runs              attribution fixtures missing
test_entry_points_importable       tomllib unavailable
```

**所有依赖真实结构数据的路径都没有测试。** 而被测的 17 个都是纯函数级测试。这解释了为什么第 4 节那些"已修"的缺陷能长期存在：诱饵去重（4.1）、DA 切分（4.2）、消融路径（4.3）**没有任何一条被测试触及**。

另外，README:519 声称"tests/ 回归测试（16 项）"，实际是 23 项；AUDIT.md:230 给的 pytest 路径 `src/meta_acdc/tests` 不存在（实际在根 `tests/`）。

### 10.4 `kn5_submission_list.tsv` 的列名与内容仍相反

README:534 记录过"提交清单的 `tcr_a/tcr_b` 两列写反 → 曾把 48/73 条结构标签张冠李戴；改用链对匹配 + 晶体链指纹回贴修复"。

实测该清单 5brz 行：

```
tcr_a 列 == RCSB β 链（AGVTQTPRYLIKTRGQ...）？ True
tcr_b 列 == RCSB α 链（AQEVTQIPAALSVPEG...）？ True
```

**列标签与内容仍然相反。** 修复是在下游做了链对匹配（实际打分的 CIF 链身份是正确的，见 3.1 表），但**源清单的错误没修**。任何直接读这张清单的下游代码都会拿到反的 α/β。

### 10.5 `kn5_submission_list.tsv` 有重复行

同一 (pdb, tcr_a, tcr_b, peptide) 组合出现多次（5brz/5bs0/5eu6 各出现 2 次以上，2252 行）。README:556 记录过"job_id 判重系统性漏判……batch1 曾因此废掉 24/60 个提交"。清单里的重复未被清理。

---

## 11. 需要更正的我自己的错误结论

审查过程中我一度得出过一个**错误**判断，必须明确撤回：

> **撤回**：我最初用一组凭记忆书写的 HLA 参考序列做比对，得出"5brz 的 MHC 链是 A\*02:01（identity 0.891），而临床扫描声称 HLA-A\*01，构成等位基因替换"的结论。

**这个结论是错的。** 用 RCSB 官方 FASTA（`https://www.rcsb.org/fasta/entry/5BRZ`）逐字比对后确认：项目 CIF 的 MHC 链与 5BRZ Chain A（RCSB 标注为"HLA class I histocompatibility antigen, **A-1** alpha chain"）**字符完全相同**，五条链全部逐字匹配。

错误原因：HLA-A 各等位基因之间本就有 89–95% 同一性，我记忆中那份"HLA-A\*01:01"参考序列与真实的 5BRZ A1 链不符（真实链含 `...FVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAH...` 这类 A1 特有片段），导致一个大致正确的 A\*02:01 参考反而胜出。**教训正是这个项目自己也踩过的坑：不能用记忆里的参考序列做身份判定。**

需要区分的是：
- ❌ **不成立**：MHC 等位基因被替换成 A\*02:01。结构忠实于 5BRZ/5BS0（HLA-A1）。
- ✅ **成立**：结构里的受体是 **MAG-IC3**，而 `kn11_clinical_scan.tsv` 的 `tcr` 列标注为 **A3A**（临床致死的那个受体）。受体替换问题依然成立，只是层次不同（不是 MHC，是 TCR）。

同样，我核实了 `0.285 Å` 这个引用：Raman 2016 原文确有"the conformation of the Cα backbone and the orientation of the peptide side chains were very similar (root mean square deviation, or rmsd, based on Cα atoms at **0.285 Å**)"。**这个引用准确**，我此前子代理将其标为存疑是不必要的。同文还有"an rmsd between the two complexes of **0.418 Å**, based on Cα atoms in the whole model"，项目未引用。

---

## 12. 问题清单汇总

| # | 问题 | 位置 | 级别 |
|---|---|---|---|
| 1 | 两张结果表对同一 job 分数矛盾，定性结论反转 | `kn11_clinical_scan.tsv` vs `prediction_scores_ensemble.relabeled.tsv` | 致命 |
| 2 | 被解释的 delta(0.02–0.06) < 分数自身 SD(中位 0.2512) | 同上，21 个 job 中 13 个 SD>0.20 | 致命 |
| 3 | 受体身份替换：标 A3A，结构是 MAG-IC3 | `kn11_clinical_scan.tsv` tcr 列 vs 5BRZ | 致命 |
| 4 | KVAKELVHFL 跨案例拼接（A\*02:01 案例肽 → A\*01 结构） | `clinical_gold_standard.tsv` 行2 vs `kn11` | 致命 |
| 5 | QC 闸门未落实：4/7 临床行带 flag 却有分数，含 hazard 判定 | `kn11_clinical_scan.tsv` + `af3_qc.tsv` | 致命 |
| 6 | AUDIT #2 未真修：诱饵去重仍用 `resname[:1]` | `dataset.py:260,270`; `train_domain_adapt.py:132,145` | 严重 |
| 7 | AUDIT #4 未真修：DA 路径 `pdb_group()` 把 1 复合物切 5 组 | `train_domain_adapt.py:194-195` | 严重 |
| 8 | AUDIT #20 未真修：`--ablation` 路径 NameError（缺 `import math`） | `simulate_al.py:290` | 严重 |
| 9 | AUDIT #1/#7/#13/#22 归属模块整体缺失，3 个模块导入失败 | `meta_acdc.data` 包不存在 | 严重 |
| 10 | AUDIT #6 未落实到数据：错误 PMID 23999405 仍在交付表 | `clinical_gold_standard.tsv` | 严重 |
| 11 | AUDIT #9 未真修：bootstrap 按 PDB(72) 聚类而非受体(59) | `da_stability.py:88-89` | 严重 |
| 12 | 嫁接诱饵按索引截断，无锚点重排（面板肽长 8–13） | `dataset.py:119` | 严重 |
| 13 | 嫁接静默回退但 resname 照写 → 身份/几何不一致样本 | `dataset.py:79-97, 125` | 严重 |
| 14 | 位移诱饵沿固定全局 +z，兜底静默吞掉失败 | `dataset.py:146`; `graph.py:360-362` | 严重 |
| 15 | 诱饵从未被证明不结合 → 0.933 与研究目标无桥梁 | AUDIT 自认，README 仍正面引用 | 严重 |
| 16 | 面板冗余：188 候选 = 59 受体 / 25 MHC；347 FASTA = 73 受体簇 | 实测 | 严重 |
| 17 | 序列伪影（NA/GS/MG 前缀）当天然序列提交，影响链分类 | 53+279/347 文件 | 严重 |
| 18 | 图非派生产物：消融序列与 0.644 硬编码 | `make_figures.py:101-111, 363` | 严重 |
| 19 | 无任何同数据集的已发表工具对照 | `benchmarks.md:1240-1246` | 严重 |
| 20 | 无随机肽基线 → 绝对分数无法解释 | 全项目 | 严重 |
| 21 | 复现产物全缺（5 个 per-seed 表、checkpoint、dataset.pt） | `data/processed/` | 严重 |
| 22 | git 仅 1 个 squash commit，AUDIT 引用的提交号全不存在 | 数据仓 | 严重 |
| 23 | 测试盲区正好覆盖所有关键路径（6 项 skip 全是数据相关） | `tests/` | 严重 |
| 24 | MHC 预筛 n=50+50，未报 CI[42%,70%]，外推 465 万肽 | paper2 §2.4 | 中等 |
| 25 | 计划书 80–90% 灵敏度无任何实验支撑 | `committee_proposal.md` | 中等 |
| 26 | 湿实验/文库/动物全未做，paper2 仍称"Dry-Wet Loop" | paper2 标题与摘要 | 中等 |
| 27 | 0.5 门限自设但叙事当客观标准（"判决不变""危机解除"） | README:301-302,377 | 中等 |
| 28 | CI 下界 0.495 < 0.5 仍判 STABLE；n=21 | `benchmarks.md:641-643` | 中等 |
| 29 | best-of-3 种子报告（0.920/0.741/0.692 → 只报 0.933） | paper1 L165-167 | 中等 |
| 30 | 种子过滤规则事后制定，无预注册标准 | paper1 §2.6 | 中等 |
| 31 | "单调改善"跳过 NaN 退化点，且与 README:305-308 自相冲突 | paper1 L462 | 中等 |
| 32 | Fisher 合并非独立，"family-stratified 版"从未实现 | `rank_stats.py:82-90` | 中等 |
| 33 | 面板规模四个口径并存（154/172/188/347） | README:137,302,317 | 中等 |
| 34 | `kn5_submission_list.tsv` 列名与内容仍相反 + 重复行 | 实测 5brz/5bs0/5eu6 | 中等 |
| 35 | AUDIT 给的测试入口路径不存在（`src/meta_acdc/tests`） | `AUDIT.md:230` | 轻微 |
| 36 | README 称 16 项测试，实际 23 项 | `README.md:519` | 轻微 |
| 37 | 撤回声明只在文档顶部，正文段落级仍正面引用已撤回数字 | README:345,349,378 | 中等 |
| 38 | Cameron/Linette/Morgan 三篇并列引用，Morgan 讲的是 MAGE-A12 非 titin | README:65, paper1 L87 | 轻微 |
| 39 | ZeroBind 等 bib 条目标 UNVERIFIED / DOI 占位 `xxxxx` | `references-al-meta.bib:52,97` | 轻微 |

---

## 13. 如果要继续，最小可行修复路径

按依赖顺序，前面的不做，后面的都无意义。

### 第一优先：让结果可追溯（否则一切讨论都是空的）

1. **交付 `prediction_scores_ensemble.tsv` 与 5 个 `da_scores_seed*_ensemble.tsv`**，以及产生它们的模型 checkpoint 和 `egnn_dataset.pt`。
2. **恢复真实 git 历史**，或至少给出一张"数字 → 提交 → 运行日志"的对照表。
3. **解决第 1 节的矛盾**：明确 kn11 与 relabeled 表哪个是当前口径，废弃另一个并从所有文档中删除其数字。
4. **把图改成派生产物**：`fig1_ablation`、`fig6_domain_adapt`、`figp2_*` 从结果表读取，删除所有硬编码数字。

### 第二优先：修掉真正影响结论的代码缺陷

5. `dataset.py:260,270` 与 `train_domain_adapt.py:132,145` 改用 `AA3TO1`（同 #2 的修法，四处）。
6. `train_domain_adapt.py` 接入 `train_egnn.pdb_group()` + `family_map`，加 `--group-by family`，三个 DA 脚本同步传参。**这是 0.644 有效性的前提。**
7. `simulate_al.py` 加 `import math`；补 `meta_acdc.data` 包（或从代码里移除对它的依赖并撤回相关结论）。
8. `da_stability.py` 的聚类键改为 TCR 链序列（或 family id），不用 PDB id。
9. `make_decoy` 改为：长度不匹配时按 MHC 锚点位置对齐或直接拒绝；`_transplanted_atoms` 回退时抛错而非静默返回，且不写供体 resname。
10. `make_displaced_peptide_decoy` 改为沿"肽质心 → MHC 沟槽质心"的方向向量平移，并断言肽仍在界面选择内；去掉 `graph.py:360` 的静默全复合物兜底。
11. `clinical_gold_standard.tsv` 的 PMID 23999400/23999405 改为 23926201；`kn11` 的 `tcr` 列改为 MAG-IC3，或换用真实 a3a 结构（若不存在，则明确声明"评的是同源 ImmTAC 受体，非临床受体"）。
12. KVAKELVHFL 从 A3A/5brz 案例中移除，归回它自己的 A\*02:01 案例。
13. **重跑 QC 闸门后的全流程**，废弃所有闸前口径数字（包括 0.644/188 和整张 kn11 表）。

### 第三优先：补上使结论可解释的对照

14. **随机肽基线**：给模型 1000 条随机 9 肽，报告分数分布。没有这个，任何绝对分数都无意义。
15. **分数分辨率下界**：用同一 job 的 5 个 AF3 模型 + 5 个种子，给出"可分辨的最小 delta"。若该值 > 0.1，则第 2 节说明所有临床解读都必须撤回。
16. **至少一个已发表工具的同数据集 head-to-head**（ERGO-II 或 NetTCR），否则删除"序列天花板"论证。
17. **诱饵有效性验证**：取 20 个嫁接诱饵，用 IEDB/文献查是否有实验结合证据。若有任何一个是已知结合肽，"诱饵=负样本"的前提即告破产。
18. **冻结面板**：确定唯一的候选集合（建议按受体去重后的 59 个），删除所有冗余口径，重算全部指标。

### 第四优先：表述收缩

19. 段落级标注 `[WITHDRAWN]`，不能只靠文档顶部横幅。
20. 删除"金标准""判决""危机解除"这类暗示客观验证的措辞，改为"自设阈值下的复现性检查"。
21. paper2 标题与摘要去掉 "Dry-Wet Loop"，明确定位为纯计算方法研究。
22. 报告 best-of-N 时同时给出全部种子的数值与中位数。
23. 临床案例分析明确标注 **n=1 受体**，删除一切暗示一般性的措辞。

---

## 14. 值得肯定的部分

审查不应只列问题。这个项目有几处做法在同类工作中是罕见的：

1. **主动撤回自己全部核心数字**。README 顶部的撤回声明明确列出"0.933、0.644、0.515、FATAL-mimicry、2 倍召回、179/281 一律不得引用"。多数项目在发现同类缺陷时会选择静默修改。
2. **阴性结果被完整记录**：v9 首训 0.988 自查为泄漏后回退；v10/v11 静电特征两代阴性并明确"方向重定向"；KN-8 蛋白组清单零富集如实报告；ILAKFLHWL"部分漏检"如实报告；Kimmtrak"安全"判读在归因修复后主动撤销（README:352-355）。
3. **Rashomon 效应的方法学观察本身有价值**。"同性能模型给出互不相关的排名"是真实现象，"外推打分必须先过跨实例稳定性检验"这个原则可推广。问题不在于提出这个门，而在于门限无外部校准、且通过后被当成了正确性证据。
4. **CIF 覆盖事故、job_id 判重漏判、pdb 标签写反**等数据事故都被记录并追查了根因，`ensemble_mean.py` 还加了"区分 std=0 与无法估计"的正确处理。
5. **`clinical_scan.py` 重写后的口径是对的**：不再用 0.1 分差判"统计不可区分"，hazard 仅在证据不含 unverified 时输出，无结构行显式报告不评分。代码层面的修复方向是正确的——问题是交付的数据文件没有跟着重跑。

**核心矛盾就在这里**：项目的方法学自觉（知道自己哪里错了）明显强于它的工程落实（改完没有重跑、重跑没有交付、交付了没有对账）。第 4 节那 8 条"声称已修实际未修"，全部是这一类问题。

---

## 15. 最终判定

- **代码层面**：`graph.py`、`egnn.py`、`train_egnn.py`、`da_stability.py`、`ensemble_mean.py`、`score_predictions.py`、`clinical_scan.py` 的修复是真实的、方向正确的，注释里对缺陷成因的记录质量很高。但 `dataset.py` 和 `train_domain_adapt.py` —— **恰恰是产出全部报告数字的两个文件** —— 仍带着审计声称已修的缺陷。
- **结果层面**：**没有任何一个数字是可引用的。** 这不是因为审计过于严苛，而是因为：(a) 唯一交付的结果表与临床扫描表矛盾；(b) 产生报告数字的输入（模型、数据集、per-seed 打分）全部缺失；(c) 被解释的效应量小于测量噪声；(d) 临床案例的受体身份标错且 QC 未落实。
- **诚实性层面**：作者的自我审计是真诚且高质量的，这一点应当明确肯定。但审计的**执行完成度**被系统性高估了——28 条中至少 8 条的"已修"状态与交付代码不符，且这些不符之处没有被任何测试捕获。

**如果这是一个待投稿的稿件，我的建议是 Major Revision 起步，且必须先完成第 13 节第一、二优先级（结果可追溯 + 真正影响结论的代码缺陷），否则不具备评审条件。** 目前它更接近一份高质量的研究过程记录，而非一份可评审的研究成果。
