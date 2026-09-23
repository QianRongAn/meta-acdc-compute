# SOP:AI 驱动的 ACDC 文库设计(Deliverable 7)

> 本 SOP 定义 **AI 引导的 ACDC 肽库设计**、**实验规模计算**与**模型标定**
> 三个工作流,作为实验室的长期技术规范。计算侧工具全部在本仓库
> (`meta_acdc` 包 + CLI,见 [`API.md`](API.md));湿实验执行由实验室完成。
>
> 复现环境见 [`SOP.md`](SOP.md);数据见 [`DATASET.md`](DATASET.md)。

---

## 1. AI 引导的 ACDC 文库设计

### 1.1 目标

给定一个**治疗 TCR**(α/β CDR3 + MHC 限制),从人蛋白组(>10⁷ 肽)中
选出一个 **~5 万条肽的合成文库**,最大化"包含真实脱靶"的概率,同时
用最少的湿实验轮次完成。

### 1.2 三级漏斗(推荐流程)

```
人蛋白组 (20,431 蛋白, ~465 万唯一 9/10-mer)
   │  ① MHC 呈递预筛(MHCflurry,本地)
   ▼
候选池 ~465 万(阈值 top-2%,保留 56% 已知结合肽 / 压缩 ~50×)
   │  ② 序列代理 EIG 排序(+ 低复杂度过滤)
   ▼
top-N 清单(默认 5 万)
   │  ③ 结构重排(AF3 提交 + 稳定域适应 EGNN 打分)
   ▼
最终文库(带风险分与排名)
```

**为什么三级而非直接结构打分**:结构预测成本高;①先按"能否被呈递"
压缩 ~50×(已本地验证),③只对压缩后的清单做结构重排。

### 1.3 命令

```bash
# ① 预筛:蛋白组 -> 5M 池
curl -sL 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(organism_id:9606)%20AND%20(reviewed:true)' \
    -o data/raw/human_sprot.fasta
meta-acdc-select --stage prefilter

# ② EIG 排序(+复杂度过滤,剔除 Q/E 低复杂度伪影)
meta-acdc-select --stage rank --tcr <CDR3_beta> --max-aa-fraction 0.5

# ③ 结构重排:为治疗 TCR 结构 + 候选肽生成 AF3 提交批
python src/meta_acdc/structure/af3_scan_batch.py \
    --pdb <TCR结构id> --peptides data/processed/kn8_top50k_<tag>.tsv \
    --n 200 --out data/processed/af3_scan_<tag>
#    (用户在 alphafoldserver.com 手动提交 -> folds_* 下载)
#    导入 + 稳定 DA 打分 + 排名
meta-acdc-native-manifest --src folds_<tag> --out <native_dir> \
    --manifest <manifest> --candidate-out data/raw/af3_predictions
meta-acdc-score --cifs data/raw/af3_predictions \
    --model data/processed/da_seed0.model.pt --mask-plddt --out <scores>
```

### 1.4 质量闸门(报告前必须通过)

| 闸门 | 工具 | 判据 |
|---|---|---|
| **跨实例稳定性** | `meta-acdc-stability` | mean pairwise Spearman ≥ 0.5 |
| 退化实例 | 同上 | 常数输出实例自动剔除 |
| 结构可靠性 | `structure/af3_qc.py` | ipTM ≥ 0.75 且 TCR Cα ≥ 400 |
| 复杂度 | `select_candidates` | 最高频残基占比 ≤ 0.5 |
| 标定 | `models/calibrate.py` | ECE 报告 |

**规则**:任何排名在通过稳定性闸门之前,只作为"研究性观察",不得写入
订购清单或安全结论。

---

## 2. 实验规模计算

### 2.1 预算公式

单轮 ACDC 实验成本 ≈ 文库肽数 × 每肽合成/递呈成本;总预算 =

```
总实验量 ≈ 轮数 R × 批大小 B × 重复数 K
```

### 2.2 批量与轮次(基于本项目消融,固定预算 5,964 样本)

| 批大小 B | 轮数 R | 最终阳性召回 |
|---|---|---|
| **3,600** | 2 | **24.5% ± 2.3** |
| 1,800 | 3 | 19.7% ± 1.4 |
| 900 | 7 | 20.8% ± 0.7 |
| 600 | 10 | 21.9% ± 0.5 |
| 300 | 20 | 21.6% ± 1.0 |

**结论**:**大批少轮更优**(与文献"小批多轮"相反);弱模型下每批需
足够样本才能学到区分信号,**批量存在下界**。推荐 **2–3 轮、大批**
(≥1,800–3,600/批)。

### 2.3 采集函数增益(阈值依赖)

| 预筛严格度 | EIG 增益(低采样区) |
|---|---|
| 宽松(0.5 阈值池) | **1.76×**(7.4% 采样) |
| 严格(0.9 阈值池) | ≈1.0×(增益消失) |

**规则**:增益是**预筛阈值的函数**;越宽松的预筛,EIG 越有用。报告须
给出增益-阈值曲线,而非单一数字。熵采样≈随机(不推荐);EIG > 方差。

### 2.4 规模估算示例

- 若目标:在预算内达到 ~25% 阳性召回 → 采 **2 轮 × 3,600/批**。
- 若需覆盖更宽表位空间 → 增大预筛池(放宽阈值)并增加批数,接受
  单轮增益下降。
- 每次新增实验后**重训模型**(闭环),下一轮 EIG 更准(数字孪生协议)。

---

## 3. 模型标定工作流

### 3.1 步骤

```bash
# 1) 域适应训练(5 seed,含 val AUROC 元数据)
meta-acdc-train-da --af3 <native_dir> --af3-manifest <manifest> \
    --seeds 0 1 2 3 4 --epochs 150

# 2) 温度标定(晶体留出集,最小化 BCE)
python -m meta_acdc.models.calibrate --model data/processed/da_seed0.model.pt \
    --out data/processed/da_seed0.temperature.json

# 3) 打分(应用温度 -> 标定概率)
meta-acdc-score --cifs <cifs> --model <ckpt> --mask-plddt \
    --temperature <T> --out <scores>
```

### 3.2 报告规则

- **必须**报告:多实例集成分(3-seed 均值)+ 跨实例稳定性(mean r + CI);
- **必须**报告:标定后的概率与 ECE(当前 ECE 0.02–0.06);
- **不得**报告:单实例分数作为结论(仅作轶事);
- **不得**报告:未过稳定性闸门的排名。

### 3.3 已知边界(诚实记录)

- 模型在**诱饵任务**上强(嫁接 AUROC 0.93),在**真实 AF3 结构**上
  排名需 ≥~144 天然正样本才稳定;
- 相容性判定是**分布式**的(无单点热点),消除误伤可能需多位点协同;
- 绝对分**不是**结合概率,仅家族内排序可信;阈值需真实标签标定。

---

## 4. 闭环协议(数字孪生 → 真实)

1. **设计**(第 1 节)→ 输出文库;
2. **湿实验**(实验室):ACDC 平台合成/递呈/FACS/NGS → 回读阳性;
3. **更新**:用回读数据替换模拟器 oracle(`simulate_al.py`),重训模型;
4. **迭代**:下一轮 EIG 基于更新后的模型 → 定向采样盲区。

模型侧重训**零改动**(接口固定);详见
[`milestones/KEY-NODES.md`](milestones/KEY-NODES.md) 技术路线表。

---

## 5. 交付物清单(对应计划书)

| 计划书交付物 | 本仓库对应 |
|---|---|
| D3 软件套件 | `meta_acdc` 包 + CLI + [`API.md`](API.md) + Docker |
| D4 交互看板 | `meta-acdc-dashboard`(TCR-Safety-Radar) |
| D5 策展数据集 | [`DATASET.md`](DATASET.md) |
| D7 SOP | **本文件** |
| D1/D2 论文 | [`paper1_manuscript.md`](paper1_manuscript.md)、[`paper2_manuscript.md`](paper2_manuscript.md) |
| D6 验证肽库 | 湿实验(实验室) |
