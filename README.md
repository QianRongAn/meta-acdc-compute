> 2026-10-18 更新：[临床原文、工程构型和结局整理](docs/CLINICAL_PRIMARY_CURATION_20261018.md)。临床资料仍按患者/组合保留；新SPR开发集51条，临床候选队列独立预留，旧版本数据量和结果仅作历史。

> **最新处置：[N2／N8 原子数据修复](docs/STRUCTURE_ATOMIC_REPAIR_20261017.md)。旧 588 张合成控制已退役；168 份无标签天然结构资料核验通过。新控制仅 5 对，未覆盖校准／测试分区，结构重训仍不放行。**

# Meta-ACDC：TCR–pMHC 实验数据研究与结构合成控制原型

**2026-10-17 当前状态：暂停旧结构数据重训；尚无获准使用的真实结合预测或临床风险模型。**

[N1–N14 逐项复核与修复](docs/N_REVIEW_20261017.md)是当前状态依据。原有“序列→形状→化学”的数字比较、稳定性放行和主动学习增益已撤回，不构成同任务方法比较。原 README、基准、稿件和日志完整保存在[历史归档](docs/history/20261017-before-N-review/README.md)，附原文件哈希。

## 当前三个研究终点

| 任务 | 可用依据 | 当前结论 |
|---|---|---|
| 结构合成控制判别 | 天然结构与指定规则生成的扰动 | 不是结合预测；旧诱饵有几何伪影，须重建验收 |
| 实验亲和力与功能研究 | SPR 亲和力记录；MIRA 激活／未激活观察 | 可以按各自终点研究；不能混作同一二分类真值 |
| 临床脱靶风险 | 病例、受体构建体线索、试验安全性信息 | 尚无足够独立且可归因的受体验收集 |

## 本轮新增可复现发现

对旧 v4 的 888 张图从原始原子重建检查：300 个天然图中 9 个缺失界面重原子；288 个嫁接图中 146 个出现肽／环境跨链原子距离不超过 1.2 Å；300 个固定方向位移图全部触发该重叠筛查。未把这些统计当作独立生物学样本数，也未声称它们证明模型实际使用了伪影。完整明细见 [atom audit](docs/n_review/structural_atom_audit.json)。

当前构建与验收代码会拒绝缺失界面原子和触发构造筛查的图。筛查通过不证明合理 rotamer、结合或不结合。旧数据、checkpoint 和日志均未覆盖。旧放行决定已被本轮发现撤销，代码指纹检查会拒绝旧 manifest。

## 实验数据库工作

[已完成的数据收集和新基线实验](docs/DATABASE_BINDING_VALIDATION.md)：外部保留排除后 91 条 SPR 亲和力记录；10,000 条 MIRA 功能观察；外部功能排名未通过有效性验收。该外部面板已使用，后续用于模型选择或训练时必须另留新的未使用验收来源。

[临床独立受体可行性核查](docs/CLINICAL_RECEPTOR_FEASIBILITY.md)：新增取得 TRAIT 公开临床表，正在按原文追溯工程构建体和结局；不同名字、PDB、试验编号或 CDR3 配对均不直接等于独立临床样本。

## 可复现检查

复用现有 Python 环境，不重复安装。在仓库目录运行：

```bash
PYTHONPATH=src /home/administrator/docs/meta-acdc/.venv/bin/python -m unittest discover -s tests -v
PYTHONPATH=src /home/administrator/docs/meta-acdc/.venv/bin/python scripts/audit_structural_artifacts.py \
  --release-dir data/releases/20261017-structural-v4 --out NEW_ATOM_AUDIT.json
```

输出必须使用新路径。训练入口要求通过源文件、代码、数据哈希及冻结切分检查的独立版本。当前没有新结构版本获准重训。hy4 仍按用户要求暂不作为验收环节，本轮未新增调度或派发训练任务。
