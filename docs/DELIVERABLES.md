# 计划书交付物对照表(Deliverables)

> 计划书(proposal §VIII)列出 7 项交付物。本表把它们映射到本仓库的
> 实际产物与状态。**干实验范围**:D1–D5、D7 可在本仓库完成;
> D6(物理肽库)属湿实验,由实验室执行。

## 状态总览

| # | 交付物 | 状态 | 主要产物 |
|---|---|---|---|
| D1 | 论文 1(方法学) | 草稿 | [`paper1_manuscript.md`](paper1_manuscript.md) |
| D2 | 论文 2(系统应用) | 草稿 | [`paper2_manuscript.md`](paper2_manuscript.md) |
| D3 | 软件套件 | **可用** | `meta_acdc` 包 + 11 CLI + [`API.md`](API.md) + Docker |
| D4 | 交互看板 | **可用** | `meta-acdc-dashboard`(TCR-Safety-Radar) |
| D5 | 策展数据集 | **可用** | [`DATASET.md`](DATASET.md) + `data/processed/*` |
| D6 | 验证肽库 | 湿实验 | 由实验室执行(技术路线见 KEY-NODES) |
| D7 | SOP | **可用** | [`SOP.md`](SOP.md) + [`SOP_ACDC_LIBRARY.md`](SOP_ACDC_LIBRARY.md) |

---

## D1 — 论文 1(Methodology)

**标题**:"Meta-TCR-GNN: A Geometric Meta-learning Framework for Few-shot
TCR Specificity Prediction"

**内容**:Cα 不可分诊断 → 侧链接触特征突破(0.500→0.933)→ 严格诱饵
协议 → 跨实例稳定性判决门 → 域适应规模-稳定性曲线(危机解除)。

**产物**:[`docs/paper1_manuscript.md`](paper1_manuscript.md);
图 `docs/figures/fig1–7`;数据 [`DATASET.md`](DATASET.md);
数字 [`benchmarks.md`](benchmarks.md)。

## D2 — 论文 2(Systematic Application)

**标题**:"An AI-driven 'Dry-Wet' Loop for Proteome-wide Safety Assessment
of Therapeutic TCRs"

**内容**:主动学习闭环(EIG 2.0× 增益 + 阈值依赖边界 + 批大小消融)→
MHC 预筛(50× 压缩)→ 蛋白组 5M 池 → 结构重排 → 临床金标准扫描
(致死模拟肽不放过)。

**产物**:[`docs/paper2_manuscript.md`](paper2_manuscript.md)。

## D3 — 软件套件(Software Suite)

- `meta_acdc` Python 包(数据 / 结构 / 模型 / 主动学习 / 看板)
- **11 个 CLI 入口**(`pyproject.toml [project.scripts]`):score /
  stability / ensemble / clinical / attribute / mutate / native-manifest /
  select / train-egnn / train-da / dashboard
- **API 文档**:[`docs/API.md`](API.md)
- **Docker**:`Dockerfile`(CPU,已构建验证,18 测试通过)
- 回归测试:`tests/`(19 项)+ `scripts/verify.sh`

## D4 — 交互看板(TCR-Safety-Radar)

- `meta-acdc-dashboard`(零依赖 Web)
- 面板:热度图 / 每 TCR 排名 / 临床金标准 / 交叉反应证据 / 结构可靠性
  闸门 / **域适应状态** / **提交进度** / **KN-11 临床扫描** / 诚信横幅

## D5 — 策展数据集(Curated Dataset)

- 结构↔特异性图谱(281 复合物,2,308 命中,1,517 交叉反应)
- VDJdb/IEDB 清洗表(脚本复现)
- 临床金标准、KN-5 清单、AF3 天然清单(144)
- 数据卡片:[`docs/DATASET.md`](DATASET.md)

> 注:计划书的"百万级 ACDC 验证条目"需实验室产出;当前为
> **公开数据策展 + 结构派生**版本(如实标注)。

## D6 — 验证肽库(物理库)— 湿实验

高-risk 脱靶肽的实体/质粒库,由实验室按
[`KEY-NODES.md`](milestones/KEY-NODES.md) 技术路线执行;本仓库提供
**候选清单与风险排序**作为输入。

## D7 — SOP(Technical Standards)

- [`SOP.md`](SOP.md):计算管线从零复现
- [`SOP_ACDC_LIBRARY.md`](SOP_ACDC_LIBRARY.md):AI 引导文库设计 /
  实验规模计算 / 模型标定 / 闭环协议

---

## 里程碑(KN-0..16)

见 [`milestones/KEY-NODES.md`](milestones/KEY-NODES.md)。截至 2026-09-23:
KN-1/3/4/6/7/7+ 完成;KN-5 进行中(用户提交);KN-11 首个结构空间结果;
KN-13 看板 v1;KN-14/15 论文草稿;KN-16 打包(Docker + SOP + API)。
