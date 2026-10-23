# 困难家族错误分析与2×2机制消融（2026-10-24，训练前冻结）

协议[protocol.json](protocol.json)与[字节锁](protocol.sha256)。先做[旧产物错误分析](diagnostics/README.md)，再用PDB原始COMPND/SEQRES核查标签：[证据](identity_audit/identity_audit.json)及原始头信息在identity_audit/headers/。发现19组完整αβ构建体被拆为44个观察序列ID：链角色正确不等于完整变量链身份正确。原始标签、图、权重和历史结果不改写，不按指标合并受体。

**准入被阻断**：当前数据只能用于观察序列哈希标签的故障诊断，不能用于生物学受体性能放行。即使新消融验证分数提高，也不会推广为真实受体预测增益。完整αβ变量链、标签/常数区边界及AF3构建体对应尚未逐条核验。新构建器默认拒绝把此类标签当已验证身份；必须显式--diagnostic-observed-labels，输出语义及未验证标记。

## 固定2×2开发实验

复用原折，无重洗牌：原折2/3/4共186图训练，原折1共86图开发验证，原折0的59图隔离于新模型之外。它们此前全部已暴露，不能称为新独立测试。验证困难口径只有21查询/3家族。实测独立验收集完全不读。

| 因素 | 对照 | 处理 |
|---|---|---|
| 身份特征 | 原逐图zscore | 仅还原原始AA one-hot列，其他变换不变 |
| 批次组织 | 随机图排序 | 家族/受体/视图分块后拼接，仍按64/64/58切批 |

4个组合×3个固定种子，均120epoch、同一网络与损失，不选epoch、不选种子。每图每epoch出现一次，批次数和每批图数不变；结构化批次同时改变正样本和困难负样本可见性，不能仅归因于负样本。统计训练中两类配对的实际覆盖率。

分别计算因素在另一因素固定时的配对差值及交互项；主口径为观察序列标签下的家族内AUROC，按家族等权和条件bootstrap报告。另排除被完整双链SEQRES证明属于同构建体的晶体负对，记录删对数、无定义查询及公共查询敏感性；这不证明剩余负对是真的不同生物学受体。

训练后仅在开发验证集做身份列置零、显式几何通道置零及肽池化置零的输入敏感性诊断。它们是分布外输入扰动，不能当作重新训练消融或生物学机制证据；几何置零仍保留图拓扑，不叫完全去几何。

参考首种子完整重放一次；保留全部模型和数值，包括负结果。任何生物学候选放行都被标签准入闸门禁止。

## 来源与解释

[RCSB 5E9D](https://www.rcsb.org/structure/5E9D)把D/I归为同一α实体、E/J归为同一β实体；本地原始SEQRES也相同，但观察坐标序列末端缺失不同，被哈希为两个类别。1AO7和1QRN均描述A6受体；当前观察序列还混入常数区差异。未用相似度自动修复这些跨PDB标签，不能只根据名称或短CDR3强行合并。核验引用：[1AO7](https://www.rcsb.org/structure/1AO7)、[1QRN](https://www.rcsb.org/structure/1QRN)。

冻结前36项测试通过，无跳过。复用现有CPU环境，无安装。运行：

```bash
PYTHONPATH=src:. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 /home/administrator/docs/meta-acdc/.venv/bin/python -u scripts/run_fingerprint_factorial.py --protocol docs/fingerprint_hard_family_20261024/protocol.json --out-dir docs/fingerprint_hard_family_20261024/run
```


## 已完成结果

[完整结果与限制](RESULTS.md)、[13套模型重计算验收](verification.json)、[身份核验队列](IDENTITY_REPAIR_QUEUE.json)。组合组开发验证0.4841，组成0.5321；不挑种子，不放行生物学候选。
