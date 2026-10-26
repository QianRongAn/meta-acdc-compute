# 2026-10-25 审计推进交付

先读 AUDIT_DRAFT.md：包含摘要、方法、结果、讨论、图表设计与收尾条件。是工作稿，不是投稿终稿。

本次实际执行：

1. 从冻结30套真实／置换模型的相似度和嵌入独立复算主要AUROC、排名Spearman及家族条件区间；与历史JSON在1e-12容差内一致。
2. 对旧15实例实际补跑训练折拟合的Procrustes；五折拟合／评价索引不交叉、覆盖331图。正交旋转下内部余弦变化小于3e-15。
3. 用精确索引进行事后跨域诊断，比较共同查询和同一负例集合，避免不同查询人口直接作差。
4. 整理项目中已有的最新跨来源取证摘要；本次没有重新跑42图的原子QC。

independent_metric_replay.json 和 cross_domain_diagnostic.json 保存完整数值、人口计数及源哈希。新增稳定性差值和跨域对比为探索性诊断，不包装成预注册确认性结果。

## 本机复现

需要原项目的冻结数据、矩阵及既有WSL环境，不是仅靠此压缩包即可重训的独立数据发布。没有安装新依赖。

从原项目根目录运行（请把NEW_OUTPUT换成不存在的新目录；不要覆盖本交付）：

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2
PY=/home/administrator/docs/meta-acdc/.venv/bin/python
DEST=/mnt/c/Users/Administrator/Documents/ChatGPT/meta-acdc/NEW_OUTPUT
mkdir "$DEST"
cp /mnt/c/Users/Administrator/Documents/ChatGPT/meta-acdc/audit_followup_20261025/{verify_metrics,cross_domain_diagnostic}.py "$DEST/"
PYTHONPATH=src:. "$PY" scripts/rashomon_fingerprint.py --dataset docs/fingerprint_retrain_20261023/input/fingerprint_dataset.pt --out-dir "$DEST/procrustes_train_only"
"$PY" "$DEST/verify_metrics.py"
"$PY" "$DEST/cross_domain_diagnostic.py"
```

跨域脚本使用与新目录同级的 audit_handoff_20261024/records_index_exact.json；verify_metrics.py中的ROOT固定为本机原项目路径。两个分析脚本拒绝覆盖已有结果。原对齐脚本副本存于sources，便于核对实际执行版本。

## 解释边界

旧Procrustes补跑仍针对缺陷旧模型，不能用于恢复旧0.946的科学主张。新模型的排名稳定性不依赖坐标对齐。条件区间不覆盖训练和研究来源总体不确定性。旧标签跨域诊断不能证明物理受体的域效应或泄漏效应。

没有训练新模型、调整数据准入、读取独立验收集或设置调度。原项目中其他工作的文件未修改。
