# 论文图表审阅材料

入口：MANUSCRIPT.md；方法与证据：METHODS_AND_EVIDENCE.md。

4张图均有PNG和可编辑SVG，已逐张目视检查。figure_data.csv保存区间图数值；family_detail.csv完整保留11家族；leave_one_family_out.csv为事后描述性影响检查，不替代主结果。

本轮从冻结输入重算组成和图大小基线，从15套真实标签矩阵重建逐家族表现。图大小两口径点估计及区间、逐家族主差值与历史JSON一致。其余图引用上轮已独立复算JSON，输入哈希见figure_verification.json。本轮未重训、未修改原图准入，也未访问独立验收。

复现需要相邻的audit_followup_20261025和audit_handoff_20261024目录及本机WSL冻结项目/环境。单独此包不含原始图或矩阵，不能声称可在干净机器上直接重演。build_figures.py拒绝覆盖既有figure_data.csv，请在新的同级目录复制脚本后运行：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 /home/administrator/docs/meta-acdc/.venv/bin/python /mnt/c/Users/Administrator/Documents/ChatGPT/meta-acdc/NEW_MATERIALS/build_figures.py
```

审阅状态：数值与图像检查完成，文献定位和投稿格式尚未完成，不是可直接投稿终稿。冻结时间的表述来自既有报告，未独立重建服务器时间线。
