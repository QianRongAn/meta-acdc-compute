#!/usr/bin/env python3
"""Append a ColabFold run record to docs (idempotent, marker-guarded).

Usage:
    .venv/bin/python scripts/record_colabfold_run.py run2 2026-10-14 85
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def load_scores() -> dict[str, float]:
    d = {}
    with open(REPO / "data/processed/colabfold_scores_ensemble.tsv",
              newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("score"):
                d[r["job_id"]] = float(r["score"])
    return d


def main() -> int:
    run, date, expected = sys.argv[1], sys.argv[2], sys.argv[3]
    jobs_file = Path(f"/tmp/colabfold_{run}_jobs.txt")
    jobs = [l.strip() for l in jobs_file.read_text().splitlines() if l.strip()]
    scores = load_scores()
    vals = sorted(scores[j] for j in jobs if j in scores)
    n = len(vals)
    if n == 0:
        print("no scores found; nothing recorded")
        return 1
    mean = sum(vals) / n
    lo, hi = vals[0], vals[-1]
    n_low = sum(1 for v in vals if v < 0.1)
    marker = f"ColabFold run{run[-1]} 回收({n}"

    # ---- benchmarks.md ----
    bm = REPO / "docs/benchmarks.md"
    txt = bm.read_text(encoding="utf-8")
    if f"## 2026-10-14:ColabFold 本地批量结构 {run} 回收" not in txt:
        sec = f"""
## {date}:ColabFold 本地批量结构 {run} 回收({n}/{expected})

- 与 run1 同参数(AF2-multimer,单模型,max-msa 128:256,3 recycle);
  `scripts/fetch_colabfold_run2.py` 自动轮询 → 远端打包(单 tar,tar -tzf 校验)
  → SFTP 下载 → 导入打分(5-seed DA,mask-plddt,T=1.0,与 da_rescore 同口径)。
- 打分:{n}/{expected} 成功,均值 **{mean:.3f}**,范围 {lo:.4f}–{hi:.4f},
  其中 <0.1 的 {n_low} 个。
- 沿用 run1 结论:ColabFold(AF2-multimer)结构与 AF3 不可互换用于 DA 打分,
  本批分数仍只作结构 QC/粗筛层,**不进主风险面板**。
- 明细:`data/processed/colabfold_scores_ensemble.tsv`、
  `colabfold_{run}_import_report.tsv`。
"""
        anchor = "## 待做基准"
        txt = txt.replace(anchor, sec + "\n" + anchor, 1)
        bm.write_text(txt, encoding="utf-8")
        print(f"benchmarks.md: appended {run} section")

    # ---- KEY-NODES.md ----
    kn = REPO / "docs/milestones/KEY-NODES.md"
    txt = kn.read_text(encoding="utf-8")
    if f"ColabFold 的 **{run}({n} 个)回收完成" not in txt:
        line = (f"- [x] **2026-10-14 会话**:待传队列本地 ColabFold **{run}"
                f"({n}/{expected} 个)回收完成**——自动拉取+导入+5-seed DA 打分"
                f"(均值 {mean:.3f});沿用 run1 结论:AF2-multimer 结构与 AF3 "
                f"不可互换用于 DA 打分,只作 QC/粗筛层,不进主风险面板\n")
        anchor = "- 项目仓库:`/home/administrator/docs/meta-acdc/`"
        txt = txt.replace(anchor, line + anchor, 1)
        kn.write_text(txt, encoding="utf-8")
        print("KEY-NODES.md: appended session line")
    return 0


if __name__ == "__main__":
    sys.exit(main())
