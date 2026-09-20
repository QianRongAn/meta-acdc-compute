"""Dataset statistics report for the unified TCR-pMHC datasets (KN-3).

Usage:
    python stats_report.py --data data/processed/*.clean.tsv --out docs/dataset-report.md
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def load(path: Path) -> list[dict]:
    rows = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
    return rows


def top(counter: Counter, n: int = 10) -> str:
    items = counter.most_common(n)
    total = sum(counter.values())
    return ", ".join(f"{k} ({v}, {100 * v / total:.1f}%)" for k, v in items)


def report(paths: list[Path], out: Path) -> str:
    sections = ["# 数据集统计报告\n"]
    for p in paths:
        rows = load(p)
        src = rows[0]["source"] if rows else p.stem
        n = len(rows)
        pep_len = Counter(len(r["peptide"]) for r in rows if r["peptide"])
        alleles = Counter(r["mhc_allele"] for r in rows if r["mhc_allele"])
        labels = Counter(r["label"] for r in rows if r["label"])
        scores = [float(r["score"]) for r in rows if r["score"]]
        n_scored = len(scores)

        sections.append(f"## {src}({n:,} 条记录)\n")
        if any(r["cdr3_beta"] for r in rows):
            cdr3_len = Counter(len(r["cdr3_beta"]) for r in rows if r["cdr3_beta"])
            sections.append(f"- 独特 CDR3-β: {len({r['cdr3_beta'] for r in rows if r['cdr3_beta']}):,}")
            sections.append(f"- CDR3-β 长度分布(氨基酸数): {top(cdr3_len, 5)}")
            sections.append(f"- 独特肽: {len({r['peptide'] for r in rows if r['peptide']}):,}")
        sections.append(f"- 肽长分布(氨基酸数): {top(pep_len, 6)}")
        sections.append(f"- MHC 等位基因数: {len(alleles):,}")
        sections.append(f"- 最常见等位基因: {top(alleles, 8)}")
        sections.append(f"- 标签分布: {top(labels, 4)}")
        if n_scored:
            sections.append(f"- score 范围: [{min(scores):.2f}, {max(scores):.2f}](n={n_scored:,})")
        sections.append("")

    text = "\n".join(sections)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    return text


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", nargs="+", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("docs/dataset-report.md"))
    args = ap.parse_args()
    print(report(args.data, args.out))
    return 0


if __name__ == "__main__":
    main()
