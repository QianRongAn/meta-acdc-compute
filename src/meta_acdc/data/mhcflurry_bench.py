"""MHCflurry pre-filter benchmark (paper 2, experiment 1 — local, no queue).

Same protocol as netmhcpan_bench.py: 50 known HLA-A*02:01 9-mer binders
(IEDB positives) + 50 shuffled decoys -> presentation scores -> recall and
compression at percentile thresholds.

Usage:
    .venv/bin/python src/meta_acdc/data/mhcflurry_bench.py --n 50
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

from mhcflurry import Class1PresentationPredictor

from meta_acdc.data.netmhcpan_bench import load_binders, shuffle_peptides

ALLELE = "HLA-A*02:01"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--iedb", type=Path, default=Path("data/processed/iedb.clean.tsv"))
    ap.add_argument("--n", type=int, default=50)
    args = ap.parse_args()

    binders = load_binders(args.iedb, ALLELE, args.n)
    decoys = shuffle_peptides(binders)
    print(f"loading MHCflurry models (first run downloads weights) ...",
          flush=True)
    predictor = Class1PresentationPredictor.load()

    print(f"predicting {len(binders) * 2} peptides ...", flush=True)
    df = predictor.predict(
        peptides=binders + decoys,
        alleles=[ALLELE],
        verbose=0,
    )
    scores = df["presentation_score"].to_numpy()
    s_b = scores[: len(binders)]
    s_d = scores[len(binders):]

    # percentile thresholds (approximating %rank semantics):
    # threshold = quantile of the decoy distribution (null reference)
    import numpy as np
    for q in (0.5, 2.0, 5.0, 10.0):
        null_q = np.quantile(s_d, 1 - q / 100)
        recall = (s_b > null_q).mean()
        print(f"decoy-quantile top {q}%: binder recall {recall:.2f} "
              f"(null threshold {null_q:.3f})")

    # rank-based summary
    all_s = np.concatenate([s_b, s_d])
    ranks = (all_s[:, None] > all_s[None, :]).sum(0) + 1
    mean_rank_b = ranks[: len(binders)].mean()
    mean_rank_d = ranks[len(binders):].mean()
    print(f"mean rank: binders {mean_rank_b:.0f} vs decoys {mean_rank_d:.0f} "
          f"of {len(all_s)} (lower = better)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
