"""Apply the coordinate-based relabel map to legacy score tables.

The legacy `prediction_scores*.tsv` were produced before the 2026-09-23
labelling fix (kn5_submission_list.tsv had tcr_a/tcr_b swapped), so their
job_ids carry the old, wrong pdb attributions. This rewrites the `job_id`
column using data/processed/af3_relabel_map.tsv so figures/dashboard show the
corrected TCR families.

Usage:
    .venv/bin/python src/meta_acdc/structure/relabel_scores.py \
        --scores data/processed/prediction_scores.tsv \
                 data/processed/prediction_scores_ensemble.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def load_map(path: Path, instance: str | None) -> dict[str, str]:
    m: dict[str, str] = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if instance and row["instance"] != instance:
                continue
            m[row["old_job"]] = row["new_job"]
    return m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path, nargs="+", required=True)
    ap.add_argument("--map", type=Path,
                    default=Path("data/processed/af3_relabel_map.tsv"))
    ap.add_argument("--instance", default="af3_predictions",
                    help="which instance's map to apply ('' = all instances)")
    ap.add_argument("--out-suffix", default=".relabeled")
    args = ap.parse_args()

    mapping = load_map(args.map, args.instance or None)
    for src in args.scores:
        if not src.exists():
            print(f"skip missing {src}", file=sys.stderr)
            continue
        with open(src, newline="", encoding="utf-8", errors="replace") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        n = 0
        for r in rows:
            jid = r.get("job_id", "")
            base, sep, tag = jid.rpartition("_model_")
            if sep and base in mapping:
                r["job_id"] = mapping[base] + sep + tag
                n += 1
            elif jid in mapping:
                r["job_id"] = mapping[jid]
                n += 1
        dst = src.with_name(src.stem + args.out_suffix + src.suffix)
        with open(dst, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"{src.name}: relabeled {n}/{len(rows)} -> {dst.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
