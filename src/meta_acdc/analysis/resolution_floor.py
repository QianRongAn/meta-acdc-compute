"""Score-resolution floor report (review item #2 machinery).

Given per-seed and per-model score tables plus a reference job, report the
DESCRIPTIVE dispersion and the deltas against the reference, so that a
future run can state the smallest delta the measurement can plausibly
resolve. This tool makes no significance or equivalence claims:

- delta < SD is a descriptive observation, NOT a test of "indistinguishable";
- the seed axis and the AF3-model axis are separate variance components and
  are reported separately, never summed into a fake "calibrated error";
- no clinical interpretation is emitted.

Usage:
    .venv/bin/python -m meta_acdc.analysis.resolution_floor \
        --seed-scores data/processed/da_scores_seed0_ensemble.tsv \
                       data/processed/da_scores_seed1_ensemble.tsv \
        --model-scores data/processed/prediction_scores.relabeled.tsv \
        --reference 5brz_EVDPIGHLY \
        --out data/processed/resolution_floor.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def load_score_table(path: Path, key_col: str = "job_id") -> dict[str, float]:
    rows = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if key_col not in r or "score" not in r or r["score"] == "":
                continue
            try:
                v = float(r["score"])
            except ValueError:
                continue
            rows.setdefault(r[key_col], []).append(v)
    return {k: vs[0] for k, vs in rows.items()}


def stats(values: list[float]) -> dict:
    import math
    n = len(values)
    if n == 0:
        return {"n": 0, "mean": None, "sd": None}
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / n if n > 1 else None
    sd = math.sqrt(var) if var is not None else None
    return {"n": n, "mean": mean, "sd": sd}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed-scores", type=Path, nargs="+", default=[],
                    help="one score table per training seed (job_id<TAB>score)")
    ap.add_argument("--model-scores", type=Path, default=None,
                    help="per-AF3-model score table (job_id<TAB>score); "
                         "multiple rows per job give the model axis")
    ap.add_argument("--reference", type=str, required=True,
                    help="job_id of the target peptide the deltas are "
                         "measured against")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.out.exists():
        raise ValueError("refusing to overwrite an existing report")

    seed_tables = [load_score_table(p) for p in args.seed_scores]
    seed_by_job: dict[str, list[float]] = {}
    if seed_tables:
        jobs = set.union(*(set(t) for t in seed_tables))
        for j in sorted(jobs):
            vs = [t[j] for t in seed_tables if j in t]
            if vs:
                seed_by_job[j] = vs

    model_by_job: dict[str, list[float]] = {}
    if args.model_scores is not None:
        with open(args.model_scores, newline="", encoding="utf-8",
                  errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("score") == "" or "job_id" not in r:
                    continue
                try:
                    model_by_job.setdefault(r["job_id"], []).append(
                        float(r["score"]))
                except ValueError:
                    continue

    jobs = sorted(set(seed_by_job) | set(model_by_job))
    if args.reference not in jobs:
        raise ValueError(
            f"reference job {args.reference} not present in the score tables")
    ref_seed = stats(seed_by_job.get(args.reference, []))
    ref_model = stats(model_by_job.get(args.reference, []))

    rows = []
    for j in jobs:
        s_seed = stats(seed_by_job.get(j, []))
        s_model = stats(model_by_job.get(j, []))
        delta = (s_seed["mean"] - ref_seed["mean"]
                 if s_seed["n"] and ref_seed["n"]
                 else s_model["mean"] - ref_model["mean"]
                 if s_model["n"] and ref_model["n"] else None)
        rows.append({
            "job_id": j,
            "seed_mean": s_seed["mean"], "seed_sd": s_seed["sd"],
            "seed_n": s_seed["n"],
            "model_mean": s_model["mean"], "model_sd": s_model["sd"],
            "model_n": s_model["n"],
            "delta_vs_reference": delta,
            "abs_delta_vs_reference": abs(delta) if delta is not None else None,
        })

    deltas = [r["abs_delta_vs_reference"] for r in rows
              if r["abs_delta_vs_reference"] is not None and
              r["job_id"] != args.reference]
    seed_sds = [r["seed_sd"] for r in rows
                if r["seed_sd"] is not None]
    model_sds = [r["model_sd"] for r in rows
                 if r["model_sd"] is not None]
    import statistics as st
    report = {
        "reference": args.reference,
        "reference_seed_stats": ref_seed,
        "reference_model_stats": ref_model,
        "n_jobs": len(jobs),
        "median_abs_delta_vs_reference":
            st.median(deltas) if deltas else None,
        "median_seed_sd": st.median(seed_sds) if seed_sds else None,
        "median_model_sd": st.median(model_sds) if model_sds else None,
        "n_jobs_with_abs_delta_below_own_seed_sd": sum(
            1 for r in rows
            if r["abs_delta_vs_reference"] is not None
            and r["seed_sd"] is not None
            and r["abs_delta_vs_reference"] < r["seed_sd"]),
        "rows": rows,
        "notice": (
            "DESCRIPTIVE ONLY. delta < SD is not a significance or "
            "equivalence test; the seed and AF3-model axes are separate "
            "variance components and are never summed into a calibrated "
            "measurement error. No clinical interpretation follows. Any "
            "claim about 'the smallest resolvable delta' must come from a "
            "pre-registered paired design, not from this table."),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "rows"},
                     indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
