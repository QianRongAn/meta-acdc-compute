"""All-planned-seed reproducibility diagnostic, never a binding-accuracy gate.

A constant output means absent ranking signal, not demonstrated uncertainty.
Seed failures and missing jobs block a passing verdict; no survivor filtering.
"""
from __future__ import annotations
import argparse
import csv
import json
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr


def load_jobs(path: Path) -> dict[str, float]:
    jobs = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            job = row["job_id"]
            if not job or job in jobs:
                raise ValueError("empty or duplicate job ID: " + job)
            value = float(row["score"])
            if not np.isfinite(value):
                raise ValueError("non-finite score: " + job)
            jobs[job] = value
    if not jobs:
        raise ValueError("empty seed result: " + str(path))
    return jobs


def analyse(per_seed, families, *, expected_seeds, resamples=2000):
    if expected_seeds < 2 or len(per_seed) != expected_seeds:
        raise ValueError("all predeclared seeds, including failures, must be accounted for")
    common = sorted(set.intersection(*(set(s) for s in per_seed)))
    union = set.union(*(set(s) for s in per_seed))
    result = dict(n_planned=expected_seeds, n_received=len(per_seed),
                  n_common_jobs=len(common), n_union_jobs=len(union),
                  coverage_complete=all(set(s) == union for s in per_seed),
                  task="ranking_reproducibility_only", threshold=0.5,
                  threshold_validated_against_binding=False)
    if len(common) < 5:
        return dict(result, verdict="INSUFFICIENT_JOBS", mean_r=None)
    M = np.array([[s[j] for j in common] for s in per_seed])
    degenerate = [i for i, row in enumerate(M) if np.std(row) < 1e-6]
    pairs = []
    for i, j in combinations(range(expected_seeds), 2):
        r = None if i in degenerate or j in degenerate else float(spearmanr(M[i], M[j]).statistic)
        pairs.append(dict(seed_i=i, seed_j=j, r=r))
    valid = [x['r'] for x in pairs if x['r'] is not None]
    all_defined = len(valid) == len(pairs)
    result.update(degenerate_seeds=degenerate, pairs=pairs,
                  mean_r=float(np.mean(valid)) if all_defined else None,
                  conditional_nonconstant_mean_r=float(np.mean(valid)) if valid else None,
                  conditional_mean_is_not_primary=True)
    groups = {}
    for k, job in enumerate(common):
        if job not in families: raise ValueError("missing receptor family: " + job)
        groups.setdefault(families[job], []).append(k)
    result['n_families'] = len(groups)
    result['verdict'] = ('INCOMPLETE_COVERAGE' if not result['coverage_complete'] else
                         'DEGENERATE_SEED' if degenerate else
                         'INSUFFICIENT_FAMILIES' if len(groups) < 2 else
                         'STABLE_QUALITY_ONLY' if result['mean_r'] >= .5 else 'UNSTABLE')
    # Undefined resamples remain in the denominator; never silently discard them.
    if all_defined and len(groups) >= 2:
        rng = np.random.default_rng(0); keys = sorted(groups); boot = []; invalid = 0
        for _ in range(resamples):
            idx = [k for n in rng.integers(0, len(keys), len(keys)) for k in groups[keys[n]]]
            if any(np.std(row[idx]) < 1e-6 for row in M):
                invalid += 1; continue
            boot.append(float(np.mean([spearmanr(M[i,idx], M[j,idx]).statistic
                                       for i,j in combinations(range(expected_seeds),2)])))
        result['bootstrap'] = dict(requested=resamples, undefined=invalid, defined=len(boot),
                                   conditional_ci95=np.percentile(boot,[2.5,97.5]).tolist() if boot else None,
                                   interval_conditions_on_defined_resamples=True)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scores',type=Path,nargs='+',required=True)
    ap.add_argument('--release',type=Path,required=True)
    ap.add_argument('--seed-plan',type=Path,required=True,
                    help='pre-run JSON: {"score_files": ["seed0.tsv", ...]}')
    ap.add_argument('--out',type=Path,required=True)
    args = ap.parse_args()
    from meta_acdc.data.release import load_release, complex_id, sha256
    load_release(args.release)
    plan = json.loads(args.seed_plan.read_text())['score_files']
    received = [x.name for x in args.scores]
    if len(set(plan)) != len(plan) or sorted(plan) != sorted(received):
        raise ValueError('score files do not match the predeclared seed plan')
    if args.out.exists(): raise ValueError('refusing to overwrite stability report')
    family_by_pdb = {}
    for row in json.loads((args.release.parent/'structure_records.json').read_text()):
        key = complex_id(row['pdb'])
        if key in family_by_pdb and family_by_pdb[key] != row['family_id']:
            raise ValueError('ambiguous receptor family for PDB')
        family_by_pdb[key] = row['family_id']
    seeds = [load_jobs(x) for x in args.scores]
    jobs = set.union(*(set(x) for x in seeds))
    families = {j:family_by_pdb[complex_id(j)] for j in jobs}
    result = analyse(seeds, families, expected_seeds=len(plan))
    result.update(seed_plan_sha256=sha256(args.seed_plan),
                  score_files={str(x):sha256(x) for x in args.scores})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))
    return 0 if result['verdict']=='STABLE_QUALITY_ONLY' else 1


if __name__=='__main__':
    raise SystemExit(main())
