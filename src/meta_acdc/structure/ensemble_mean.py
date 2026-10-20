"""All-planned-seed structural diagnostics with explicit score provenance.

Between-seed SD is not AF3 spread, calibrated uncertainty or clinical risk.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics
from pathlib import Path
from meta_acdc.structure.artifact_status import assert_not_withdrawn,validate_score_receipt

def load_seed(path):
    data={}
    with Path(path).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if not row.get('score'):continue
            job=row['job_id'];value=float(row['score'])
            if not job or job in data or not math.isfinite(value):raise ValueError('duplicate, empty, or nonfinite seed score')
            if row.get('qc')!='qc-pass':raise ValueError('score lacks QC-pass status')
            data[job]=value
    if not data:raise ValueError('empty seed table')
    return data

def combine(seeds,planned_jobs):
    if len(seeds)<2:raise ValueError('at least two predeclared seeds required')
    common=set.intersection(*(set(s) for s in seeds));planned=set(planned_jobs)
    if not set.union(*(set(s) for s in seeds))<=planned:raise ValueError('unplanned job in scores')
    degenerate=[i for i,s in enumerate(seeds) if len(s)<2 or statistics.pstdev(s.values())<1e-6]
    rows=[]
    for job in sorted(common):
        v=[s[job] for s in seeds]
        rows.append(dict(job_id=job,score=statistics.mean(v),std=statistics.pstdev(v),std_axis='between_seed_means',std_estimator='population',n_seeds=len(seeds),uncertainty_known='yes',quality_eligible='no' if degenerate or common!=planned else 'yes',n_degenerate_seeds=len(degenerate),clinical_interpretation_allowed='no'))
    return rows,dict(n_planned_jobs=len(planned),n_common_jobs=len(common),dropped_jobs=sorted(planned-common),degenerate_seeds_retained=degenerate,all_planned_seeds_retained=True,within_AF3_spread_not_included_in_std=True)

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores',type=Path,nargs='+',required=True)
    p.add_argument('--release',type=Path,required=True)
    p.add_argument('--seed-plan',type=Path,required=True,help='registered JSON score_files list; all seeds required')
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(argv)
    if a.out.exists() or Path(str(a.out)+'.provenance.json').exists():raise ValueError('refusing to overwrite ensemble')
    for path in a.scores:assert_not_withdrawn(path)
    from meta_acdc.data.release import load_release,sha256
    manifest,_=load_release(a.release)
    plan=json.loads(a.seed_plan.read_text())['score_files'];names=[p.name for p in a.scores]
    if len(set(plan))!=len(plan) or sorted(plan)!=sorted(names):raise ValueError('all predeclared seed files must match exactly')
    rh=sha256(a.release);receipts=[validate_score_receipt(path,rh) for path in a.scores];ids=[x['seed_ids'][0] for x in receipts]
    if len(set(ids))!=len(ids):raise ValueError('duplicate seed identity')
    if len({tuple(sorted(x['planned_jobs'])) for x in receipts})!=1:raise ValueError('different predeclared structure panels across seeds')
    rows,diagnostic=combine([load_seed(path) for path in a.scores],receipts[0]['planned_jobs'])
    if not rows:raise ValueError('no shared QC-passing jobs; no ensemble written')
    diagnostic['all_planned_models_scored']=all(x.get('all_planned_models_scored') is True for x in receipts)
    if not diagnostic['all_planned_models_scored']:
        for row in rows:row['quality_eligible']='no'
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    receipt=dict(scope='structural synthetic-control score diagnostic only; no binding or clinical validation',release_sha256=rh,output_sha256=sha256(a.out),seed_plan_sha256=sha256(a.seed_plan),seed_ids=ids,std_axis='between_seed_means',inputs={str(p):sha256(p) for p in a.scores},input_receipts={str(p)+'.provenance.json':sha256(Path(str(p)+'.provenance.json')) for p in a.scores},diagnostic=diagnostic)
    with Path(str(a.out)+'.provenance.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    return 0
if __name__=='__main__':raise SystemExit(main())
