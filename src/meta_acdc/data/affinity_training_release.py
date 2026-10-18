"""Versioned development-only release for the independent SPR baseline."""
from __future__ import annotations
import argparse,json,math,sys,platform,os
from pathlib import Path
from datetime import datetime,timezone
from .release import sha256
from .independent_spr import ROOT,FREEZE,verify_freeze

REPORT=ROOT/'docs/affinity_baseline_20261018'
ACCESS=REPORT/'test_access_ledger.json'
CODE=['src/meta_acdc/data/affinity_training_release.py','src/meta_acdc/models/frozen_affinity_baseline.py','src/meta_acdc/models/affinity_acceptance.py','src/meta_acdc/models/affinity_baseline.py','src/meta_acdc/data/independent_spr.py','src/meta_acdc/data/experimental_evidence.py','src/meta_acdc/data/build_release.py','src/meta_acdc/data/evidence_completion.py','src/meta_acdc/data/release.py']

def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def write(p,value):
    with Path(p).open('x') as f:json.dump(value,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

def build(output):
    output=Path(output)
    if output.exists():raise ValueError('refusing release overwrite')
    if ACCESS.exists():raise ValueError('frozen test already accessed; new training protocol required')
    gate=verify_freeze()
    dev=read(FREEZE/'development_records.json')
    for x in dev:
        if x['fold'] not in ('train','validation'):raise ValueError('test in development release')
        if x['measurement']['status']!='exact' or x['training_label'] is not None:raise ValueError('wrong objective')
        if not math.isclose(x['pkd'],6-math.log10(x['measurement']['kd_um']),rel_tol=1e-12):raise ValueError('pKD target mismatch')
    plan=dict(task='experimental_SPR_pKD_regression',model='Ridge; existing sequence/interaction/assay feature map',alpha_grid=[1.,10.,100.],selection='validation connected-block-equal weighted RMSE; lower alpha breaks exact ties',training_weights='equal total weight per inherited training connected block',validation_weights='equal total weight per inherited validation connected block',feature_fit='DictVectorizer fitted on train only; validation/test unseen categories ignored; no other normalization',assay_feature='instrument from verified primary_assay for development; frozen instrument for test; missing or unassigned instrument remains explicit',refit_on_validation=False,baseline='training connected-block-weighted mean pKD',test_weights='equal source total; equal lineage within source; equal observation within lineage',primary='14 metric_role=primary observations; original15 preserved; excluded extrapolated S8WT/A11 sensitivity only',sensitivity='all15; each source omitted in turn for metric sensitivity only, never refit',acceptance_predeclared='primary source/lineage-weighted RMSE and MAE both strictly lower than train-mean comparator AND RMSE lower in all3 leave-one-source-out metric subsets; no significance or clinical claim',uncertainty='3 source components insufficient for precise population CI; no bootstrap CI or p-value asserted',clinical_validation_allowed=False,binary_validation_allowed=False,structural_EGNN_validation_allowed=False,test_access_policy='exclusive persisted ledger before reading test labels; no new fits on this profile after test access; no retry/overwrite of completed output',created_utc=now())
    output.mkdir(parents=True,exist_ok=False)
    write(output/'development.json',dev);write(output/'plan.json',plan);write(output/'preflight.json',gate)
    manifest=dict(schema_version=1,task=plan['task'],training_allowed=True,scope='research-only affinity baseline; not structural, binary or clinical training',created_utc=now(),parent_freeze=str(FREEZE.relative_to(ROOT)),parent_manifest_sha256=sha256(FREEZE/'manifest.json'),test_path=str((FREEZE/'frozen_test.json').relative_to(ROOT)),test_sha256=sha256(FREEZE/'frozen_test.json'),test_access_ledger=str(ACCESS.relative_to(ROOT)),fit_access_ledger=str((REPORT/'fit_access_ledger.json').relative_to(ROOT)),code_files={n:sha256(ROOT/n) for n in CODE},files={p.name:sha256(p) for p in output.iterdir() if p.is_file()},counts={k:v for k,v in gate.items() if k not in ('training_allowed','model_scored','optimizer_started')},approval_basis='new user-authorized development-only profile; parent inspection state is immutable')
    write(output/'manifest.json',manifest);return manifest

def verify(release,*,for_fit=False):
    release=Path(release);m=read(release/'manifest.json')
    if m.get('task')!='experimental_SPR_pKD_regression' or m.get('training_allowed') is not True:raise ValueError('not released for affinity baseline')
    if for_fit and (ROOT/m['test_access_ledger']).exists():raise ValueError('frozen test already accessed; no additional fitting on this profile')
    if for_fit and (ROOT/m['fit_access_ledger']).exists():raise ValueError('development fit already claimed; duplicate dispatch forbidden')
    for n,h in m['files'].items():
        p=(release/n).resolve()
        if not p.is_relative_to(release.resolve()) or sha256(p)!=h:raise ValueError('release artifact changed '+n)
    for n,h in m['code_files'].items():
        p=(ROOT/n).resolve()
        if not p.is_relative_to(ROOT) or sha256(p)!=h:raise ValueError('release code changed '+n)
    freeze=ROOT/m['parent_freeze']
    if sha256(freeze/'manifest.json')!=m['parent_manifest_sha256']:raise ValueError('parent freeze changed')
    # Integrity reads only. No test JSON is decoded in the fitting path.
    if sha256(ROOT/m['test_path'])!=m['test_sha256']:raise ValueError('holdout changed')
    f=read(freeze/'manifest.json')
    for group in ('artifacts','code_files','previous_freezes','reference_files'):
        for n,h in f[group].items():
            if sha256(ROOT/n)!=h:raise ValueError('parent contract changed '+n)
    for n,h in read(freeze/'required_sources.json').items():
        if sha256(ROOT/n)!=h:raise ValueError('primary evidence changed '+n)
    dev=read(release/'development.json')
    if dev!=read(freeze/'development_records.json'):raise ValueError('development recipe changed')
    if any(x['fold'] not in ('train','validation') for x in dev):raise ValueError('test injected into development')
    ids={f:{x['split_block'] for x in dev if x['fold']==f} for f in ('train','validation')}
    if ids['train']&ids['validation'] or min(map(len,ids.values()))<3:raise ValueError('development blocks')
    return m

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args();print(json.dumps(build(a.out),indent=2))
if __name__=='__main__':main()
