"""Frozen 2x2 developmental diagnostic. Independent acceptance data remain unopened."""
import argparse
from collections import defaultdict
import fcntl
import hashlib
import json
from pathlib import Path
import time
import traceback
import numpy as np
import torch
from meta_acdc.models.fingerprint_factorial import ARMS,train_model,encode,auc_rows,tensor_digest
from scripts.fingerprint_shortcut_audit import cosine_rows
from scripts.rashomon_fingerprint import bootstrap_ci,finite_json
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def validate(path):
    if sha(path)!=path.with_suffix('.sha256').read_text().split()[0]:raise ValueError('protocol hash mismatch')
    p=json.loads(path.read_text())
    for rel,expected in p['frozen_sha256'].items():
        if sha(ROOT/rel)!=expected:raise ValueError('frozen input/source changed: '+rel)
    p["protocol_sha256"]=sha(path)
    return p

def save_instance(p,data,arm,seed,out):
    marker=out/'metrics.json'
    if marker.exists():
        saved=json.loads(marker.read_text())
        if saved.get('protocol_sha256')!=p['protocol_sha256'] or saved['protocol_id']!=p['protocol_id'] or any(sha(out/k)!=v for k,v in saved['artifact_sha256'].items()):
            raise ValueError('foreign or modified completed instance')
        return saved
    if out.exists() and any(out.iterdir()):raise ValueError('partial instance cannot be overwritten')
    out.mkdir(parents=True,exist_ok=True)
    model,initial,trace,coverage=train_model(data['graphs'],data['records'],p['fit_indices'],seed,arm,p['recipe'])
    torch.save(dict(model=model.state_dict(),arm=arm,seed=seed,protocol_id=p['protocol_id'],
                    diagnostic_only=True,identity_truth_validated=False,protocol_sha256=p['protocol_sha256'],initial_tensor_sha256=initial),out/'model.pt')
    E=encode(model,data['graphs'],data['records'],list(range(len(data['graphs']))),ARMS[arm]['raw_identity'])
    S=(E@E.T).numpy();np.save(out/'embeddings.npy',E.numpy());np.save(out/'similarities.npy',S)
    (out/'trace.json').write_text(json.dumps(trace,indent=2)+'\n')
    result=dict(protocol_id=p['protocol_id'],protocol_sha256=p['protocol_sha256'],arm=arm,seed=seed,initial_tensor_sha256=initial,coverage=coverage,partitions={},interventions={})
    for partition,indices in [('fit',p['fit_indices']),('validation',p['validation_indices'])]:
        result['partitions'][partition]={}
        for name,mask,hard in [('raw_hard',[],True),('known_same_excluded',p['known_same_pairs'],True),('all',[],False)]:
            summary,rows=auc_rows(S,data['records'],indices,mask,hard)
            result['partitions'][partition][name]=dict(summary=summary,query_auc=rows.tolist())
    for intervention in p['interventions']:
        indices=p['validation_indices'];e=encode(model,data['graphs'],data['records'],indices,ARMS[arm]['raw_identity'],intervention)
        np.save(out/(intervention+'.npy'),e.numpy())
        r=[data['records'][i] for i in indices];summary,rows=auc_rows((e@e.T).numpy(),r,list(range(len(r))))
        result['interventions'][intervention]=dict(summary=summary,query_auc=rows.tolist())
    result['artifact_sha256']={x.name:sha(x) for x in out.iterdir() if x.is_file()}
    marker.write_text(json.dumps(finite_json(result),indent=2,allow_nan=False)+'\n')
    return result

def report(p,data,results,out):
    records=data['records'];baseline={};base_rows={}
    h=np.array([np.asarray(g[0])[:,4:24].mean(0) for g in data['graphs']]);S=cosine_rows(h)
    for partition,idx in [('fit',p['fit_indices']),('validation',p['validation_indices'])]:
        baseline[partition]={};base_rows[partition]={}
        for name,mask in [('raw_hard',[]),('known_same_excluded',p['known_same_pairs'])]:
            summary,rows=auc_rows(S,records,idx,mask);baseline[partition][name]=summary;base_rows[partition][name]=rows
    aggregate={};rows_by_arm={};effects={}
    for arm in ARMS:
        aggregate[arm]={};rows_by_arm[arm]={}
        for part,idx in [('fit',p['fit_indices']),('validation',p['validation_indices'])]:
            aggregate[arm][part]={};rows_by_arm[arm][part]={}
            groups=[records[i]['family_id'] for i in idx]
            for metric in ['raw_hard','known_same_excluded','all']:
                a=np.array([results[arm][s]['partitions'][part][metric]['query_auc'] for s in p['seeds']],dtype=float).T
                rows_by_arm[arm][part][metric]=a
                aggregate[arm][part][metric]=bootstrap_ci(a,groups)
    contrasts={'identity_at_random':('identity_preserved','reference'),
               'identity_at_structured':('identity_and_family','family_batches'),
               'batching_at_legacy':('family_batches','reference'),
               'batching_at_raw_identity':('identity_and_family','identity_preserved')}
    for part,idx in [('fit',p['fit_indices']),('validation',p['validation_indices'])]:
        groups=[records[i]['family_id'] for i in idx];effects[part]={}
        for metric in ['raw_hard','known_same_excluded']:
            effects[part][metric]={}
            for name,(a,b) in contrasts.items():
                effects[part][metric][name]=bootstrap_ci(rows_by_arm[a][part][metric]-rows_by_arm[b][part][metric],groups)
            interaction=(rows_by_arm['identity_and_family'][part][metric]-rows_by_arm['family_batches'][part][metric]
                         -rows_by_arm['identity_preserved'][part][metric]+rows_by_arm['reference'][part][metric])
            effects[part][metric]['interaction']=bootstrap_ci(interaction,groups)
    changes={}
    groups=[records[i]['family_id'] for i in p['validation_indices']]
    for arm in ARMS:
        changes[arm]={}
        for intervention in p['interventions']:
            a=np.array([results[arm][s]['interventions'][intervention]['query_auc'] for s in p['seeds']],dtype=float).T
            changes[arm][intervention]=bootstrap_ci(a-rows_by_arm[arm]['validation']['raw_hard'],groups)
    seed_checks={str(s):len({results[a][s]['initial_tensor_sha256'] for a in ARMS})==1 for s in p['seeds']}
    if not all(seed_checks.values()):raise RuntimeError('arms did not share initial weights')
    result=dict(protocol_id=p['protocol_id'],scope=p['scope'],biological_receptor_candidate_promotable=False,
                promotion_blocker='observed-chain labels split identical deposited constructs; paired complete variable identity unverified',
                same_seed_initial_weights_equal=seed_checks,baselines=baseline,aggregate=aggregate,
                factorial_effects=effects,input_intervention_effects=changes,
                intervention_caveat='post-training feature corruption, not retraining ablation; topology retained in numeric-geometry mask; not used for selection',
                all_instance_metrics=results,independent_acceptance_sets_read=False)
    (out/'evaluation.json').write_text(json.dumps(finite_json(result),indent=2,allow_nan=False)+'\n')
    lines=['# Frozen 2x2 diagnostic results','',p['scope'],'',
           'NO biological receptor-model promotion: original labels are observed-sequence hashes with confirmed false negatives.',
           '', '| Arm | Fit hard AUROC | Validation hard AUROC | Validation known-same-excluded AUROC |',
           '|---|---:|---:|---:|']
    fmt=lambda x:'undefined' if x is None else f'{x:.4f}'
    for arm,a in aggregate.items():lines.append(f"| {arm} | {fmt(a['fit']['raw_hard']['mean'])} | {fmt(a['validation']['raw_hard']['mean'])} | {fmt(a['validation']['known_same_excluded']['mean'])} |")
    lines+=['',f"Composition validation hard AUROC: {baseline['validation']['raw_hard']['mean']}",
            f"Family-equal conditional intervals, all seeds retained; eligible validation families: {aggregate['reference']['validation']['raw_hard']['n_defined_family_groups']}.",
            'Known-same exclusion does not establish remaining negatives as true biological differences. No canonical identity labels were imputed.']
    (out/'README.md').write_text('\n'.join(lines)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--protocol',type=Path,required=True);ap.add_argument('--out-dir',type=Path,required=True);a=ap.parse_args()
    p=validate(a.protocol);out=a.out_dir;out.mkdir(parents=True,exist_ok=True)
    with (out/'execution.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        data=torch.load(ROOT/p['dataset'],map_location='cpu',weights_only=False)
        fit=set(p['fit_indices']);val=set(p['validation_indices'])
        if fit&val or fit|val!=set(range(len(data['graphs']))):raise ValueError('invalid train/validation partition')
        if {data['records'][i]['family_id'] for i in fit}&{data['records'][i]['family_id'] for i in val}:raise ValueError('family leakage')
        if set(data['source_indices'])&set(p['quarantined_source_indices']):raise ValueError('quarantined graphs entered development data')
        results=defaultdict(dict)
        def status(state,**extra):(out/'status.json').write_text(json.dumps(dict(state=state,protocol_id=p['protocol_id'],**extra),indent=2)+'\n')
        try:
            for seed in p['seeds']:
                for arm in ARMS:
                    validate(a.protocol);status('running',seed=seed,arm=arm)
                    results[arm][seed]=save_instance(p,data,arm,seed,out/arm/str(seed))
            status('replaying')
            replay=save_instance(p,data,'reference',p['seeds'][0],out/'replay')
            original=out/'reference'/str(p['seeds'][0]);replaydir=out/'replay'
            c=torch.load(original/'model.pt',weights_only=False)['model'];d=torch.load(replaydir/'model.pt',weights_only=False)['model']
            if not all(torch.equal(c[k],d[k]) for k in c) or not np.array_equal(np.load(original/'embeddings.npy'),np.load(replaydir/'embeddings.npy')):
                raise RuntimeError('determinism replay failed')
            (out/'replay_check.json').write_text(json.dumps(dict(weights_bitwise_equal=True,embeddings_bitwise_equal=True))+'\n')
            validate(a.protocol);report(p,data,results,out);status('complete',n_main_models=len(ARMS)*len(p['seeds']),n_replays=1)
        except BaseException as error:
            status('failed',error=repr(error));traceback.print_exc();raise

if __name__=='__main__':main()
