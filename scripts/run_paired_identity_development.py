"""Frozen corrected-label DEVELOPMENT ablations; never accesses acceptance sets."""
import argparse,hashlib,json,random,fcntl,copy
from pathlib import Path
from collections import Counter
import numpy as np
import torch
from meta_acdc.models.fingerprint_factorial import train_model,encode,auc_rows,tensor_digest
from meta_acdc.models.train_fingerprint_v2 import EmbedderV2,set_seed
from scripts.fingerprint_shortcut_audit import cosine_rows
from scripts.rashomon_fingerprint import bootstrap_ci,finite_json
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def validate(path):
    if sha(path)!=path.with_suffix('.sha256').read_text().split()[0]:raise ValueError('protocol changed')
    p=json.loads(path.read_text())
    for name,h in p['frozen_sha256'].items():
        if sha(ROOT/name)!=h:raise ValueError('frozen file changed: '+name)
    p['protocol_sha256']=sha(path);return p

def training_records(records,fit,seed,arm):
    rows=copy.deepcopy(records)
    if arm.endswith('_shuffled'):
        ids=[records[i]['receptor_id'] for i in fit];random.Random(seed+9173).shuffle(ids)
        for i,label in zip(fit,ids):rows[i]['receptor_id']=label
        assert Counter(rows[i]['receptor_id'] for i in fit)==Counter(records[i]['receptor_id'] for i in fit)
    return rows

def baseline_matrices(data,fit):
    graphs=data['graphs'];records=data['records']
    composition=np.array([np.asarray(g[0])[:,4:24].mean(0) for g in graphs])
    size=np.log1p(np.array([[len(g[0]),len(g[2]),sum(r['peptide_mask'])] for g,r in zip(graphs,records)],float))
    mean=size[fit].mean(0);std=size[fit].std(0);std[std<1e-12]=1
    z=(size-mean)/std
    return {'composition':cosine_rows(composition),'graph_size':-np.sqrt(((z[:,None]-z[None,:])**2).sum(-1)),'constant':np.zeros((len(graphs),len(graphs)))},dict(mean=mean.tolist(),std=std.tolist(),fit_indices=fit)

def metric(S,data,fold):
    return {part:{'summary':summary,'query_auc':rows.tolist()} for part,idx in [('fit',fold['fit_indices']),('validation',fold['validation_indices'])] for summary,rows in [auc_rows(S,data['records'],idx)]}

def instance(p,data,fold,seed,arm,out,verify=False):
    marker=out/'metrics.json';basearm='reference' if arm=='reference' else 'identity_preserved';raw=basearm!='reference'
    if marker.exists():
        saved=json.loads(marker.read_text())
        if saved['protocol_sha256']!=p['protocol_sha256'] or any(sha(out/k)!=h for k,h in saved['artifact_sha256'].items()):raise ValueError('foreign or modified completed output')
        if not verify:return saved
        ck=torch.load(out/'model.pt',map_location='cpu',weights_only=False);model=EmbedderV2(26,16,emb_dim=p['recipe']['emb_dim']);model.load_state_dict(ck['model'])
        E=encode(model,data['graphs'],data['records'],list(range(len(data['graphs']))),raw).numpy();S=E@E.T
        if not np.array_equal(E,np.load(out/'embeddings.npy')) or not np.array_equal(S,np.load(out/'similarities.npy')):raise ValueError('weight to embedding/similarity mismatch')
        if finite_json(metric(S,data,fold))!=saved['metrics']:raise ValueError('weight to metrics mismatch')
        expected=training_records(data['records'],fold['fit_indices'],seed,arm)
        labels=json.loads((out/'training_labels.json').read_text())
        if labels!=[dict(index=i,receptor_id=expected[i]['receptor_id'],canonical_receptor_id=data['records'][i]['receptor_id']) for i in fold['fit_indices']]:raise ValueError('label control mismatch')
        if ck['scope']!='complete deposited paired-variable sequence diagnostic only':raise ValueError('promotion gate lost')
        return saved
    if verify:raise FileNotFoundError(marker)
    if out.exists() and any(out.iterdir()):raise ValueError('partial output cannot be overwritten')
    out.mkdir(parents=True,exist_ok=True);records=training_records(data['records'],fold['fit_indices'],seed,arm)
    model,initial,trace,coverage=train_model(data['graphs'],records,fold['fit_indices'],seed,basearm,p['recipe'])
    torch.save(dict(model=model.state_dict(),arm=arm,seed=seed,fold_id=fold['fold_id'],initial_tensor_sha256=initial,protocol_sha256=p['protocol_sha256'],scope='complete deposited paired-variable sequence diagnostic only'),out/'model.pt')
    E=encode(model,data['graphs'],data['records'],list(range(len(data['graphs']))),raw).numpy();S=E@E.T
    np.save(out/'embeddings.npy',E);np.save(out/'similarities.npy',S)
    (out/'training_labels.json').write_text(json.dumps([dict(index=i,receptor_id=records[i]['receptor_id'],canonical_receptor_id=data['records'][i]['receptor_id']) for i in fold['fit_indices']],indent=2)+'\n')
    (out/'trace.json').write_text(json.dumps(trace,indent=2)+'\n')
    saved=dict(arm=arm,fold_id=fold['fold_id'],seed=seed,protocol_sha256=p['protocol_sha256'],initial_tensor_sha256=initial,coverage=coverage,metrics=metric(S,data,fold),artifact_sha256={x.name:sha(x) for x in out.iterdir() if x.is_file()})
    marker.write_text(json.dumps(finite_json(saved),indent=2,allow_nan=False)+'\n');return saved

def report(p,data,results,out):
    eligible=sorted({i for f in p['folds'] for i in f['validation_indices']});groups=[data['records'][i]['family_id'] for i in eligible];baselines={};baseline_details={};by_arm={}
    index={i:j for j,i in enumerate(eligible)}
    for fold in p['folds']:
        matrices,norm=baseline_matrices(data,fold['fit_indices']);baseline_details[fold['fold_id']]={'size_training_normalization':norm}
        for name,S in matrices.items():
            summary,rows=auc_rows(S,data['records'],fold['validation_indices']);baseline_details[fold['fold_id']][name]=finite_json(dict(summary=summary,query_auc=rows.tolist()))
            baselines.setdefault(name,np.full(len(eligible),np.nan))
            for i,v in zip(fold['validation_indices'],rows):baselines[name][index[i]]=v
    for arm in p['arms']:
        a=np.full((len(eligible),len(p['seeds'])),np.nan)
        for fold in p['folds']:
            for j,seed in enumerate(p['seeds']):
                rows=results[fold['fold_id']][arm][str(seed)]['metrics']['validation']['query_auc']
                for i,v in zip(fold['validation_indices'],rows):a[index[i],j]=np.nan if v is None else v
        by_arm[arm]=a
    aggregate={a:bootstrap_ci(v,groups) for a,v in by_arm.items()};base={a:bootstrap_ci(v,groups) for a,v in baselines.items()}
    comparisons={a:{b:bootstrap_ci(v-baselines[b][:,None],groups) for b in baselines} for a,v in by_arm.items()}
    effects={'identity_preservation_minus_reference':bootstrap_ci(by_arm['identity_preserved']-by_arm['reference'],groups),'identity_preserved_minus_shuffled':bootstrap_ci(by_arm['identity_preserved']-by_arm['identity_preserved_shuffled'],groups)}
    value=dict(protocol_id=p['protocol_id'],protocol_sha256=p['protocol_sha256'],scope=p['scope'],aggregate=aggregate,baselines=base,paired_comparisons=comparisons,effects=effects,baseline_details=baseline_details,all_instances=results,eligible_validation_indices=eligible,n_independent_contexts=len(p['folds']),promotion=False,independent_acceptance_accessed=False,ci_caveat='descriptive conditional cluster bootstrap on only two V-gene contexts, not generalization confidence or confirmation')
    (out/'evaluation.json').write_text(json.dumps(finite_json(value),indent=2,allow_nan=False)+'\n')
    return value

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--protocol',type=Path,required=True);ap.add_argument('--out-dir',type=Path,required=True);ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=validate(args.protocol);data=torch.load(ROOT/p['dataset'],map_location='cpu',weights_only=False);out=args.out_dir;out.mkdir(parents=True,exist_ok=True);torch.set_num_threads(p['recipe']['threads'])
    with (out/'execution.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);results={};checks=[]
        for fold in p['folds']:
            fit,val=set(fold['fit_indices']),set(fold['validation_indices'])
            assert not fit&val and fit|val==set(range(len(data['records'])))
            for key in ['receptor_id','family_id','alpha_variable','beta_variable']:
                assert not {data['records'][i][key] for i in fit}&{data['records'][i][key] for i in val},key
            results[fold['fold_id']]={}
            for arm in p['arms']:
                results[fold['fold_id']][arm]={}
                for seed in p['seeds']:
                    validate(args.protocol);label=fold['fold_id']+'/'+arm+'/'+str(seed);print('VERIFY' if args.verify else 'INSTANCE',label,flush=True)
                    result=instance(p,data,fold,seed,arm,out/label,args.verify);results[fold['fold_id']][arm][str(seed)]=result;checks.append(label)
            for seed in p['seeds']:
                assert len({results[fold['fold_id']][a][str(seed)]['initial_tensor_sha256'] for a in p['arms']})==1
        first=p['folds'][0];seed=p['seeds'][0];replay=instance(p,data,first,seed,'reference',out/'replay',args.verify)
        original=out/first['fold_id']/'reference'/str(seed)
        aa=torch.load(original/'model.pt',map_location='cpu',weights_only=False)['model'];bb=torch.load(out/'replay/model.pt',map_location='cpu',weights_only=False)['model']
        assert all(torch.equal(aa[k],bb[k]) for k in aa) and np.array_equal(np.load(original/'embeddings.npy'),np.load(out/'replay/embeddings.npy'))
        value=report(p,data,results,out)
        if args.verify:(out.parent/'verification.json').write_text(json.dumps(dict(status='pass',n_instances=len(checks)+1,weight_embedding_similarity_metrics_exact=True,shuffle_mapping_exact=True,same_seed_initial_weights_equal=True,fit_validation_disjoint=True,replay_bitwise_equal=True,independent_acceptance_accessed=False),indent=2)+'\n')
        (out/'status.json').write_text(json.dumps(dict(state='verified' if args.verify else 'complete',n_instances=len(checks)+1),indent=2)+'\n');print('COMPLETE',finite_json(value['aggregate']),flush=True)
if __name__=='__main__':main()
