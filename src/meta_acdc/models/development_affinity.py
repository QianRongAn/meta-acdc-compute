"""Bounded development-only affinity comparison. Consumed holdouts are unreadable."""
from __future__ import annotations
import argparse,json,math,os,sys,hashlib
from pathlib import Path
from collections import Counter
from contextlib import contextmanager
import numpy as np
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import Ridge
from .affinity_baseline import features as legacy_features
from .frozen_affinity_baseline import weighted_metrics
from meta_acdc.structure.graph import AA_PHYSICOCHEM
from meta_acdc.data.release import sha256

ROOT=Path(__file__).resolve().parents[3]
AA='ACDEFGHIKLMNPQRSTVWY'
_GUARD=None

def _audit(event,args):
    if event!='open' or _GUARD is None or not isinstance(args[0],(str,bytes,os.PathLike)):return
    p=Path(os.fsdecode(args[0])).resolve();mode=args[1] or ''
    if any(x in mode for x in 'wax+'):return
    if not any(p.is_relative_to(ROOT/n) for n in ('docs','data')):return
    if p in _GUARD['allowed'] or p.is_relative_to(_GUARD['output']):_GUARD['reads'].add(str(p.relative_to(ROOT)));return
    raise PermissionError('development-only data firewall rejects '+str(p))
sys.addaudithook(_audit)

@contextmanager
def data_firewall(allowed,output):
    global _GUARD
    if _GUARD is not None:raise ValueError('nested data firewall')
    state=dict(allowed={Path(p).resolve() for p in allowed},output=Path(output).resolve(),reads=set());_GUARD=state
    try:yield state
    finally:_GUARD=None

def _properties(seq):
    return np.asarray([[AA_PHYSICOCHEM[a][0],AA_PHYSICOCHEM[a][1]/5,AA_PHYSICOCHEM[a][2]/200] for a in seq],dtype=float)

def chemical_features(row,words=False):
    out={};whole={}
    for key,seq in [('a',row['cdr3_alpha']),('b',row['cdr3_beta']),('p',row['peptide'])]:
        if not seq or set(seq)-set(AA):raise ValueError('invalid CDR3/peptide')
        regions={'all':seq,'left':seq[:max(1,len(seq)//3)],'middle':seq[len(seq)//3:max(len(seq)//3+1,2*len(seq)//3)],'right':seq[2*len(seq)//3:]}
        out[key+':length']=len(seq)/20;whole[key]=_properties(seq).mean(axis=0)
        for region,part in regions.items():
            p=_properties(part)
            for i,prop in enumerate(('charge','hydrophobicity','volume')):
                out[f'{key}:{region}:mean:{prop}']=float(p[:,i].mean());out[f'{key}:{region}:sd:{prop}']=float(p[:,i].std())
            for a in AA:out[f'{key}:{region}:AA:{a}']=part.count(a)/len(part)
        if words:
            c=Counter(seq[i:i+2] for i in range(len(seq)-1))
            for word,n in c.items():out[f'{key}:word2:{word}']=n/max(1,len(seq)-1)
    for key in ('a','b'):
        for i in range(3):
            for j in range(3):out[f'interaction:{key}:{i}:{j}']=float(whole[key][i]*whole['p'][j])
    return out

def feature_map(row,policy):
    if policy in ('cdr3_chemical','cdr3_words'):return chemical_features(row,words=policy=='cdr3_words')
    x=dict(row,method=row['primary_assay']['instrument']);f=legacy_features(x)
    if policy=='anchor_raw':return f
    if policy=='sequence_balanced':return {k:v for k,v in f.items() if not k.startswith(('mhc=','instrument=','temperature'))}
    raise ValueError('unknown feature policy')

def group(key):
    if key.startswith('interaction'):return 'interaction'
    if key.startswith(('a','b','p')):return key[0]
    return 'metadata'

def weights(rows,normalize=True):
    n=Counter(x['split_block'] for x in rows);w=np.array([1/n[x['split_block']] for x in rows])
    return w*len(rows)/w.sum() if normalize else w

def train(rows,policy,alpha):
    vec=DictVectorizer(sparse=False);X=vec.fit_transform([feature_map(x,policy) for x in rows]);w=weights(rows,normalize=policy!='anchor_raw');y=np.array([x['pkd'] for x in rows]);names=list(vec.feature_names_)
    center=np.zeros(X.shape[1]);scale=np.ones(X.shape[1]);traces={}
    if policy!='anchor_raw':
        center=np.average(X,axis=0,weights=w)
        for g in sorted({group(n) for n in names}):
            ids=[i for i,n in enumerate(names) if group(n)==g];trace=float(np.average((X[:,ids]-center[ids])**2,axis=0,weights=w).sum());traces[g]=trace;scale[ids]=0 if trace<=1e-12 else 1/math.sqrt(trace)
    model=Ridge(alpha=alpha,solver='svd').fit((X-center)*scale,y,sample_weight=w);coef=model.coef_*scale;intercept=float(model.intercept_-np.dot(coef,center))
    return dict(feature_policy=policy,alpha=alpha,feature_names=names,coefficients=coef.tolist(),intercept=intercept,n_fit=len(rows),fit_blocks=sorted({x['split_block'] for x in rows}),preprocessing=dict(weighted_center=center.tolist(),group_scale=scale.tolist(),group_variance_traces=traces),target='pKD',external_validated=False)

def predict(model,rows):
    if len(model['feature_names'])!=len(model['coefficients']):raise ValueError('model dimensions')
    result=[]
    for x in rows:
        f=feature_map(x,model['feature_policy']);result.append(model['intercept']+math.fsum(c*f.get(n,0.) for n,c in zip(model['feature_names'],model['coefficients'])))
    return result

def inner_cv(rows,policy,alphas):
    blocks=sorted({x['split_block'] for x in rows});results=[]
    for alpha in alphas:
        folds=[]
        for b in blocks:
            fit=[x for x in rows if x['split_block']!=b];held=[x for x in rows if x['split_block']==b];m=train(fit,policy,alpha);p=predict(m,held);metric=weighted_metrics([x['pkd'] for x in held],p,weights(held));folds.append(dict(held_block=b,fit_ids=[x['evidence_id'] for x in fit],held_ids=[x['evidence_id'] for x in held],metrics=metric))
        results.append(dict(alpha=alpha,block_equal_rmse=math.sqrt(sum(x['metrics']['rmse_pkd']**2 for x in folds)/len(folds)),block_equal_mae=sum(x['metrics']['mae_pkd'] for x in folds)/len(folds),folds=folds))
    return min(results,key=lambda x:(x['block_equal_rmse'],x['alpha']))['alpha'],results

def run(plan_path,output):
    plan_path=Path(plan_path).resolve();output=Path(output).resolve()
    if output.exists():raise ValueError('refusing duplicate development experiment')
    if not plan_path.is_relative_to(ROOT/'docs/affinity_development_20261018') or not plan_path.name.startswith('plan'):raise PermissionError('only registered development plan may be loaded')
    plan=json.loads(plan_path.read_text());dev_path=(ROOT/plan['development_path']).resolve()
    if sha256(dev_path)!=plan['development_sha256']:raise ValueError('development recipe changed')
    if plan['development_sha256']!=plan['parent_development_sha256'] and not plan.get('source_excluded_ids'):raise ValueError('unregistered development subset')
    output.mkdir(parents=True,exist_ok=False)
    with data_firewall([plan_path,dev_path],output) as access:
        rows=json.loads(dev_path.read_text());tr=[x for x in rows if x['fold']=='train'];va=[x for x in rows if x['fold']=='validation']
        if plan.get('source_excluded_ids'):
            expected=set(plan['parent_record_sha256'])-set(plan['source_excluded_ids'])
            if {x['evidence_id'] for x in rows}!=expected:raise ValueError('source-safe subset coverage')
            for x in rows:
                if hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()!=plan['parent_record_sha256'][x['evidence_id']]:raise ValueError('parent row changed')
        if len(tr)+len(va)!=len(rows) or len({x['evidence_id'] for x in rows})!=len(rows):raise ValueError('invalid development folds')
        if {x['split_block'] for x in tr}&{x['split_block'] for x in va}:raise ValueError('development block leakage')
        reports=[];models=[];predictions=[]
        for spec in plan['candidates']:
            policy=spec['name'];alpha,cv=(spec['alphas'][0],[]) if policy=='anchor_raw' else inner_cv(tr,policy,spec['alphas']);m=train(tr,policy,alpha);p=predict(m,va);score=weighted_metrics([x['pkd'] for x in va],p,weights(va));reports.append(dict(name=policy,alpha=alpha,training_cv=cv,validation=score,n_coefficients=len(m['coefficients']),validation_prediction_sd=float(np.std(p)),validation_prediction_range=[min(p),max(p)]));models.append(m);predictions.append(p)
        best=min(range(len(reports)),key=lambda i:(reports[i]['validation']['rmse_pkd'],reports[i]['n_coefficients'],reports[i]['name']));anchor=next(i for i,x in enumerate(reports) if x['name']=='anchor_raw');mean=float(np.average([x['pkd'] for x in tr],weights=weights(tr)));baseline=weighted_metrics([x['pkd'] for x in va],[mean]*len(va),weights(va));selected=models[best]
        result=dict(selected=reports[best]['name'],selected_alpha=reports[best]['alpha'],validation=reports[best]['validation'],anchor_validation=reports[anchor]['validation'],train_mean_validation=baseline,improves_both_validation_metrics_over_anchor=all(reports[best]['validation'][k]<reports[anchor]['validation'][k] for k in ('rmse_pkd','mae_pkd')),candidates=reports,training_rows=len(tr),validation_rows=len(va),test_used_for_selection=False,consumed_test_evaluated=False,new_test_evaluated=False,independent_acceptance_pending=True)
        for name,value in [('comparison.json',result),('selected_model.json',selected),('validation_predictions.json',[dict(evidence_id=x['evidence_id'],block=x['split_block'],measured_pkd=x['pkd'],selected_prediction=p,anchor_prediction=a,train_mean=mean) for x,p,a in zip(va,predictions[best],predictions[anchor])])]:
            (output/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
        result['data_reads']=sorted(access['reads']);(output/'data_access.json').write_text(json.dumps(dict(allowed_inputs=[str(plan_path.relative_to(ROOT)),str(dev_path.relative_to(ROOT))],observed_reads=sorted(access['reads']),consumed_test_accesses=0),indent=2)+'\n')
        code=['src/meta_acdc/models/development_affinity.py','src/meta_acdc/models/affinity_baseline.py','src/meta_acdc/models/frozen_affinity_baseline.py','src/meta_acdc/structure/graph.py'];lock=dict(plan_sha256=sha256(plan_path),development_sha256=sha256(dev_path),model_sha256=sha256(output/'selected_model.json'),selection='train-only LOBO for candidatealpha, original validation for candidate choice; no test',code_files={n:sha256(ROOT/n) for n in code},artifacts={p.name:sha256(p) for p in output.iterdir() if p.is_file()},test_evaluated=False);(output/'model_lock.json').write_text(json.dumps(lock,indent=2)+'\n')
    return {k:v for k,v in result.items() if k!='candidates'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args();print(json.dumps(run(a.plan,a.out),indent=2))
if __name__=='__main__':main()
