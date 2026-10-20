"""Fit on development only; publish a model lock before independent acceptance."""
from __future__ import annotations
import argparse,math,os,sys,platform
from pathlib import Path
from collections import Counter
import numpy as np
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import Ridge
from importlib.metadata import version
from .affinity_baseline import features
from meta_acdc.data.affinity_training_release import ROOT,verify,read,write,now
from meta_acdc.data.release import sha256


def model_features(row):
    x=dict(row);x['method']=row['primary_assay']['instrument'] if 'primary_assay' in row else row['instrument']
    return features(x)


def block_weights(rows):
    c=Counter(x['split_block'] for x in rows)
    return np.array([1/c[x['split_block']] for x in rows],dtype=float)


def weighted_metrics(truth,prediction,weights):
    y=np.asarray(truth,dtype=float);p=np.asarray(prediction,dtype=float);w=np.asarray(weights,dtype=float)
    if len(y)==0 or y.shape!=p.shape or y.shape!=w.shape or not np.isfinite(y).all() or not np.isfinite(p).all() or not np.isfinite(w).all() or (w<=0).any():raise ValueError('invalid metric input')
    return dict(rmse_pkd=float(np.sqrt(np.average((y-p)**2,weights=w))),mae_pkd=float(np.average(np.abs(y-p),weights=w)))


def predict_export(model,rows):
    names=model['feature_names'];coef=model['coefficients']
    if len(names)!=len(coef) or len(set(names))!=len(names) or not all(math.isfinite(x) for x in coef+[model['intercept']]):raise ValueError('invalid locked linear model')
    result=[]
    for row in rows:
        f=model_features(row);result.append(model['intercept']+math.fsum(c*f.get(n,0.) for n,c in zip(names,coef)))
    return result


def fit(release,output):
    release=Path(release);output=Path(output)
    if output.exists():raise ValueError('refusing duplicate training run')
    m=verify(release,for_fit=True);plan=read(release/'plan.json')
    output.mkdir(parents=True,exist_ok=False)
    fit_claim=ROOT/m['fit_access_ledger'];fit_claim.parent.mkdir(parents=True,exist_ok=True)
    try:write(fit_claim,dict(claimed_utc=now(),pid=os.getpid(),release_manifest_sha256=sha256(release/'manifest.json'),output=str(output)))
    except FileExistsError:raise ValueError('development fit already claimed; duplicate dispatch forbidden') from None
    write(output/'fit_started.json',dict(started_utc=now(),pid=os.getpid(),test_labels_loaded=False,release_manifest_sha256=sha256(release/'manifest.json')))
    dev=read(release/'development.json');tr=[x for x in dev if x['fold']=='train'];va=[x for x in dev if x['fold']=='validation']
    vec=DictVectorizer(sparse=False);xt=vec.fit_transform([model_features(x) for x in tr]);xv=vec.transform([model_features(x) for x in va])
    yt=np.array([x['pkd'] for x in tr]);yv=np.array([x['pkd'] for x in va]);wt=block_weights(tr);wv=block_weights(va);mean=float(np.average(yt,weights=wt));candidates=[];models=[]
    for alpha in plan['alpha_grid']:
        model=Ridge(alpha=alpha,solver='svd');model.fit(xt,yt,sample_weight=wt);score=weighted_metrics(yv,model.predict(xv),wv)
        candidates.append(dict(alpha=alpha,validation=score));models.append(model)
    k=min(range(len(candidates)),key=lambda i:(candidates[i]['validation']['rmse_pkd'],candidates[i]['alpha']));model=models[k]
    exported=dict(task=m['task'],feature_names=list(vec.feature_names_),coefficients=model.coef_.tolist(),intercept=float(model.intercept_),selected_alpha=candidates[k]['alpha'],train_mean_pkd=mean,fit_samples=len(tr),fit_connected_blocks=len({x['split_block'] for x in tr}),release_manifest_sha256=sha256(release/'manifest.json'))
    pred=predict_export(exported,va)
    if not np.allclose(pred,model.predict(xv),rtol=1e-11,atol=1e-11):raise ValueError('export does not reproduce sklearn')
    write(output/'model.json',exported)
    report=dict(created_utc=now(),n_train=len(tr),n_validation=len(va),candidate_selection=candidates,selected_alpha=exported['selected_alpha'],selected_validation=weighted_metrics(yv,pred,wv),validation_train_mean=weighted_metrics(yv,[mean]*len(va),wv),training=weighted_metrics(yt,model.predict(xt),wt),selection='validation only; train-only vectorizer and fit; no validation refit',test_labels_loaded=False,test_predictions_created=False,uses_pretrained_checkpoint=False,pid=os.getpid(),runtime=dict(python=sys.version,platform=platform.platform(),packages={n:version(n) for n in ['numpy','scipy','scikit-learn']}))
    write(output/'fit_report.json',report);write(output/'validation_predictions.json',[dict(evidence_id=x['evidence_id'],split_block=x['split_block'],measured_pkd=x['pkd'],predicted_pkd=p,train_mean_pkd=mean) for x,p in zip(va,pred)])
    lock=dict(locked_utc=now(),pid=os.getpid(),fit_access_claim_sha256=sha256(fit_claim),model_sha256=sha256(output/'model.json'),release_manifest_sha256=sha256(release/'manifest.json'),plan_sha256=sha256(release/'plan.json'),test_labels_loaded=False,files={p.name:sha256(p) for p in output.iterdir() if p.is_file()},code_files=m['code_files'],no_further_fitting=True)
    write(output/'model_lock.json',lock);return report


def verify_lock(release,run):
    m=verify(release);run=Path(run);lock=read(run/'model_lock.json')
    if lock['release_manifest_sha256']!=sha256(Path(release)/'manifest.json') or lock['plan_sha256']!=sha256(Path(release)/'plan.json') or lock['code_files']!=m['code_files'] or lock['test_labels_loaded'] is not False:raise ValueError('model lock provenance mismatch')
    if sha256(ROOT/m['fit_access_ledger'])!=lock['fit_access_claim_sha256']:raise ValueError('training claim changed')
    for n,h in lock['files'].items():
        p=(run/n).resolve()
        if not p.is_relative_to(run.resolve()) or sha256(p)!=h:raise ValueError('locked artifact changed '+n)
    if lock['model_sha256']!=sha256(run/'model.json'):raise ValueError('model changed after lock')
    dev=read(Path(release)/'development.json');va=[x for x in dev if x['fold']=='validation'];model=read(run/'model.json');report=read(run/'fit_report.json')
    actual=weighted_metrics([x['pkd'] for x in va],predict_export(model,va),block_weights(va))
    if any(not math.isclose(actual[k],report['selected_validation'][k],abs_tol=1e-10,rel_tol=1e-10) for k in actual):raise ValueError('locked validation reproduction failed')
    return lock


def main():
    p=argparse.ArgumentParser();p.add_argument('--release',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args();print(__import__('json').dumps(fit(a.release,a.out),indent=2))
if __name__=='__main__':main()
