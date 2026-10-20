"""Independent process: validate lock, predict once, evaluate the frozen SPR panel.

No model fitting or hyperparameter search occurs in this entry point. Metrics are
recomputed with math (not the fitting process's numpy metrics implementation).
"""
from __future__ import annotations
import argparse,math,os,json
from pathlib import Path
from collections import Counter
from .affinity_baseline import features
from .frozen_affinity_baseline import verify_lock
from meta_acdc.data.affinity_training_release import ROOT,verify,read,write,now
from meta_acdc.data.release import sha256


def hierarchical_weights(rows):
    if not rows:raise ValueError('no evaluation observations')
    sources={x['research_source'] for x in rows}
    lineages={s:{x['lineage'] for x in rows if x['research_source']==s} for s in sources}
    n=Counter((x['research_source'],x['lineage']) for x in rows)
    return [1/(len(sources)*len(lineages[x['research_source']])*n[x['research_source'],x['lineage']]) for x in rows]


def metrics(y,p,w):
    if not y or len(y)!=len(p) or len(y)!=len(w) or not all(math.isfinite(v) for v in list(y)+list(p)+list(w)) or any(x<=0 for x in w):raise ValueError('invalid acceptance metric inputs')
    den=math.fsum(w)
    return dict(rmse_pkd=math.sqrt(math.fsum(a*(b-c)**2 for a,b,c in zip(w,y,p))/den),mae_pkd=math.fsum(a*abs(b-c) for a,b,c in zip(w,y,p))/den)


def locked_predict(model,rows):
    result=[]
    for row in rows:
        f=features(dict(row,method=row['instrument']))
        result.append(float(model['intercept'])+math.fsum(float(c)*f.get(n,0.) for n,c in zip(model['feature_names'],model['coefficients'])))
    if not all(math.isfinite(x) for x in result):raise ValueError('nonfinite locked prediction')
    return result


def claim_test(path,claim):
    # Exclusive create is the permanent access event, retained even after a crash.
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    try:write(path,claim)
    except FileExistsError:raise ValueError('frozen test already claimed; duplicate acceptance/refit forbidden') from None


def compare(rows,predictions,mean):
    w=hierarchical_weights(rows);y=[x['pkd'] for x in rows]
    model=metrics(y,predictions,w);constant=metrics(y,[mean]*len(rows),w)
    return dict(n=len(rows),n_sources=len({x['research_source'] for x in rows}),n_lineages=len({x['lineage'] for x in rows}),model=model,train_mean=constant,rmse_improvement=constant['rmse_pkd']-model['rmse_pkd'],mae_improvement=constant['mae_pkd']-model['mae_pkd'])


def acceptance_rule(primary,leave_one_source_out):
    return (primary['rmse_improvement']>0 and primary['mae_improvement']>0
            and all(x['rmse_improvement']>0 for x in leave_one_source_out.values()))


def accept(release,run,output):
    release,run,output=map(Path,(release,run,output))
    if output.exists():raise ValueError('refusing duplicate acceptance output')
    m=verify(release);lock=verify_lock(release,run);plan=read(release/'plan.json')
    if os.getpid()==lock['pid']:raise ValueError('acceptance must run in a separate process')
    # Reproduce the comparator and candidate choice using development only.
    dev=read(release/'development.json');tr=[x for x in dev if x['fold']=='train'];counts=Counter(x['split_block'] for x in tr);w=[1/counts[x['split_block']] for x in tr];mean=math.fsum(x['pkd']*a for x,a in zip(tr,w))/math.fsum(w)
    model=read(run/'model.json');fit_report=read(run/'fit_report.json');candidates=fit_report['candidate_selection'];best=min(candidates,key=lambda x:(x['validation']['rmse_pkd'],x['alpha']))
    if model['selected_alpha']!=best['alpha'] or [x['alpha'] for x in candidates]!=plan['alpha_grid'] or not math.isclose(mean,model['train_mean_pkd'],abs_tol=1e-12):raise ValueError('development selection/comparator inconsistent')
    output.mkdir(parents=True,exist_ok=False)
    claim=dict(claimed_utc=now(),acceptance_pid=os.getpid(),fit_pid=lock['pid'],model_lock_sha256=sha256(run/'model_lock.json'),model_sha256=sha256(run/'model.json'),parent_freeze_manifest_sha256=m['parent_manifest_sha256'],release_manifest_sha256=sha256(release/'manifest.json'),status='test access consumed; no further fitting on this profile',model_locked_before_access=True)
    claim_test(ROOT/m['test_access_ledger'],claim);write(output/'test_access_claim.json',claim)
    # First decode of held-out records in this entire model fitting/acceptance flow.
    test=read(ROOT/m['test_path']);pred=locked_predict(model,test);primary_ids=[i for i,x in enumerate(test) if x['metric_role']=='primary'];main=[test[i] for i in primary_ids];p=[pred[i] for i in primary_ids]
    primary=compare(main,p,mean);all15=compare(test,pred,mean);loo={};by_source={};by_lineage={};by_study={}
    for source in sorted({x['research_source'] for x in main}):
        ids=[i for i,x in enumerate(main) if x['research_source']!=source];loo[source]=compare([main[i] for i in ids],[p[i] for i in ids],mean)
        ids=[i for i,x in enumerate(main) if x['research_source']==source];by_source[source]=compare([main[i] for i in ids],[p[i] for i in ids],mean)
    for key,dest in [('lineage',by_lineage),('pmid',by_study)]:
        for name in sorted({x[key] for x in main}):
            ids=[i for i,x in enumerate(main) if x[key]==name];dest[name]=compare([main[i] for i in ids],[p[i] for i in ids],mean)
    orientation={'39576860':'pMHC immobilized','35484264':'TCR immobilized','33468649':'TCR immobilized','39042676':'pMHC immobilized','39500929':'pMHC immobilized'};by_orientation={}
    for label in sorted(set(orientation.values())):
        ids=[i for i,x in enumerate(main) if orientation[x['pmid']]==label];by_orientation[label]=compare([main[i] for i in ids],[p[i] for i in ids],mean)
    result=dict(created_utc=now(),acceptance_process_pid=os.getpid(),fit_process_pid=lock['pid'],independent_process=True,model_lock_verified=True,model_lock_sha256=claim['model_lock_sha256'],model_sha256=claim['model_sha256'],release_manifest_sha256=claim['release_manifest_sha256'],test_sha256=m['test_sha256'],plan_sha256=sha256(release/'plan.json'),selected_alpha=model['selected_alpha'],train_mean_pkd=mean,primary=primary,all15_sensitivity=all15,leave_one_source_out_metric_sensitivity=loo,by_source=by_source,by_lineage=by_lineage,by_study=by_study,by_assay_orientation=by_orientation,acceptance_predeclared=plan['acceptance_predeclared'],exploratory_performance_passed=acceptance_rule(primary,loo),general_affinity_validation_passed=False,binary_binding_validation_passed=False,clinical_validation_passed=False,structural_EGNN_validation_passed=False,confidence_interval=None,uncertainty=plan['uncertainty'],test_labels_used_for_selection=False,post_test_fit_allowed=False,scope='one frozen source-disjoint test of a10-observation-trained sequence/assay Ridge baseline; not population or clinical validation')
    weights=hierarchical_weights(main);wm={x['evidence_id']:a for x,a in zip(main,weights)}
    write(output/'test_predictions.json',[dict(evidence_id=x['evidence_id'],lineage=x['lineage'],pmid=x['pmid'],research_source=x['research_source'],metric_role=x['metric_role'],measured_pkd=x['pkd'],predicted_pkd=a,train_mean_pkd=mean,primary_weight=wm.get(x['evidence_id'],0),assay_orientation=orientation[x['pmid']]) for x,a in zip(test,pred)])
    write(output/'report.json',result);write(output/'receipt.json',dict(completed_utc=now(),test_access_claim_sha256=sha256(ROOT/m['test_access_ledger']),report_sha256=sha256(output/'report.json'),predictions_sha256=sha256(output/'test_predictions.json'),locked_model_unchanged=sha256(run/'model.json')==lock['model_sha256'],no_refit_performed=True))
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--release',required=True,type=Path);p.add_argument('--run',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args();print(json.dumps(accept(a.release,a.run,a.out),indent=2))
if __name__=='__main__':main()
