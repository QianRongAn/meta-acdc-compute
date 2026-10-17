"""Fresh sequence/assay regression baseline for a frozen experimental SPR dataset.

No EGNN checkpoint, synthetic label, assay-classification AUROC or clinical claim.
Hyperparameters use validation only; test is evaluated once per fresh output.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter
from pathlib import Path
import numpy as np
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error,mean_squared_error
from meta_acdc.data.release import sha256
from meta_acdc.data.experimental_evidence import AA


def features(row):
    result={'mhc='+row['mhc_allele']:1.,'instrument='+row['method']:1.}
    for key,prefix in [('alpha_variable','a'),('beta_variable','b'),('peptide','p')]:
        seq=row[key];result[prefix+'_length']=len(seq)/100
        for width in [1,2]:
            counts=Counter(seq[i:i+width] for i in range(len(seq)-width+1));denom=max(1,len(seq)-width+1)
            for term,n in counts.items():result[prefix+str(width)+':'+term]=n/denom
    for aa in AA:
        for bb in AA:
            result['interaction:'+aa+bb]=((row['alpha_variable']+row['beta_variable']).count(aa)/len(row['alpha_variable']+row['beta_variable']))*(row['peptide'].count(bb)/len(row['peptide']))
    try:temperature=float(row['temperature_K'])
    except (ValueError,TypeError):temperature=None
    if temperature is not None and not math.isfinite(temperature):temperature=None
    result['temperature_known']=float(temperature is not None)
    if temperature is not None:result['temperature_centered']=(temperature-298.15)/10
    return result


def metrics(y,p):return {'rmse_pkd':float(np.sqrt(mean_squared_error(y,p))),'mae_pkd':float(mean_absolute_error(y,p))}


def run(release,output):
    if output.exists():raise ValueError('refusing duplicate run/output overwrite')
    from meta_acdc.data.atlas_binding import verify_release
    m=verify_release(release)
    rows=json.loads((release/'affinity_dataset.json').read_text());splits=json.loads((release/'splits.json').read_text());tr,va,te=[splits[f] for f in ['train','validation','test']]
    if set(tr)&set(va) or set(tr)&set(te) or set(va)&set(te):raise ValueError('overlapping folds')
    if set(tr+va+te)!=set(range(len(rows))):raise ValueError('fold coverage')
    if any(x['training_label'] is not None or x['measurement']['status']!='exact' for x in rows):raise ValueError('invalid objective labels')
    groups=[{rows[i]['split_block'] for i in ids} for ids in [tr,va,te]]
    if any(groups[i]&groups[j] for i in range(3) for j in range(i)):raise ValueError('block leakage')
    vec=DictVectorizer(sparse=False);xt=vec.fit_transform([features(rows[i]) for i in tr]);xv=vec.transform([features(rows[i]) for i in va]);xe=vec.transform([features(rows[i]) for i in te]);y=np.array([x['pkd'] for x in rows]);candidates=[]
    weights=Counter(rows[i]['pair_id'] for i in tr);sample_weight=np.array([1/weights[rows[i]['pair_id']] for i in tr])
    for alpha in [1.,10.,100.]:
        model=Ridge(alpha=alpha);model.fit(xt,y[tr],sample_weight=sample_weight);candidates.append((metrics(y[va],model.predict(xv))['rmse_pkd'],alpha,model))
    _,alpha,model=min(candidates,key=lambda x:(x[0],x[1]));prediction=model.predict(xe);constant=np.full(len(te),np.average(y[tr],weights=sample_weight));ridge=metrics(y[te],prediction);baseline=metrics(y[te],constant)
    # Bootstrap connected test blocks rather than treating within-family points as independent.
    block_ids=sorted(groups[2]);positions={b:[j for j,i in enumerate(te) if rows[i]['split_block']==b] for b in block_ids};rng=np.random.default_rng(20261017);improvements=[]
    for _ in range(2000):
        sampled=rng.choice(block_ids,len(block_ids),replace=True);idx=[j for b in sampled for j in positions[b]];truth=y[te][idx]
        improvements.append(metrics(truth,constant[idx])['rmse_pkd']-metrics(truth,prediction[idx])['rmse_pkd'])
    ci=list(map(float,np.quantile(improvements,[.025,.975])));passed=ci[0]>0 and ridge['mae_pkd']<baseline['mae_pkd']
    output.mkdir(parents=True,exist_ok=False)
    exported={'feature_names':list(vec.feature_names_),'coefficients':model.coef_.tolist(),'intercept':float(model.intercept_),'selected_alpha':alpha,'release_manifest_sha256':sha256(release/'manifest.json'),'task':m['task'],'fresh_fit':True}
    (output/'model.json').write_text(json.dumps(exported,indent=2)+'\n')
    report={'model_sha256':sha256(output/'model.json'),'task':m['task'],'release_manifest_sha256':sha256(release/'manifest.json'),'objective':'pKD = -log10(KD in molar); not probability of binding','n_train':len(tr),'n_validation':len(va),'n_test':len(te),'n_test_connected_blocks':len(block_ids),'hyperparameter_selection':'Ridge alpha 1,10,100; validation RMSE only; no refit on test','selected_alpha':alpha,'test_ridge':ridge,'test_train_mean_baseline':baseline,'rmse_improvement_block_bootstrap_95ci':ci,'acceptance_predeclared':'lower 95% block-bootstrap RMSE improvement > 0; lower MAE; exploratory comparison only, not a clinical/industry standard or broad affinity validation','exploratory_comparison_passed':passed,'general_affinity_validation_passed':False,'binary_binding_validation_passed':False,'clinical_validation_passed':False,'uses_pretrained_checkpoint':False,'external_primary_panel_used':False,'scope':'fresh sequence/assay exploratory model; does not validate structural EGNN or unseen nonbinders'}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');(output/'test_predictions.json').write_text(json.dumps([{'evidence_id':rows[i]['evidence_id'],'split_block':rows[i]['split_block'],'measured_pkd':float(y[i]),'predicted_pkd':float(p),'train_mean_pkd':float(c)} for i,p,c in zip(te,prediction,constant)],indent=2)+'\n');return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();print(json.dumps(run(a.release,a.out),indent=2))
if __name__=='__main__':main()
