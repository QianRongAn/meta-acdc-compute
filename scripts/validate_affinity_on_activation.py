"""One-shot external MIRA ranking check for a freshly fitted SPR affinity model.

Cross-endpoint association only: pKD predictions are not activation probabilities.
The external panel must be reserved from the affinity training release first.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from meta_acdc.data.release import sha256
from meta_acdc.data.atlas_binding import verify_release
from meta_acdc.data.build_release import domain
from meta_acdc.models.affinity_baseline import features

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--affinity-release',type=Path,required=True);p.add_argument('--model-run',type=Path,required=True);p.add_argument('--functional-release',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise ValueError('refusing duplicate external evaluation')
    m=verify_release(a.affinity_release);fit=json.loads((a.model_run/'report.json').read_text());model=json.loads((a.model_run/'model.json').read_text())
    if fit['model_sha256']!=sha256(a.model_run/'model.json') or model['release_manifest_sha256']!=sha256(a.affinity_release/'manifest.json'):raise ValueError('model/release mismatch')
    fm=json.loads((a.functional_release/'manifest.json').read_text())
    for name,digest in fm['files'].items():
        if sha256(a.functional_release/name)!=digest:raise ValueError('functional artifact changed')
    needed=json.loads((a.functional_release/'reservation_registry.json').read_text())['receptors'];reserved=json.loads(Path(m['registry_path']).read_text())['receptors'];keys=lambda rs:{(x['alpha_variable'],x['beta_variable']) for x in rs}
    if not keys(needed)<=keys(reserved):raise ValueError('external receptors were not reserved before fitting')
    coeff=dict(zip(model['feature_names'],model['coefficients']));records=json.loads((a.functional_release/'functional_observations.json').read_text());by_pep={};predictions=[]
    for x in records:
        query=dict(x,alpha_variable=domain(x['alpha_variable'])[0],beta_variable=domain(x['beta_variable'])[0],method='Biacore 3000',temperature_K='298.15')
        f=features(query);score=model['intercept']+sum(coeff.get(k,0)*v for k,v in f.items());key=x['mhc_allele']+':'+x['peptide'];by_pep.setdefault(key,[]).append((x['observed_activation'],score));predictions.append({'evidence_id':x['evidence_id'],'receptor_pair_id':x['receptor_pair_id'],'peptide':x['peptide'],'hla':x['mhc_allele'],'observed_activation':x['observed_activation'],'spr_reference_pkd_score':score})
    metrics=[]
    for key,values in sorted(by_pep.items()):
        y,score=np.array(values).T
        if set(y)!={0,1}:raise ValueError('one-class external peptide')
        metrics.append({'peptide_hla':key,'n':len(y),'n_positive':int(y.sum()),'auroc':float(roc_auc_score(y,score)),'average_precision':float(average_precision_score(y,score))})
    rng=np.random.default_rng(20261017);aucs=np.array([x['auroc'] for x in metrics]);ci=np.quantile(np.mean(rng.choice(aucs,(2000,len(aucs)),replace=True),axis=1),[.025,.975]).tolist()
    report={'records':len(records),'n_peptides':len(metrics),'model_run':str(a.model_run),'model_sha256':fit['model_sha256'],'external_manifest_sha256':sha256(a.functional_release/'manifest.json'),'macro_peptide_auroc':float(aucs.mean()),'macro_peptide_average_precision':float(np.mean([x['average_precision'] for x in metrics])),'macro_auroc_peptide_bootstrap_95ci':ci,'ci_scope':'Peptide resampling conditional on these 1000 TCRs/platform; not a receptor-family/study-independent clinical CI','score':'predicted pKD under a fixed hypothetical Biacore 3000 / 298.15 K reference profile; not an activation probability','external_panel_used_for_training_or_tuning':False,'physical_binding_validation_passed':False,'clinical_validation_passed':False,'limitations':['Cross-endpoint SPR-to-MIRA association only','One selected viral-panel source; HLA partly inferred; multiplex nonactivation conditional','No clinical outcomes and no probability calibration','No test-driven model changes or threshold tuning permitted after this evaluation']}
    a.out.mkdir(parents=True,exist_ok=False)
    for name,x in [('report.json',report),('per_peptide_metrics.json',metrics),('predictions.json',predictions)]: (a.out/name).write_text(json.dumps(x,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
