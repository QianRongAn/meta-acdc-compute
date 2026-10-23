"""Post-hoc error descriptions on frozen DEVELOPMENT outputs; no model selection."""
from pathlib import Path
import json
import numpy as np
import torch
from scripts.run_paired_identity_development import baseline_matrices
ROOT=Path(__file__).resolve().parents[1]
def main():
    base=ROOT/'docs/fingerprint_paired_identity_20261024';p=json.loads((base/'development/protocol.json').read_text());d=torch.load(ROOT/p['dataset'],map_location='cpu',weights_only=False);out=[];geometry=[]
    for f in p['folds']:
        idx=f['validation_indices'];matrices,_=baseline_matrices(d,f['fit_indices'])
        for q in idx:
            record=d['records'][q];pos=[j for j in idx if j!=q and d['records'][j]['receptor_id']==record['receptor_id']];neg=[j for j in idx if d['records'][j]['receptor_id']!=record['receptor_id']]
            geometry.append(dict(index=q,source_index=record['source_index'],complex_id=record['complex_id'],pdb=record['pdb'],family_id=record['family_id'],n_nodes=len(d['graphs'][q][0]),n_edges=len(d['graphs'][q][2]),n_positive=len(pos),n_negative=len(neg)))
            if not pos or not neg:continue
            for arm in p['arms']:
                for seed in p['seeds']:
                    s=np.load(base/'development/run'/f['fold_id']/arm/str(seed)/'similarities.npy')
                    a=np.asarray(s[q,pos],float);b=np.asarray(s[q,neg],float);all_scores=np.concatenate((a,b))
                    auc=float((a[:,None]>b[None,:]).mean()+.5*(a[:,None]==b[None,:]).mean())
                    out.append(dict(fold_id=f['fold_id'],query_index=q,complex_id=record['complex_id'],arm=arm,seed=seed,auroc=auc,mean_positive_cosine=float(a.mean()),mean_negative_cosine=float(b.mean()),mean_positive_minus_negative=float(a.mean()-b.mean()),score_range=float(np.ptp(all_scores)),exactly_constant=bool(np.ptp(all_scores)==0),positive_complex_ids=[d['records'][j]['complex_id'] for j in pos],negative_complex_ids=[d['records'][j]['complex_id'] for j in neg],composition_mean_positive_minus_negative=float(matrices['composition'][q,pos].mean()-matrices['composition'][q,neg].mean()),size_mean_positive_minus_negative=float(matrices['graph_size'][q,pos].mean()-matrices['graph_size'][q,neg].mean())))
    result=dict(scope='post-hoc descriptive error analysis only; no changes to labels, partitions, hyperparameters or seed selection',n_query_seed_arm_rows=len(out),n_exactly_constant_rows=sum(x['exactly_constant'] for x in out),row_metrics=out,graph_metadata=geometry)
    (base/'corrected_error_analysis.json').write_text(json.dumps(result,indent=2)+'\n');print('ERROR_ROWS',len(out),'EXACT_CONSTANT',result['n_exactly_constant_rows'])
    for f in p['folds']:
        print(f['fold_id'],[(a,[x['auroc'] for x in out if x['fold_id']==f['fold_id'] and x['arm']==a]) for a in p['arms']])
if __name__=='__main__':main()
