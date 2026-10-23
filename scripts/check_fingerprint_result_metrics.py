"""Post-training metric fragments must match already weight-verified embeddings."""
import argparse
import json
from pathlib import Path
import numpy as np
import torch
from meta_acdc.models.train_fingerprint import evaluate
from scripts.run_fingerprint_retrain import validate_protocol,ROOT

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--protocol',type=Path,required=True)
    ap.add_argument('--run',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.out.exists(): raise FileExistsError(a.out)
    p=validate_protocol(a.protocol);torch.set_num_threads(p['recipe_cli']['cpu-threads'])
    d=torch.load(p['dataset'],map_location='cpu',weights_only=False);records=d['records']
    folds=json.loads((ROOT/p['fold_assignment']).read_text())['per_graph_fold']
    checks=[]
    for arm in ['true','shuffled','replay']:
        for fold in ([0] if arm=='replay' else range(p['n_folds'])):
            for seed in (p['seeds'][:1] if arm=='replay' else p['seeds']):
                tag=f'fold{fold}_seed{seed}';root=a.run/arm
                E=torch.from_numpy(np.load(root/'model_run_v2'/f'emb_{tag}.npy',allow_pickle=False))
                expected=evaluate(E,records,[q for q,f in enumerate(folds) if f==fold],None)
                metric=json.loads((root/f'metrics_{tag}.json').read_text())
                if metric['fold']!=fold or metric['seed']!=seed:
                    raise ValueError('metric identity mismatch')
                for key in ['retrieval','family','domain']:
                    if expected[key]!=metric[key]: raise ValueError('metric reproduction failed: '+arm+'/'+tag+'/'+key)
                checks.append(dict(arm=arm,fold=fold,seed=seed,metrics_exact_equal=True))
                print('METRICS_VERIFIED',len(checks),arm,tag,flush=True)
    a.out.write_text(json.dumps(dict(status='pass',n_models=len(checks),checks=checks,
                                    scope='metrics recomputed from weight-verified embeddings; no fitting or new test data'),indent=2)+'\n')

if __name__=='__main__':main()
