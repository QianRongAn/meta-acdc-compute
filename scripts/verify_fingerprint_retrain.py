"""Post-training artifact acceptance; no fitting, selection, or new test data."""
import argparse
from collections import Counter
import fcntl
import hashlib
import json
from pathlib import Path
import numpy as np
import torch
from meta_acdc.models import train_fingerprint_v2 as trainer
from scripts.run_fingerprint_retrain import validate_protocol, ROOT

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--protocol',type=Path,required=True)
    ap.add_argument('--run',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists(): raise FileExistsError(args.out)
    p=validate_protocol(args.protocol)
    if torch.__version__ != p['environment']['torch']: raise ValueError('torch version changed')
    torch.set_num_threads(p['recipe_cli']['cpu-threads'])
    torch.use_deterministic_algorithms(True)
    d=torch.load(p['dataset'],map_location='cpu',weights_only=False)
    graphs,records=d['graphs'],d['records']
    fa=json.loads((ROOT/p['fold_assignment']).read_text());fold_of=fa['per_graph_fold']
    if fa['graph_records'] != [r['complex_id'] for r in records]: raise ValueError('graph identity mismatch')
    checks=[]
    with (args.run/'execution.lock').open('r') as lock:
        fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
        if json.loads((args.run/'status.json').read_text())['state'] != 'complete':
            raise ValueError('run is incomplete')
        for arm in ['true','shuffled','replay']:
            for fold in ([0] if arm=='replay' else range(p['n_folds'])):
                for seed in (p['seeds'][:1] if arm=='replay' else p['seeds']):
                    tag=f'fold{fold}_seed{seed}';base=args.run/arm;rd=base/'model_run_v2'
                    metric=json.loads((base/f'metrics_{tag}.json').read_text())
                    for rel,expected in metric['artifact_sha256'].items():
                        if hashlib.sha256((base/rel).read_bytes()).hexdigest()!=expected:
                            raise ValueError('artifact hash mismatch')
                    checkpoint=torch.load(rd/f'{tag}.pt',map_location='cpu',weights_only=False)
                    config=checkpoint['config']
                    if (checkpoint['fold']!=fold or checkpoint['seed']!=seed or config['dataset_sha256']!=p['dataset_sha256']
                            or config['recipe_v2']['epochs']!=120 or config['device']!='cpu'):
                        raise ValueError('checkpoint identity/config mismatch')
                    if any(not torch.isfinite(v).all() for v in checkpoint['model'].values()):
                        raise ValueError('nonfinite checkpoint tensor')
                    train_idx=[i for i,f in enumerate(fold_of) if f!=fold]
                    mode='train-shuffle' if arm=='shuffled' else 'true'
                    effective_seed=p['recipe_cli']['permutation-seed']+1000003*fold+seed
                    training,mapping=trainer.training_supervision(records,train_idx,mode,effective_seed)
                    labels=json.loads((rd/f'labels_{tag}.json').read_text())
                    if labels['mapping']!=mapping or labels['mode']!=mode:
                        raise ValueError('supervision mapping mismatch')
                    if Counter((training[i]['receptor_id'],training[i]['family_id']) for i in train_idx)!=Counter(
                            (records[i]['receptor_id'],records[i]['family_id']) for i in train_idx):
                        raise ValueError('joint label counts changed')
                    if any(training[i]!=records[i] for i,f in enumerate(fold_of) if f==fold):
                        raise ValueError('test records changed')
                    model=trainer.EmbedderV2(26,16,emb_dim=p['recipe_cli']['emb-dim']).to('cpu')
                    model.load_state_dict(checkpoint['model'],strict=True)
                    embedded=trainer.embed_all(graphs,records,'cpu',model)
                    recomputed=embedded.numpy().astype(np.float32)
                    stored=np.load(rd/f'emb_{tag}.npy',allow_pickle=False)
                    S=(embedded@embedded.T).numpy().astype(np.float64)
                    stored_S=np.load(rd/f'sim_{tag}.npy',allow_pickle=False)
                    if not np.array_equal(recomputed,stored) or not np.array_equal(S,stored_S):
                        raise ValueError('checkpoint/input to embedding/similarity reproduction failed: '+arm+'/'+tag)
                    ranking=trainer.per_query_rankings(embedded,[i for i,f in enumerate(fold_of) if f==fold])
                    if ranking!=json.loads((base/f'rankings_{tag}.json').read_text()):
                        raise ValueError('ranking reproduction failed')
                    checks.append(dict(arm=arm,fold=fold,seed=seed,embedding_bitwise_equal=True,
                                       similarity_bitwise_equal=True,ranking_exact_equal=True,
                                       labels_training_only=True,joint_label_counts_preserved=True,
                                       n_changed=labels['n_changed_receptor_labels']))
                    print(f'VERIFIED {len(checks)}/31 {arm} {tag}',flush=True)
        validate_protocol(args.protocol)
    args.out.write_text(json.dumps(dict(status='pass',n_models=len(checks),checks=checks,
                                        frozen_protocol_sha256=hashlib.sha256(args.protocol.read_bytes()).hexdigest(),
                                        scope='post-training independent code path; same exposed graphs, no new independent test',
                                        measured_acceptance_sets_read=False),indent=2)+'\n')

if __name__=='__main__':main()
