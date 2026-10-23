"""Execute the frozen corrected fingerprint experiment, one process at a time."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import time
import traceback
import torch
from meta_acdc.models import train_fingerprint_v2 as trainer

ROOT = Path(__file__).resolve().parents[1]

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def validate_protocol(path):
    expected = path.with_suffix('.sha256').read_text().split()[0]
    if digest(path) != expected:
        raise ValueError('protocol lock mismatch')
    p = json.loads(path.read_text())
    for rel, sha in p['source_sha256'].items():
        if digest(ROOT/rel) != sha:
            raise ValueError('frozen source changed: '+rel)
    if digest(p['dataset']) != p['dataset_sha256']:
        raise ValueError('dataset hash mismatch')
    if digest(ROOT/p['fold_assignment']) != p['fold_assignment_sha256']:
        raise ValueError('fold assignment changed')
    return p

def arguments(p, out, mode, fold, seeds):
    a = ['--dataset',p['dataset'],'--out-dir',str(out),'--device','cpu',
         '--label-mode',mode,'--only-fold',str(fold),'--resume','--seeds'] + [str(s) for s in seeds]
    for key,value in p['recipe_cli'].items():
        a += ['--'+key,str(value)]
    return a

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--protocol',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True)
    args = ap.parse_args()
    p = validate_protocol(args.protocol)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    with (args.out_dir/'execution.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        identity = {'protocol_sha256':digest(args.protocol),'pid':os.getpid()}
        binding = args.out_dir/'binding.json'
        if binding.exists() and json.loads(binding.read_text())['protocol_sha256'] != identity['protocol_sha256']:
            raise ValueError('foreign output directory')
        binding.write_text(json.dumps(identity,indent=2)+'\n')
        status_path = args.out_dir/'status.json'
        def status(state, **kwargs):
            status_path.write_text(json.dumps(dict(state=state,wall_time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                                                  protocol_sha256=identity['protocol_sha256'],**kwargs),indent=2)+'\n')
        try:
            for fold in range(p['n_folds']):
                for arm,mode in p['arms'].items():
                    validate_protocol(args.protocol)
                    status('running',arm=arm,fold=fold)
                    start=time.monotonic()
                    trainer.main(arguments(p,args.out_dir/arm,mode,fold,p['seeds']))
                    with (args.out_dir/'phase_timing.jsonl').open('a') as log:
                        log.write(json.dumps(dict(arm=arm,fold=fold,seconds=time.monotonic()-start))+'\n')
            validate_protocol(args.protocol)
            status('determinism_replay')
            trainer.main(arguments(p,args.out_dir/'replay', 'true',0,[p['seeds'][0]]))
            tag=f"fold0_seed{p['seeds'][0]}"
            a=torch.load(args.out_dir/'true/model_run_v2'/f'{tag}.pt',map_location='cpu',weights_only=False)['model']
            b=torch.load(args.out_dir/'replay/model_run_v2'/f'{tag}.pt',map_location='cpu',weights_only=False)['model']
            import numpy as np
            same_weights = all(torch.equal(a[k],b[k]) for k in a)
            same_emb = np.array_equal(np.load(args.out_dir/'true/model_run_v2'/f'emb_{tag}.npy'),
                                      np.load(args.out_dir/'replay/model_run_v2'/f'emb_{tag}.npy'))
            replay={'weights_bitwise_equal':same_weights,'embeddings_bitwise_equal':same_emb,
                    'scope':'full fold0 first seed on same CPU/software; not cross-device guarantee'}
            (args.out_dir/'replay_check.json').write_text(json.dumps(replay,indent=2)+'\n')
            if not same_weights or not same_emb:
                raise RuntimeError('full determinism replay failed')
            status('evaluating')
            from scripts.evaluate_fingerprint_retrain import evaluate_run
            evaluate_run(p,args.out_dir)
            validate_protocol(args.protocol)
            status('complete',n_main_instances=len(p['arms'])*p['n_folds']*len(p['seeds']),n_replays=1)
        except BaseException as error:
            status('failed',error=repr(error))
            traceback.print_exc()
            raise

if __name__ == '__main__':
    main()
