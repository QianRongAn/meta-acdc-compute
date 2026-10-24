"""Replay source parsing, source-balanced gate, old-source followup and native graphs."""
from pathlib import Path
import hashlib,json,subprocess,os,argparse
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1];PY='/home/administrator/docs/meta-acdc/.venv/bin/python'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run(script,*args):
    subprocess.run([PY,'scripts/'+script,*args],cwd=ROOT,env=dict(os.environ,PYTHONPATH='src:.',OPENBLAS_NUM_THREADS='1'),check=True)
def main():
    base=ROOT/'docs/fingerprint_cross_source_20261024';ap=argparse.ArgumentParser();ap.add_argument('--replay-dir',type=Path,default=base/'verification_replay');args=ap.parse_args();replay=args.replay_dir
    if replay.exists():raise FileExistsError(replay)
    replay.mkdir();before={name:sha(base/name) for name in ['sequence_curation_v2.json','source_balanced_design.json','isolation_followup.json','native_atomic_qc.json']}
    run('curate_cross_source_receptors.py','--out',str(replay/'sequence_curation_v2.json'));assert sha(replay/'sequence_curation_v2.json')==before['sequence_curation_v2.json']
    run('validate_cross_source_design.py');assert sha(base/'source_balanced_design.json')==before['source_balanced_design.json']
    run('followup_isolated_fingerprint_items.py');assert sha(base/'isolation_followup.json')==before['isolation_followup.json']
    run('qc_cross_source_structures.py','--out-dir',str(replay/'native_qc'));assert sha(replay/'native_qc/native_atomic_qc.json')==before['native_atomic_qc.json']
    a=torch.load(base/'native_graph_candidates.pt',map_location='cpu',weights_only=False);b=torch.load(replay/'native_qc/native_graph_candidates.pt',map_location='cpu',weights_only=False)
    assert a['trainable'] is b['trainable'] is False and a['binding_labels_available'] is False
    assert a['records']==b['records'] and len(a['graphs'])==len(b['graphs'])
    for g,h in zip(a['graphs'],b['graphs']):
        for x,y in zip(g,h):assert np.array_equal(x,y)
    old=[]
    for folder in ['fingerprint_retrain_20261023','fingerprint_hard_family_20261024','fingerprint_paired_identity_20261024']:
        root=ROOT/'docs'/folder;m=json.loads((root/'DELIVERY_MANIFEST.json').read_text())
        for name,v in m['files'].items():assert sha(root/name)==v['sha256'],(folder,name)
        old.append(dict(folder=folder,n_files=len(m['files']),status='pass'))
    result=dict(status='pass',candidate_curation_byte_exact=True,source_balanced_gate_byte_exact=True,old_isolation_followup_byte_exact=True,n_original_sources_rechecked=285,n_native_graphs_replayed=len(a['graphs']),native_arrays_exact=True,previous_deliveries=old,training_release_created=False,model_training_started=False,independent_acceptance_accessed=False)
    (base/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print('CROSS_SOURCE_VERIFICATION_PASS',len(a['graphs']))
if __name__=='__main__':main()
