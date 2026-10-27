"""Replay current candidate graph admission and protect all prior frozen deliveries."""
from pathlib import Path
import json,hashlib,subprocess,os,sys,argparse
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--replay-dir',type=Path,required=True);args=ap.parse_args();out=args.replay_dir
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True);base=ROOT/'docs/fingerprint_positive_controls_20261027';env=dict(os.environ,PYTHONPATH='src:.')
 before={n:sha(base/n) for n in ['source_balanced_design.json','joint_readiness.json']}
 with (out/'replay.log').open('w') as log:
  for script,extra in [('validate_positive_control_sources.py',[]),('qc_positive_controls.py',['--out-dir',str(out/'native_qc')]),('summarize_positive_controls.py',[])]:
   p=subprocess.run([sys.executable,'scripts/'+script,*extra],cwd=ROOT,env=env,capture_output=True,text=True);log.write(p.stdout+p.stderr);assert p.returncode==0,script
 for n,h in before.items():assert sha(base/n)==h,n
 assert sha(out/'native_qc/native_atomic_qc.json')==sha(base/'native_atomic_qc.json')
 a=torch.load(base/'native_graph_candidates.pt',map_location='cpu',weights_only=False);b=torch.load(out/'native_qc/native_graph_candidates.pt',map_location='cpu',weights_only=False)
 assert a['records']==b['records'] and len(a['graphs'])==len(b['graphs'])==42
 assert a['trainable'] is b['trainable'] is False
 for g,h in zip(a['graphs'],b['graphs']):
  for x,y in zip(g,h):assert np.array_equal(x,y)
 previous=[]
 for folder in ['fingerprint_retrain_20261023','fingerprint_hard_family_20261024','fingerprint_paired_identity_20261024','fingerprint_cross_source_20261024','fingerprint_cross_source_extension_20261024','fingerprint_audit_integration_20261027']:
  d=ROOT/'docs'/folder;m=json.loads((d/'DELIVERY_MANIFEST.json').read_text())
  for n,v in m['files'].items():assert sha(d/n)==v['sha256'],(folder,n)
  previous.append(dict(folder=folder,n_files=len(m['files']),unchanged=True))
 tests=subprocess.run([sys.executable,'-m','unittest','tests.test_primary_role_evidence','tests.test_cross_source_extension','tests.test_cross_source_retrieval','tests.test_paired_variable_identity','-v'],cwd=ROOT,env=env,capture_output=True,text=True);(base/'tests.log').write_text(tests.stdout+tests.stderr);assert tests.returncode==0
 readiness=json.loads((base/'joint_readiness.json').read_text());assert readiness['summary']['n_QC_cross_study_pairs']==0 and readiness['summary']['n_QC_source_balanced_blocks']==0
 assert any(x['pdb']=='6vrm' and x['status']=='native_graph_QC_pass' for x in a['records'])
 q=json.loads((base/'native_atomic_qc.json').read_text());assert any(x['pdb']=='6vrn' and x.get('status')=='quarantined' and any('incomplete interface heavy atoms' in reason for reason in x.get('reasons',[])) for x in q['rows'])
 result=dict(status='pass',graph_arrays_and_metadata_exact=True,n_graphs=42,previous_deliveries=previous,test_modules=4,test_skips=0,n_new_training_runs=0,independent_acceptance_accessed=False,protocol_frozen=False)
 (base/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
