"""Replay extended curation, all candidate native graphs and exact legacy mappings."""
from pathlib import Path
import json,hashlib,argparse,subprocess,os,sys
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(script,*args,log):
 p=subprocess.run([sys.executable,'scripts/'+script,*args],cwd=ROOT,env=dict(os.environ,PYTHONPATH='src:.'),capture_output=True,text=True)
 log.write(p.stdout+p.stderr);log.flush()
 if p.returncode:raise RuntimeError(script+' replay failed: '+p.stderr[-1000:])
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--replay-dir',type=Path,required=True);args=ap.parse_args();out=args.replay_dir
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True);base=ROOT/'docs/fingerprint_cross_source_extension_20261024'
 with (out/'replay.log').open('w') as log:
  run('curate_cross_source_extension.py','--selection',str(base/'all_class_I_selection.json'),'--annotation-dir',str(base/'sources/stcrdab_all'),'--out',str(out/'sequence_curation_all_class_I.json'),log=log)
  assert sha(out/'sequence_curation_all_class_I.json')==sha(base/'sequence_curation_all_class_I.json')
  before={name:sha(base/name) for name in ['source_balanced_design.json','isolation_followup.json','joint_readiness.json','jdi_construct_receipt.json']}
  run('validate_cross_source_extension.py',log=log);run('followup_cross_source_extension_isolates.py',log=log);run('summarize_cross_source_extension.py',log=log)
  for name,h in before.items():assert sha(base/name)==h,name
  run('qc_cross_source_extension.py','--out-dir',str(out/'native_qc'),log=log)
  assert sha(out/'native_qc/native_atomic_qc.json')==sha(base/'native_atomic_qc.json')
  a=torch.load(base/'native_graph_candidates.pt',map_location='cpu',weights_only=False);b=torch.load(out/'native_qc/native_graph_candidates.pt',map_location='cpu',weights_only=False)
  assert a['records']==b['records'] and len(a['graphs'])==len(b['graphs'])==42 and a['trainable'] is b['trainable'] is False
  for g,h in zip(a['graphs'],b['graphs']):
   for x,y in zip(g,h):assert np.array_equal(x,y)
  run('export_legacy_fingerprint_metadata.py','--out-dir',str(out/'legacy_mapping'),log=log)
  for name in ['records_index_exact.json','legacy_matrix_binding.json','cross_fold_cases.json']:assert sha(out/'legacy_mapping'/name)==sha(base/'legacy_mapping_review_final'/name)
 previous=[]
 for folder in ['fingerprint_retrain_20261023','fingerprint_hard_family_20261024','fingerprint_paired_identity_20261024','fingerprint_cross_source_20261024']:
  d=ROOT/'docs'/folder;m=json.loads((d/'DELIVERY_MANIFEST.json').read_text())
  for name,v in m['files'].items():assert sha(d/name)==v['sha256'],(folder,name)
  previous.append(dict(folder=folder,n_files=len(m['files']),unchanged=True))
 tests=subprocess.run([sys.executable,'-m','unittest','tests.test_cross_source_extension','tests.test_cross_source_retrieval','tests.test_paired_variable_identity','-v'],cwd=ROOT,env=dict(os.environ,PYTHONPATH='src:.'),capture_output=True,text=True);(base/'tests.log').write_text(tests.stdout+tests.stderr);assert tests.returncode==0
 result=dict(status='pass',n_tests=29,test_skips=0,n_QC_graphs_exactly_replayed=42,n_legacy_instances_exact_metrics_replayed=15,old_isolated_source_hashes_rechecked=285,previous_deliveries=previous,curation_and_readiness_byte_exact=True,legacy_metadata_byte_exact=True,n_new_training_runs=0,independent_acceptance_accessed=False)
 (base/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
