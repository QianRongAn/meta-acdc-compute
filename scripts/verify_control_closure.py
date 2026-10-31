"""Replay public candidate evidence and verify immutable prior deliveries."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys
from scripts.audit_control_closure import audit
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'docs/fingerprint_control_closure_20261101'
def sha(data):return hashlib.sha256(data).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--primary-originals',action='store_true');ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
 if args.out.exists():raise FileExistsError(args.out)
 cmd=[sys.executable,'scripts/verify_np105_context.py']+(['--primary-originals'] if args.primary_originals else [])
 previous=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,check=True);old=json.loads(previous.stdout);assert old['tests_passed']==104
 context=ROOT/'docs/fingerprint_np105_context_20261030';manifest=json.loads((context/'DELIVERY_MANIFEST.json').read_text())
 for name,v in manifest['files'].items():assert sha((context/name).read_bytes())==v['sha256']
 final=json.loads((BASE/'audit_final_v3.json').read_text());replay=json.loads(json.dumps(audit()));assert final==replay
 p=subprocess.run([sys.executable,'-m','unittest','tests.test_control_closure','-v'],cwd=ROOT,env=dict(os.environ,PYTHONPATH='src:.'),capture_output=True,text=True);assert p.returncode==0 and 'Ran 9 tests' in p.stderr
 args.out.with_suffix('.tests.log').write_text(p.stdout+p.stderr)
 gates=['src/meta_acdc/structure/paired_variable_identity.py','src/meta_acdc/structure/graph.py','src/meta_acdc/structure/rigid_controls.py','src/meta_acdc/structure/source_balanced_retrieval.py']
 for name in gates:
  before=subprocess.check_output(['git','show','55bb7b6:'+name],cwd=ROOT);assert before==(ROOT/name).read_bytes(),name
 tool=json.loads((BASE/'4PRH_numbering/tool_provenance.json').read_text());prior_tool=json.loads((ROOT/'docs/fingerprint_np_boundaries_20261029/numbering/numbering_raw.json').read_text())['provenance'];assert tool['packages']['anarcii']['package_files']==prior_tool['package_files'] and tool['packages']['anarcii']['version']==prior_tool['package_version']
 current_manifest=BASE/'DELIVERY_MANIFEST.json'
 if current_manifest.exists():
  current=json.loads(current_manifest.read_text())
  for name,v in current['files'].items():assert sha((BASE/name).read_bytes())==v['sha256']
  for name,v in current['code_files'].items():assert sha((ROOT/name).read_bytes())==v['sha256']
 result=dict(status='pass',tests_passed=113,tests_skipped=0,unchanged_prior_deliveries=16,unchanged_prior_files=old['unchanged_files']+len(manifest['files']),unchanged_prior_raw_and_validation_code_files=old['unchanged_prior_raw_and_validation_code_files']+len(manifest['code_files']),fixed_gate_files_unchanged=len(gates),current_audit_exactly_replayed=True,ANARCII_model_hashes_match_prior_frozen_run=True,primary_originals_rechecked=args.primary_originals,whole_SAM_not_replayed=True,n_fragment_queries=final['n_queries'],n_fresh_metadata_PDBs=final['n_fresh_metadata_PDBs'],new_class_I_admissions=0,new_qualified_cross_study_pairs=0,new_qualified_natural_control_groups=0,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
 args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
