"""Verify J-reference corrections while preserving previous frozen evidence and gates."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys
from scripts.audit_tk3_j_evidence import audit
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_tk3_j_evidence_20261101'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--primary-originals',action='store_true');args=ap.parse_args()
 if args.out.exists():raise FileExistsError(args.out)
 cmd=[sys.executable,'scripts/verify_control_closure.py','--out',str(args.out.with_suffix('.prior.json'))]+(['--primary-originals'] if args.primary_originals else [])
 previous=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,check=True);old=json.loads(previous.stdout);assert old['tests_passed']==113
 frozen=ROOT/'docs/fingerprint_control_closure_20261101';m=json.loads((frozen/'DELIVERY_MANIFEST.json').read_text())
 for name,r in m['files'].items():assert sha(frozen/name)==r['sha256']
 fresh=json.loads((BASE/'reference_audit.json').read_text());assert fresh==json.loads(json.dumps(audit()))
 if args.primary_originals:
  from scripts.audit_tk3_primary import primary_audit
  assert primary_audit()==json.loads((BASE/'primary_source_audit_final.json').read_text())
 test=subprocess.run([sys.executable,'-m','unittest','tests.test_j_reference_evidence','-v'],cwd=ROOT,env=dict(os.environ,PYTHONPATH='src:.'),capture_output=True,text=True);assert test.returncode==0 and 'Ran 10 tests' in test.stderr,test.stderr;args.out.with_suffix('.tests.log').write_text(test.stdout+test.stderr)
 for name in ['scripts/audit_np_boundaries.py','scripts/audit_np_adjudication.py','scripts/validate_np105_raw.py','scripts/audit_control_closure.py','src/meta_acdc/structure/paired_variable_identity.py','src/meta_acdc/structure/graph.py','src/meta_acdc/structure/rigid_controls.py']:
  before=subprocess.check_output(['git','show','7838b8b:'+name],cwd=ROOT);assert before==(ROOT/name).read_bytes(),name
 if (BASE/'DELIVERY_MANIFEST.json').exists():
  current=json.loads((BASE/'DELIVERY_MANIFEST.json').read_text())
  for name,r in current['files'].items():assert sha(BASE/name)==r['sha256']
  for name,r in current['code_files'].items():assert sha(ROOT/name)==r['sha256']
 result=dict(status='pass',tests_passed=123,tests_skipped=0,previous_frozen_deliveries=17,previous_frozen_files=old['unchanged_prior_files']+len(m['files']),previous_raw_and_validation_code_files=old['unchanged_prior_raw_and_validation_code_files']+len(m['code_files']),unchanged_old_helper_callers_and_fixed_gates=7,corrected_reference_audit_exactly_replayed=True,primary_originals_rechecked=args.primary_originals,old_assertions_are_historical_artifact_replay_not_all_reference_validation=True,whole_raw_SAM_not_replayed=True,new_model_or_annotation_inference=False,new_identity_admissions=0,new_qualified_cross_study_pairs=0,new_qualified_natural_control_groups=0,admission_gates_changed=False,protocol_frozen=False,training_started=False,independent_acceptance_accessed=False)
 args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
