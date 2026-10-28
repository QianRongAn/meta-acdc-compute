"""Replay finite discovery offline and verify eight frozen deliveries without rewriting them."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'docs/fingerprint_native_search_20261029'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--replay-dir',type=Path,required=True); args=ap.parse_args()
    out=args.replay_dir.resolve(); out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,PYTHONPATH='src:.')
    def run(script,*extra):
        p=subprocess.run([sys.executable,'scripts/'+script,*map(str,extra)],cwd=ROOT,env=env,text=True,capture_output=True)
        with (out/'replay.log').open('a') as f: f.write(p.stdout+p.stderr)
        assert p.returncode==0,p.stderr
    run('prepare_native_search.py','--out',out/'plan')
    for name in ['search_plan.json','requests.json']: assert sha(out/'plan'/name)==sha(BASE/name),name
    for name,expected in json.loads((BASE/'search_plan.json').read_text())['source_files'].items(): assert sha(ROOT/name)==expected,name
    n_sources=0
    for folder in ['search','metadata']:
        for row in json.loads((BASE/'sources'/folder/'manifest.json').read_text()):
            assert row['status']==200
            assert sha(BASE/'sources'/folder/row['file'])==row['sha256']
            n_sources+=1
    run('join_native_search.py','--base',BASE,'--out',out/'search_results.json')
    run('validate_native_search.py','--out',out/'primary_validation_final.json')
    for name in ['search_results.json','primary_validation_final.json']: assert sha(out/name)==sha(BASE/name),name
    previous=[]
    folders=['fingerprint_retrain_20261023','fingerprint_hard_family_20261024','fingerprint_paired_identity_20261024','fingerprint_cross_source_20261024','fingerprint_cross_source_extension_20261024','fingerprint_audit_integration_20261027','fingerprint_positive_controls_20261027','fingerprint_construct_followup_20261027']
    for folder in folders:
        base=ROOT/'docs'/folder; manifest=json.loads((base/'DELIVERY_MANIFEST.json').read_text())
        for name,expected in manifest['files'].items(): assert sha(base/name)==expected['sha256'],(folder,name)
        previous.append(dict(folder=folder,n_files=len(manifest['files']),unchanged=True))
    gates=json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/DELIVERY_MANIFEST.json').read_text())
    for name,expected in {**gates['unchanged_gates'],'src/meta_acdc/structure/primary_role_evidence.py':gates['code']['src/meta_acdc/structure/primary_role_evidence.py']}.items(): assert sha(ROOT/name)==expected['sha256'],name
    modules=['tests.test_construct_role_adjudication','tests.test_primary_role_evidence','tests.test_cross_source_extension','tests.test_cross_source_retrieval','tests.test_paired_variable_identity','tests.test_native_search']
    p=subprocess.run([sys.executable,'-m','unittest',*modules,'-v'],cwd=ROOT,env=env,text=True,capture_output=True)
    (out/'tests.log').write_text(p.stdout+p.stderr); assert p.returncode==0,p.stderr
    assert 'Ran 47 tests' in p.stderr and 'skipped=' not in p.stderr
    report=dict(status='pass',byte_exact_plan_search_and_primary_replay=True,n_raw_sources_verified=n_sources,previous_deliveries=previous,n_frozen_files=sum(x['n_files'] for x in previous),unchanged_role_boundary_atomic_gates=True,n_tests=47,test_skips=0,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n'); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
