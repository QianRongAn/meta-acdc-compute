"""Replay NP boundary and broad-search evidence while preserving ten prior deliveries."""
from pathlib import Path
import argparse,gzip,hashlib,json,os,subprocess,sys
from scripts.summarize_np_broad_search import report_bytes
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_boundaries_20261029'
def sha(raw):return hashlib.sha256(raw).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--replay-dir',type=Path,required=True);ap.add_argument('--anarcii-path',type=Path);ap.add_argument('--anarci-path',type=Path);ap.add_argument('--hmmerpath',type=Path);a=ap.parse_args();out=a.replay_dir.resolve();out.mkdir(parents=True,exist_ok=False);env=dict(os.environ,PYTHONPATH='src:.')
    def run(script,*extra,pythonpath=None):
        localenv=dict(env)
        if pythonpath:localenv['PYTHONPATH']=str(pythonpath)+':src:.'
        p=subprocess.run([sys.executable,'scripts/'+script,*map(str,extra)],cwd=ROOT,env=localenv,capture_output=True,text=True)
        with (out/'replay.log').open('a') as f:f.write(p.stdout+p.stderr)
        assert p.returncode==0,p.stderr[-3000:]
    public=0;local=0
    for row in json.loads((BASE/'archive_layout.json').read_text()):
        if row['packaging']=='local_only_primary_or_tool_original':
            path=BASE/row['source']
            if path.exists():assert len(path.read_bytes())==row['bytes'] and sha(path.read_bytes())==row['sha256'];local+=1
            continue
        raw=(BASE/row['packaged_file']).read_bytes()
        if row['packaging']=='gzip_exact':raw=gzip.decompress(raw)
        assert len(raw)==row['bytes'] and sha(raw)==row['sha256'],row['source'];public+=1
    run('prepare_np_broad_search.py','--out',out/'plan')
    for name in ['search_plan.json','requests.json']:assert (out/'plan'/name).read_bytes()==(BASE/'broad_search'/name).read_bytes(),name
    run('join_native_search.py','--base',BASE/'broad_search','--out',out/'search_results.json')
    assert (out/'search_results.json').read_bytes()==(BASE/'broad_search/search_results.json').read_bytes()
    run('validate_np_broad_search.py','--out',out/'primary_validation_final.json')
    assert (out/'primary_validation_final.json').read_bytes()==report_bytes()
    run('reassemble_np_search_report.py','--out',out/'reassembled_validation.json');assert (out/'reassembled_validation.json').read_bytes()==report_bytes()
    run('summarize_np_broad_search.py','--out',out/'broad_summary.json');assert (out/'broad_summary.json').read_bytes()==(BASE/'broad_cross_study_summary.json').read_bytes()
    run('audit_np_boundaries.py','--out',out/'boundary_audit.json');assert (out/'boundary_audit.json').read_bytes()==(BASE/'boundary_audit_both_methods.json').read_bytes()
    final_cells_reextracted=False
    if (BASE/'sources/final/final_mmc7.xlsx').exists() and (BASE/'sources/final/final_mmc3.xlsx').exists() and (ROOT/'docs/fingerprint_expanded_controls_20261029/sources/expanded/media-5.xlsx').exists():
        run('extract_np_final_sources.py','--out',out/'source_extract')
        for name in ['final_primary_cells.json','final_repertoire_receipt.json']:assert (out/'source_extract'/name).read_bytes()==(BASE/name).read_bytes(),name
        final_cells_reextracted=True
    numbering_recomputed=False
    if a.anarcii_path or a.anarci_path or a.hmmerpath:
        assert a.anarcii_path and a.anarci_path and a.hmmerpath
        run('number_np_constructs.py','--out',out/'neural',pythonpath=a.anarcii_path)
        assert (out/'neural/numbering_raw.json').read_bytes()==(BASE/'numbering/numbering_raw.json').read_bytes()
        run('number_np_classic.py','--out',out/'classic_raw.json','--hmmerpath',a.hmmerpath,pythonpath=a.anarci_path)
        assert (out/'classic_raw.json').read_bytes()==(BASE/'numbering/classic_raw.json').read_bytes()
        numbering_recomputed=True
    previous=[]
    folders=['fingerprint_retrain_20261023','fingerprint_hard_family_20261024','fingerprint_paired_identity_20261024','fingerprint_cross_source_20261024','fingerprint_cross_source_extension_20261024','fingerprint_audit_integration_20261027','fingerprint_positive_controls_20261027','fingerprint_construct_followup_20261027','fingerprint_native_search_20261029','fingerprint_expanded_controls_20261029']
    for folder in folders:
        base=ROOT/'docs'/folder;m=json.loads((base/'DELIVERY_MANIFEST.json').read_text())
        for name,expected in m['files'].items():assert sha((base/name).read_bytes())==expected['sha256'],(folder,name)
        previous.append(dict(folder=folder,n_files=len(m['files']),unchanged=True))
    gates=json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/DELIVERY_MANIFEST.json').read_text())
    for name,expected in {**gates['unchanged_gates'],'src/meta_acdc/structure/primary_role_evidence.py':gates['code']['src/meta_acdc/structure/primary_role_evidence.py']}.items():assert sha((ROOT/name).read_bytes())==expected['sha256'],name
    modules=['tests.test_construct_role_adjudication','tests.test_primary_role_evidence','tests.test_cross_source_extension','tests.test_cross_source_retrieval','tests.test_paired_variable_identity','tests.test_native_search','tests.test_expanded_controls','tests.test_np_boundaries']
    p=subprocess.run([sys.executable,'-m','unittest',*modules,'-v'],cwd=ROOT,env=env,capture_output=True,text=True);(out/'tests.log').write_text(p.stdout+p.stderr);assert p.returncode==0 and 'Ran 59 tests' in p.stderr and 'skipped=' not in p.stderr,p.stderr[-3000:]
    result=dict(status='pass',n_tests=59,test_skips=0,n_public_source_payloads_exact=public,n_local_only_original_payloads_hash_checked=local,final_workbook_cells_reextracted_exact=final_cells_reextracted,both_numbering_methods_recomputed_byte_exact=numbering_recomputed,plan_search_primary_boundary_and_summary_byte_exact=True,previous_deliveries=previous,n_frozen_files=sum(x['n_files'] for x in previous),role_boundary_atomic_gate_code_unchanged=True,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
