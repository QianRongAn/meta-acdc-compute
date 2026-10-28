"""Replay NP source adjudication and current-coordinate QC; keep eleven prior deliveries unchanged."""
from pathlib import Path
import argparse,gzip,hashlib,json,os,subprocess,sys
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_adjudication_20261029'
def sha(raw):return hashlib.sha256(raw).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--replay-dir',type=Path,required=True);a=ap.parse_args();out=a.replay_dir.resolve();out.mkdir(parents=True,exist_ok=False);env=dict(os.environ,PYTHONPATH='src:.')
    def run(script,*extra):
        p=subprocess.run([sys.executable,'scripts/'+script,*map(str,extra)],cwd=ROOT,env=env,capture_output=True,text=True)
        with (out/'replay.log').open('a') as f:f.write(p.stdout+p.stderr)
        assert p.returncode==0,p.stderr[-3000:]
    n_sources=0;coordinate_payloads={}
    for row in json.loads((BASE/'archive_layout.json').read_text()):
        raw=(BASE/row['packaged_file']).read_bytes()
        if row['packaging']=='gzip_exact':raw=gzip.decompress(raw)
        assert len(raw)==row['bytes'] and sha(raw)==row['sha256'],row['source'];n_sources+=1
        if row['source'].endswith('.pdb'):coordinate_payloads[Path(row['source']).stem]=raw
    run('audit_np_adjudication.py','--out',out/'adjudication.json');assert (out/'adjudication.json').read_bytes()==(BASE/'adjudication_source_tiers.json').read_bytes()
    result=json.loads((out/'adjudication.json').read_text());assert result['n_chain_role_sources_supported']==12 and result['n_FR1_gap_conflicts_source_adjudicated']==4
    np105=next(r for r in result['NP105_sample_results'] if r['sample_label']=='D4');assert np105['n_unambiguous_quality_paired_barcodes']==92 and result['n_new_fixed_identity_admissions']==0
    run('prepare_adjudication_cross_QC.py','--out',out/'fresh_identity.json');assert (out/'fresh_identity.json').read_bytes()==(BASE/'fresh_cross_identity.json').read_bytes()
    identity=json.loads((out/'fresh_identity.json').read_text());(out/'coordinates').mkdir()
    for row in identity['rows']:
        path=out/'coordinates'/(row['pdb']+'.pdb');raw=coordinate_payloads[row['pdb']];assert sha(raw)==row['primary_header_sha256'];path.write_bytes(raw);row['primary_header']=str(path)
    (out/'replay_identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    run('qc_construct_followup.py','--identity',out/'replay_identity.json','--out',out/'native_QC')
    assert (out/'native_QC/native_atomic_qc.json').read_bytes()==(BASE/'fresh_cross_QC/native_atomic_qc.json').read_bytes()
    first=torch.load(BASE/'fresh_cross_QC/native_graph_candidates.pt',map_location='cpu',weights_only=False);second=torch.load(out/'native_QC/native_graph_candidates.pt',map_location='cpu',weights_only=False)
    assert first['records']==second['records'] and len(first['graphs'])==len(second['graphs'])==2
    assert first['trainable'] is second['trainable'] is False
    for x,y in zip(first['graphs'],second['graphs']):
        for p,q in zip(x,y):assert np.array_equal(p,q)
    run('summarize_adjudication_cross_QC.py','--out',out/'cross_summary.json');assert (out/'cross_summary.json').read_bytes()==(BASE/'fresh_cross_summary.json').read_bytes()
    previous=[]
    folders=['fingerprint_retrain_20261023','fingerprint_hard_family_20261024','fingerprint_paired_identity_20261024','fingerprint_cross_source_20261024','fingerprint_cross_source_extension_20261024','fingerprint_audit_integration_20261027','fingerprint_positive_controls_20261027','fingerprint_construct_followup_20261027','fingerprint_native_search_20261029','fingerprint_expanded_controls_20261029','fingerprint_np_boundaries_20261029']
    for folder in folders:
        base=ROOT/'docs'/folder;m=json.loads((base/'DELIVERY_MANIFEST.json').read_text())
        for name,expected in m['files'].items():assert sha((base/name).read_bytes())==expected['sha256'],(folder,name)
        previous.append(dict(folder=folder,n_files=len(m['files']),unchanged=True))
    gates=json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/DELIVERY_MANIFEST.json').read_text())
    for name,expected in {**gates['unchanged_gates'],'src/meta_acdc/structure/primary_role_evidence.py':gates['code']['src/meta_acdc/structure/primary_role_evidence.py']}.items():assert sha((ROOT/name).read_bytes())==expected['sha256'],name
    modules=['tests.test_construct_role_adjudication','tests.test_primary_role_evidence','tests.test_cross_source_extension','tests.test_cross_source_retrieval','tests.test_paired_variable_identity','tests.test_native_search','tests.test_expanded_controls','tests.test_np_boundaries','tests.test_np_adjudication']
    p=subprocess.run([sys.executable,'-m','unittest',*modules,'-v'],cwd=ROOT,env=env,capture_output=True,text=True);(out/'tests.log').write_text(p.stdout+p.stderr);assert p.returncode==0 and 'Ran 67 tests' in p.stderr and 'skipped=' not in p.stderr,p.stderr[-3000:]
    receipt=dict(status='pass',n_tests=67,test_skips=0,n_public_source_payloads_hash_checked=n_sources,n_NP105_strict_same_cell_pair_barcodes=92,n_FR1_gap_conflicts_source_adjudicated=4,n_current_coordinate_PDBs_replayed=5,n_existing_graphs_array_exactly_replayed=2,adjudication_and_QC_and_summary_byte_exact=True,previous_deliveries=previous,n_frozen_files=sum(r['n_files'] for r in previous),role_boundary_atomic_gate_code_unchanged=True,n_new_complete_identity_admissions=0,n_new_cross_study_qualified_pairs=0,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    (out/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
