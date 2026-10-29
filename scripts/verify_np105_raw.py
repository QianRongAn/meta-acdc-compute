"""Replay archived read-only NP105 evidence without touching acceptance data or prior gates."""
from pathlib import Path
import argparse,gzip,hashlib,json,subprocess,sys,os
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'docs/fingerprint_np105_raw_20261029'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--replay-dir',type=Path,required=True);a=ap.parse_args();out=a.replay_dir.resolve();out.mkdir(parents=True,exist_ok=False);selected=out/'pilot32_selected';selected.mkdir();sources=0
 for row in json.loads((B/'archive_layout.json').read_text()):
  data=gzip.decompress((B/row['packaged']).read_bytes());assert len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'];(selected/row['source']).write_bytes(data);sources+=1
 old=json.loads((B/'selection_receipt.json').read_text())
 for name,x in old['outputs'].items():
  if (selected/name).is_file():assert sha(selected/name)==x['sha256']
 top=old['top_three_by_preassembly_read_count'];env=dict(os.environ,NP105_SPADES_ROOT='/home/administrator/docs/np105-raw-20261029/spades-tools',LD_LIBRARY_PATH='/home/administrator/docs/np105-raw-20261029/spades-tools/usr/lib/x86_64-linux-gnu',PYTHONPATH='src:.')
 replay=[]
 for name,stem in [('pilot32_pool_assembly_v2','pooled')]+[('pilot32_cell_'+bc,bc) for bc in top]:
  dest=out/name;cmd=[sys.executable,'scripts/run_np105_spades.py','--rna','--only-assembler','-t','2','-m','4','-1',str(selected/(stem+'_R1.fastq')),'-2',str(selected/(stem+'_R2.fastq')),'-o',str(dest)]
  with (out/(name+'.log')).open('w') as f:subprocess.run(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
  assert sha(dest/'transcripts.fasta')==sha(B/'evidence'/name/'transcripts.fasta');replay.append(name)
 with (out/'validation.log').open('w') as f:subprocess.run([sys.executable,'scripts/validate_np105_raw.py','--raw-root',str(out),'--out',str(out/'validation.json')],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 oldv=json.loads((B/'raw_sequence_validation.json').read_text());new=json.loads((out/'validation.json').read_text());assert oldv['domains']==new['domains'] and oldv['alignment_counts']==new['alignment_counts']
 previous=[]
 for base in sorted((ROOT/'docs').glob('fingerprint_*/DELIVERY_MANIFEST.json')):
  if base.parent==B:continue
  m=json.loads(base.read_text())
  for name,x in m['files'].items():assert sha(base.parent/name)==x['sha256'],str(base.parent/name)
  previous.append(dict(folder=base.parent.name,n_files=len(m['files']),unchanged=True))
 gates=json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/DELIVERY_MANIFEST.json').read_text())
 for name,x in {**gates['unchanged_gates'],'src/meta_acdc/structure/primary_role_evidence.py':gates['code']['src/meta_acdc/structure/primary_role_evidence.py']}.items():assert sha(ROOT/name)==x['sha256']
 modules=['tests.test_construct_role_adjudication','tests.test_primary_role_evidence','tests.test_cross_source_extension','tests.test_cross_source_retrieval','tests.test_paired_variable_identity','tests.test_native_search','tests.test_expanded_controls','tests.test_np_boundaries','tests.test_np_adjudication','tests.test_np105_raw']
 t=subprocess.run([sys.executable,'-m','unittest',*modules,'-v'],cwd=ROOT,env=env,capture_output=True,text=True);(out/'tests.log').write_text(t.stdout+t.stderr);assert t.returncode==0 and 'Ran 75 tests' in t.stderr and 'skipped=' not in t.stderr
 hits={role:json.loads((B/'sources/fresh_search'/f'{role}.json').read_text())['result_set'] for role in ['TRA','TRB']};assert [x['identifier'] for x in hits['TRA']]==['9NVA_1'] and [x['identifier'] for x in hits['TRB']]==['9NVA_2']
 receipt=dict(status='pass',tests=75,test_skips=0,compressed_read_payloads_hash_checked=sources,assemblies_byte_exact=replay,per_base_validation_replayed=True,archived_selection_outputs_hash_checked=True,previous_deliveries=previous,n_previous_deliveries=len(previous),n_frozen_files=sum(x['n_files'] for x in previous),fixed_gate_code_unchanged=True,whole_run_MD5_verified=False,new_complete_identity_admissions=0,new_qualified_cross_study_pairs=0,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
 (out/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
