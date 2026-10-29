"""Replay all-cell assemblies and weak-variant evidence; never reads acceptance data."""
from pathlib import Path
import argparse,gzip,hashlib,json,subprocess,sys,os,tarfile,io
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'docs/fingerprint_np105_validation_20261029';P=ROOT/'docs/fingerprint_np105_raw_20261029'
def sha(data):return hashlib.sha256(data).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--replay-dir',type=Path,required=True);args=ap.parse_args();out=args.replay_dir.resolve();out.mkdir(parents=True,exist_ok=False);selected=out/'pilot32_selected';selected.mkdir();env=dict(os.environ,PYTHONPATH='src:.')
 for mate in [1,2]:
  data=gzip.decompress((P/f'evidence/pooled_R{mate}.fastq.gz').read_bytes());(selected/f'pooled_R{mate}.fastq').write_bytes(data);lines=data.decode('ascii').splitlines(keepends=True);assert len(lines)%4==0;handles={}
  for i in range(0,len(lines),4):
   header=lines[i];bc=header.split('.BC_')[1].split('.UMI_')[0];assert len(bc)==16 and set(bc)<=set('ACGT')
   if bc not in handles:handles[bc]=(selected/f'{bc}_R{mate}.fastq').open('xb')
   handles[bc].write(''.join(lines[i:i+4]).encode())
  for f in handles.values():f.close()
 selection=json.loads((P/'selection_receipt.json').read_text());(selected/'selection_receipt.json').write_bytes((P/'selection_receipt.json').read_bytes())
 for name,x in selection['outputs'].items():assert sha((selected/name).read_bytes())==x['sha256'],name
 archived=json.loads((B/'archive_layout.json').read_text());payloads={}
 for row in archived:
  data=gzip.decompress((B/row['packaged']).read_bytes())
  if row['kind']=='gzip_exact':assert len(data)==row['bytes'] and sha(data)==row['sha256'];payloads[row['source']]=json.loads(data)
  else:
   with tarfile.open(fileobj=io.BytesIO(data)) as tar:
    members=tar.getmembers();assert {m.name for m in members}==set(row['files'])
    for member in members:
     assert member.isfile() and len(Path(member.name).parts)==2 and not Path(member.name).is_absolute();content=tar.extractfile(member).read();assert sha(content)==row['files'][member.name]['sha256']
 def run(script,*argv):
  p=subprocess.run([sys.executable,'scripts/'+script,*map(str,argv)],cwd=ROOT,env=env,capture_output=True,text=True);(out/(script+'.log')).write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr[-3000:]
 run('assemble_np105_all_cells.py','--selected',selected,'--out',out/'assemblies','--report',out/'assemblies.json')
 fasta_hashes=next(x for x in archived if x['kind']=='tar_gzip_exact_files')['files']
 for bc in selection['read_pairs_by_barcode']:
  name=bc+'/transcripts.fasta';path=out/'assemblies'/name
  if name in fasta_hashes:assert sha(path.read_bytes())==fasta_hashes[name]['sha256'],bc
  else:assert not path.exists(),bc
 run('number_np105_all_cells.py','--assemblies',out/'assemblies.json','--out',out/'numbering.json');numbered=json.loads((out/'numbering.json').read_text());assert numbered['cells']==payloads['all92_numbering.json']['cells'];assert numbered['cell_counts']==payloads['all92_numbering.json']['cell_counts']
 run('validate_np105_raw.py','--raw-root',out,'--out',out/'base_validation.json')
 run('audit_np105_cell_variant.py','--numbering',out/'numbering.json','--sam',out/'pilot32_validation/observed_reads.sam','--out',out/'variant.json');assert json.loads((out/'variant.json').read_text())==json.loads((B/'variant_raw_support.json').read_text())
 run('check_np105_paired_cell_support.py','--sam',out/'pilot32_validation/observed_reads.sam','--out',out/'cell_support.json');assert json.loads((out/'cell_support.json').read_text())==json.loads((B/'paired_cells_base_support.json').read_text())
 previous=[]
 for manifest in sorted((ROOT/'docs').glob('fingerprint_*/DELIVERY_MANIFEST.json')):
  if manifest.parent==B:continue
  m=json.loads(manifest.read_text())
  for name,x in m['files'].items():assert sha((manifest.parent/name).read_bytes())==x['sha256'],str(manifest.parent/name)
  previous.append(dict(folder=manifest.parent.name,n_files=len(m['files']),unchanged=True))
 oldcode=json.loads((P/'DELIVERY_MANIFEST.json').read_text())['code']
 for name,x in oldcode.items():assert sha((ROOT/name).read_bytes())==x['sha256'],name
 gates=json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/DELIVERY_MANIFEST.json').read_text())
 for name,x in {**gates['unchanged_gates'],'src/meta_acdc/structure/primary_role_evidence.py':gates['code']['src/meta_acdc/structure/primary_role_evidence.py']}.items():assert sha((ROOT/name).read_bytes())==x['sha256']
 modules=['tests.test_construct_role_adjudication','tests.test_primary_role_evidence','tests.test_cross_source_extension','tests.test_cross_source_retrieval','tests.test_paired_variable_identity','tests.test_native_search','tests.test_expanded_controls','tests.test_np_boundaries','tests.test_np_adjudication','tests.test_np105_raw','tests.test_np105_validation']
 test=subprocess.run([sys.executable,'-m','unittest',*modules,'-v'],cwd=ROOT,env=env,capture_output=True,text=True);(out/'tests.log').write_text(test.stdout+test.stderr);assert test.returncode==0 and 'Ran 83 tests' in test.stderr and 'skipped=' not in test.stderr
 receipt=dict(status='pass',tests=83,test_skips=0,all92_transcript_outputs_and_absences_reproduced=True,n_transcript_FASTAs=sum(name.endswith('/transcripts.fasta') for name in fasta_hashes),all92_numbering_evidence_replayed=True,variant_codon_evidence_replayed=True,per_cell_every_base_Q30_support_replayed=True,original_selected_read_hashes_all_reconstructed=True,previous_deliveries=previous,n_previous_deliveries=len(previous),n_frozen_files=sum(x['n_files'] for x in previous),previous_raw_code_files_unchanged=len(oldcode),fixed_gates_unchanged=True,whole_run_MD5_verified=False,new_complete_identity_admissions=0,new_qualified_cross_study_pairs=0,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
 (out/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
