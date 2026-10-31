"""Fixed-tool full-input annotation of one rescue lead; no identity admission."""
from pathlib import Path
import argparse,hashlib,importlib.metadata,json,sys
sys.path[:0]=['/tmp/np-classic-deps-20261029','/tmp/np-numbering-deps-20261029']
import anarci
from anarcii import Anarcii
from scripts.audit_np_boundaries import checked_numbering
from meta_acdc.structure.paired_variable_identity import germline_references
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_control_closure_20261101'
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=False);plan=json.loads((BASE/'4PRH_numbering_plan.json').read_text());inputs=plan['inputs']
 fasta=''.join('>'+x['name']+'\n'+x['sequence']+'\n' for x in inputs);(out/'inputs.fasta').write_text(fasta)
 classic,details,hits=anarci.anarci([(x['name'],x['sequence']) for x in inputs],scheme='imgt',allow={'A','B'},hmmerpath='/tmp/np-hmmer-bin-20261029/usr/bin',ncpu=2,assign_germline=True,allowed_species=['human'],bit_score_threshold=80)
 model=Anarcii(seq_type='tcr',mode='accuracy',cpu=True,ncpu=2,batch_size=8);neural=model.number(str(out/'inputs.fasta'));refs_dir=ROOT/'docs/fingerprint_paired_identity_20261024/sources';refs=germline_references(list(refs_dir.glob('Homo_sapiens_TR*.fasta')));rows=[]
 for i,x in enumerate(inputs):
  expected='A' if x['role']=='alpha' else 'B';methods={}
  try:
   assert classic[i] and len(classic[i])==1;mapping,start,end=classic[i][0];methods['ANARCI']=checked_numbering(x['sequence'],mapping,start,end,details[i][0]['chain_type'],expected)
  except (ValueError,AssertionError) as error:methods['ANARCI']=dict(error=str(error) or type(error).__name__)
  result=neural[x['name']]
  try:
   if result.get('error'):raise ValueError(result['error'])
   methods['ANARCII']=checked_numbering(x['sequence'],result['numbering'],result['query_start'],result['query_end'],result['chain_type'],expected)
  except ValueError as error:methods['ANARCII']=dict(error=str(error))
  supported=[v for v in methods.values() if 'sequence' in v];sequence=supported[0]['sequence'] if supported else '';vp='TRAV' if expected=='A' else 'TRBV';jp='TRAJ' if expected=='A' else 'TRBJ'
  rows.append(dict(name=x['name'],role=x['role'],full_input_sha256=hashlib.sha256(x['sequence'].encode()).hexdigest(),methods=methods,exact_two_method_position_agreement=len(supported)==2 and supported[0]==supported[1],functional_V_first10_matches=sorted(g for (sp,g),values in refs.items() if sp=='human' and g.startswith(vp) and sequence and any(sequence.startswith(v['amino_acids'][:10]) for v in values)),functional_J_last10_matches=sorted(g for (sp,g),values in refs.items() if sp=='human' and g.startswith(jp) and sequence and any(sequence.endswith(v['amino_acids'][-10:]) for v in values)),identity_admitted=False))
 raw=dict(inputs=inputs,classic=classic,classic_details=details,classic_hits=hits,neural=neural,versions={x:importlib.metadata.version(x) for x in ['anarci','anarcii']},parameters=plan['parameters'],input_plan_sha256=hashlib.sha256((BASE/'4PRH_numbering_plan.json').read_bytes()).hexdigest(),annotation_only=True)
 (out/'raw.json').write_text(json.dumps(raw,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(dict(rows=rows,trimmed_inputs=0,imputed_residues=0,external_STCRDab_still_missing=True,fixed_identity_admissions=0,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False),indent=2)+'\n');print([(x['name'],x['exact_two_method_position_agreement'],x['functional_V_first10_matches'],x['functional_J_last10_matches']) for x in rows])
if __name__=='__main__':main()
