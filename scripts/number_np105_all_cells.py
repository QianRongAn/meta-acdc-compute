"""Number every post-assembly CDR3-bearing ORF, retaining all failures and variants."""
import sys,json,collections,hashlib,argparse
from pathlib import Path
sys.path[:0]=['.','src','/tmp/np-classic-deps-20261029']
import anarci
from scripts.audit_np_boundaries import checked_numbering
B=Path('docs/fingerprint_np105_validation_20261029');P=Path('docs/fingerprint_np105_raw_20261029')
def strict_boundary(full,numbering,start,end,detail,role):
 expected='A' if role=='TRA' else 'B';boundary=checked_numbering(full,numbering,start,end,detail['chain_type'],expected)
 if boundary['first_IMGT']!='1' or boundary['last_IMGT']!=('128' if expected=='A' else '127'):raise ValueError('strict full interval endpoints not observed')
 if detail.get('species')!='mouse':raise ValueError('numbering fallback species is not mouse')
 return boundary

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--assemblies',type=Path,default=B/'all92_assemblies.json');ap.add_argument('--out',type=Path,default=B/'all92_numbering.json');args=ap.parse_args()
 assembled=json.loads(args.assemblies.read_text());assert assembled['finished'] and len(assembled['cells'])==92
 inputs=[]
 for row in assembled['cells']:
  for i,hit in enumerate(row.get('annotation',{}).get('matches',[])):
   inputs.append(dict(name=row['barcode']+'_'+str(i),barcode=row['barcode'],**hit))
 numbered,details,hits=anarci.anarci([(x['name'],x['ORF']) for x in inputs],scheme='imgt',allow={'A','B'},hmmerpath='/tmp/np-hmmer-bin-20261029/usr/bin',ncpu=2,assign_germline=True,allowed_species=['mouse'],bit_score_threshold=80)
 prior=json.loads((P/'raw_sequence_validation.json').read_text());cells={x['barcode']:dict(barcode=x['barcode'],read_pairs=x['selected_read_pairs'],assembly_exit_code=x['exit_code'],TRA=[],TRB=[]) for x in assembled['cells']};annotations=[]
 for x,nn,dd,hh in zip(inputs,numbered,details,hits):
  result=dict(name=x['name'],barcode=x['barcode'],role=x['role'],transcript=x['transcript'],ORF=x['ORF'],ORF_nt=x['ORF_nt'],numbered=nn,details=dd,hits=hh,complete=False)
  try:
   if not nn or len(nn)!=1 or len(dd)!=1:raise ValueError('missing or multiple numbered domains')
   n,start,end=nn[0];expected='A' if x['role']=='TRA' else 'B';boundary=strict_boundary(x['ORF'],n,start,end,dd[0],x['role']);nt=x['ORF_nt'][start*3:(end+1)*3];ref=prior['domains'][x['role']]
   assert len(nt)==3*boundary['length'] and set(nt)<=set('ACGT')
   result.update(complete=True,boundary=boundary,nucleotides=nt,exact_previous_aa=boundary['sequence']==ref['amino_acids'],exact_previous_nt=nt==ref['nucleotides'],sha256_nt=hashlib.sha256(nt.encode()).hexdigest())
  except (ValueError,AssertionError,KeyError) as e:result['unresolved_reason']=str(e) or type(e).__name__
  annotations.append(result);cells[x['barcode']][x['role']].append(result)
 stats=collections.Counter();variants=[]
 for row in cells.values():
  for role in ['TRA','TRB']:
   full=[x for x in row[role] if x['complete']];row[role+'_complete_count']=len(full);row[role+'_complete_unique_nt']=len({x['nucleotides'] for x in full})
   variants.extend(x for x in full if not x['exact_previous_nt'])
  a=row['TRA_complete_count']>0;b=row['TRB_complete_count']>0;state='both_complete' if a and b else 'alpha_only_complete' if a else 'beta_only_complete' if b else 'neither_complete';row['state']=state;stats[state]+=1
 result=dict(all_cells_predeclared=True,n_cells=len(cells),n_CDR3_bearing_ORFs=len(inputs),n_complete_ORFs=sum(x['complete'] for x in annotations),cell_counts=dict(stats),n_complete_ORFs_different_from_prior=len(variants),variants=variants,cells=list(cells.values()),input_sha256=hashlib.sha256(args.assemblies.read_bytes()).hexdigest(),biological_reference_used_for_assembly=False,all_cells_are_one_receptor_source=True,new_independent_receptors=0,identity_admitted=False)
 args.out.write_text(json.dumps(result,indent=2)+'\n');print({k:result[k] for k in ['n_cells','n_CDR3_bearing_ORFs','n_complete_ORFs','cell_counts','n_complete_ORFs_different_from_prior']})
if __name__=='__main__':main()
