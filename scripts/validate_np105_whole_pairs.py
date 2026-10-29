"""Validate whole-run independent assemblies against observed Q30/UMI bases; no sequence repair."""
from pathlib import Path
import argparse,collections,hashlib,json,subprocess,os,sys
sys.path[:0]=['.','src']
from scripts.validate_np105_raw import aligned_bases
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'docs/fingerprint_np105_raw_20261029'
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as stream:
  for data in iter(lambda:stream.read(4*1024*1024),b''):h.update(data)
 return h.hexdigest()
def unique_canonical_pair(cell):
 for role in ['TRA','TRB']:
  full=[x for x in cell[role] if x['complete']]
  if not full or cell[role+'_complete_unique_nt']!=1 or any(not x['exact_previous_nt'] for x in full):return False
 return True

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--selected',type=Path,required=True);ap.add_argument('--numbering',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(exist_ok=False)
 selection=json.loads((a.selected/'selection_receipt.json').read_text());assert selection['complete'] and selection['source_whole_MD5_verified'] and selection['prefix_all186_selected_outputs_byte_exact'];numbered=json.loads(a.numbering.read_text());assert numbered['n_cells']==92;cells={c['barcode']:c for c in numbered['cells']};paired={bc for bc,c in cells.items() if c['state']=='both_complete'}
 refs=json.loads((P/'raw_sequence_validation.json').read_text())['domains'];original=json.loads((P/'raw_classic_numbering.json').read_text())['inputs'];fasta=a.out/'observed_ORFs.fasta';fasta.write_text(''.join('>'+role+'\n'+next(x['ORF_nt'] for x in original if x['role']==role)+'\n' for role in refs));tools=ROOT.parent/'np105-raw-20261029/spades-tools';env=dict(os.environ,LD_LIBRARY_PATH=str(tools/'usr/lib/x86_64-linux-gnu'));bwa=str(tools/'usr/bin/bwa')
 for mate in [1,2]:
  name=f'pooled_R{mate}.fastq';assert digest(a.selected/name)==selection['outputs'][name]['sha256']
 with (a.out/'bwa.log').open('w') as log:
  subprocess.run([bwa,'index',str(fasta)],env=env,stderr=log,check=True);cmd=[bwa,'mem','-t','2',str(fasta),str(a.selected/'pooled_R1.fastq'),str(a.selected/'pooled_R2.fastq')]
  with (a.out/'whole_reads.sam').open('w') as output:subprocess.run(cmd,env=env,stdout=output,stderr=log,check=True)
 split=a.out/'paired_cell_sam';split.mkdir();handles={};headers=[]
 with (a.out/'whole_reads.sam').open() as stream:
  for line in stream:
   if line.startswith('@'):headers.append(line);continue
   name=line.split('\t',1)[0];bc=name.split('.BC_')[1].split('.UMI_')[0]
   if bc not in paired:continue
   if bc not in handles:handles[bc]=(split/(bc+'.sam')).open('w');handles[bc].writelines(headers)
   handles[bc].write(line)
 for f in handles.values():f.close()
 rows=[];global_counts={role:[collections.Counter() for _ in ref['nucleotides']] for role,ref in refs.items()};global_cells={role:[0 for _ in ref['nucleotides']] for role,ref in refs.items()}
 for bc in sorted(paired):
  evidence=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(set)))
  with (split/(bc+'.sam')).open() as stream:
   for line in stream:
    if line.startswith('@'):continue
    f=line.rstrip().split('\t');name,flag,role,pos,mapq,cigar=f[:6];flag=int(flag)
    if flag&(4|256|2048) or int(mapq)<30:continue
    umi=name.split('.UMI_')[1];assert len(umi)==10;offset=refs[role]['ORF_start_nt0']
    for position,base,q in aligned_bases(int(pos),cigar,f[9],f[10]):
     index=position-offset
     if q>=30 and base in 'ACGT' and 0<=index<len(refs[role]['nucleotides']):evidence[role][index][base].add(umi)
  row=dict(barcode=bc,roles={})
  for role,ref in refs.items():
   counts=[]
   for index,base in enumerate(ref['nucleotides']):
    alleles=evidence[role][index];other=set().union(*(v for b,v in alleles.items() if b!=base)) if alleles else set();count=len(alleles.get(base,set())-other);counts.append(count);global_cells[role][index]+=count>0
    for allele,molecules in alleles.items():
     conflicts=set().union(*(v for b,v in alleles.items() if b!=allele));global_counts[role][index][allele]+=len(molecules-conflicts)
   row['roles'][role]=dict(minimum_unambiguous_UMI_support=min(counts),n_bases_without_Q30_support=sum(x==0 for x in counts),per_base_unambiguous_UMI_counts=counts)
  row['all_reference_bases_have_same_cell_Q30_support']=all(x['n_bases_without_Q30_support']==0 for x in row['roles'].values());c=cells[bc];row['single_full_sequence_per_role_and_exact_prior_nt']=unique_canonical_pair(c);row['observed_canonical_pair_supported']=row['all_reference_bases_have_same_cell_Q30_support'] and row['single_full_sequence_per_role_and_exact_prior_nt'];rows.append(row);print('CELL_Q30_DONE',len(rows),'OF',len(paired),flush=True) if len(rows)%12==0 else None
 result=dict(source_whole_MD5_verified=True,n_total_whitelist_cells=92,n_independently_assembled_paired_cells=len(paired),n_paired_cells_all_reference_bases_with_same_cell_Q30_support=sum(x['all_reference_bases_have_same_cell_Q30_support'] for x in rows),n_unique_canonical_pairs_with_every_base_Q30_support=sum(x['observed_canonical_pair_supported'] for x in rows),numbering_cell_counts=numbered['cell_counts'],n_complete_ORFs_different_from_prior=numbered['n_complete_ORFs_different_from_prior'],cells=rows,aggregate_canonical_support_cells_per_base=global_cells,aggregate_unambiguous_base_UMI_counts={role:[dict(x) for x in values] for role,values in global_counts.items()},postassembly_reference_is_observed_ORFs_only=True,biological_reference_used_for_assembly=False,bases_substituted=0,alignment_command=cmd,identity_admitted=False,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
 (a.out/'whole_base_support.json').write_text(json.dumps(result,indent=2)+'\n');print({k:result[k] for k in ['n_independently_assembled_paired_cells','n_paired_cells_all_reference_bases_with_same_cell_Q30_support','n_unique_canonical_pairs_with_every_base_Q30_support','n_complete_ORFs_different_from_prior']})
if __name__=='__main__':main()
