"""Audit all observed Q30 codons at the assembled alpha variant; no consensus replacement."""
from pathlib import Path
import json,collections,sys,argparse
sys.path[:0]=['.','src']
from scripts.validate_np105_raw import aligned_bases
B=Path('docs/fingerprint_np105_validation_20261029');P=Path('docs/fingerprint_np105_raw_20261029');RAW=Path('/home/administrator/docs/np105-raw-20261029')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--numbering',type=Path,default=B/'all92_numbering.json');ap.add_argument('--sam',type=Path,default=RAW/'pilot32_validation/observed_reads.sam');ap.add_argument('--out',type=Path,default=B/'variant_raw_support.json');a=ap.parse_args();result=json.loads(a.numbering.read_text());prior=json.loads((P/'raw_sequence_validation.json').read_text());targets=[]
 for variant in result['variants']:
  ref=prior['domains'][variant['role']];assert len(variant['nucleotides'])==len(ref['nucleotides'])
  positions=[i for i,(a,b) in enumerate(zip(ref['nucleotides'],variant['nucleotides'])) if a!=b];assert len(positions)==1
  pos=positions[0];codon_start=pos//3*3;aa_index=pos//3;mapping=ref['classic_mapping'];observed=[(p,a) for p,a in mapping if a!='-'];imgt=observed[aa_index][0]
  targets.append(dict(barcode=variant['barcode'],role=variant['role'],position_domain_nt1=pos+1,codon_start_domain_nt0=codon_start,position_IMGT=imgt,reference_nt=ref['nucleotides'][pos],variant_nt=variant['nucleotides'][pos],reference_codon=ref['nucleotides'][codon_start:codon_start+3],variant_codon=variant['nucleotides'][codon_start:codon_start+3],reference_aa=ref['amino_acids'][aa_index],variant_aa=variant['boundary']['sequence'][aa_index],ORF_codon_start_nt0=ref['ORF_start_nt0']+codon_start))
 for t in targets:
  observations=collections.defaultdict(list);raw_counts=collections.Counter();all_quality_target=collections.Counter();all_quality_variants=[];linked=collections.defaultdict(set);cell_molecules=collections.defaultdict(lambda:collections.defaultdict(set))
  with a.sam.open() as stream:
   for line in stream:
    if line.startswith('@'):continue
    f=line.rstrip().split('\t');name,flag,role,pos,mapq,cigar=f[:6];flag=int(flag)
    if role!=t['role'] or flag&(4|256|2048) or int(mapq)<30:continue
    bases={p:(base,q) for p,base,q in aligned_bases(int(pos),cigar,f[9],f[10])};start=t['ORF_codon_start_nt0']
    if not all(start+i in bases for i in range(3)):continue
    raw_codon=''.join(bases[start+i][0] for i in range(3));bc=name.split('.BC_')[1].split('.UMI_')[0];umi=name.split('.UMI_')[1]
    if bc==t['barcode']:
     all_quality_target[raw_codon]+=1
     if raw_codon==t['variant_codon']:all_quality_variants.append(dict(spot=name.split('.BC_')[0],UMI=umi,mate=1 if flag&64 else 2,codon=raw_codon,qualities=[bases[start+i][1] for i in range(3)]))
    if not all(bases[start+i][1]>=30 for i in range(3)):continue
    cdr3='TGTGCAGTCATGAATTATAACCAGGGGAAGCTTATCTTT';ref=prior['domains'][t['role']];offset=ref['ORF_start_nt0']+ref['nucleotides'].index(cdr3)
    if all(offset+i in bases and bases[offset+i][0]==base and bases[offset+i][1]>=30 for i,base in enumerate(cdr3)):linked[raw_codon].add(bc+':'+umi)
    codon=''.join(bases[start+i][0] for i in range(3));bc=name.split('.BC_')[1].split('.UMI_')[0];umi=name.split('.UMI_')[1];cell_molecules[bc][umi].add(codon)
    if bc==t['barcode']:raw_counts[codon]+=1;observations[umi].append(dict(spot=name.split('.BC_')[0],mate=1 if flag&64 else 2,codon=codon,min_codon_Q=min(bases[start+i][1] for i in range(3)),MAPQ=int(mapq)))
  unambiguous=collections.Counter();conflicting=0;cell_counts=collections.Counter()
  for bc,molecules in cell_molecules.items():
   for umi,codons in molecules.items():
    if len(codons)==1:
     codon=next(iter(codons));cell_counts[codon]+=1
     if bc==t['barcode']:unambiguous[codon]+=1
    elif bc==t['barcode']:conflicting+=1
  t.update(target_cell_all_quality_read_codon_counts=dict(all_quality_target),variant_codon_target_read_qualities=all_quality_variants,Q30_same_read_exact_CDR3_linked_molecules={codon:len(m) for codon,m in linked.items()},target_cell_Q30_read_codon_counts=dict(raw_counts),target_cell_unambiguous_UMI_codon_counts=dict(unambiguous),target_cell_conflicting_UMIs=conflicting,all_whitelist_unambiguous_barcode_UMI_codon_counts=dict(cell_counts),all_whitelist_cells_by_codon={codon:sum(any(v=={codon} for v in m.values()) for m in cell_molecules.values()) for codon in cell_counts},target_cell_observations_by_UMI=dict(observations),decision='quarantine assembly variant; do not replace consensus or infer an independent natural receptor')
 a.out.write_text(json.dumps(dict(scope='bounded prefix; post-assembly candidate-aligned raw observed codons, MAPQ>=30 and Q>=30 on all three bases',targets=targets,biological_reference_used_for_assembly=False),indent=2)+'\n')
 for t in targets:print({k:v for k,v in t.items() if k not in ['target_cell_observations_by_UMI']})
if __name__=='__main__':main()
