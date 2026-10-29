"""Per-cell Q30 support for every base of the 37 independently assembled pairs."""
import sys,json,collections,argparse,gzip
from pathlib import Path
sys.path[:0]=['.','src']
from scripts.validate_np105_raw import aligned_bases
B=Path('docs/fingerprint_np105_validation_20261029');P=Path('docs/fingerprint_np105_raw_20261029')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--sam',type=Path,default=Path('/home/administrator/docs/np105-raw-20261029/pilot32_validation/observed_reads.sam'));ap.add_argument('--out',type=Path,default=B/'paired_cells_base_support.json');args=ap.parse_args();summary=json.loads((B/'cell_summary.json').read_text());barcodes={x['barcode'] for x in summary['cells'] if x['state']=='both_complete'};numbered=json.loads(gzip.decompress((B/'evidence/all92_numbering.json.gz').read_bytes()));assert barcodes=={c['barcode'] for c in numbered['cells'] if c['state']=='both_complete'};refs=json.loads((P/'raw_sequence_validation.json').read_text())['domains'];evidence=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(set))))
 with args.sam.open() as stream:
  for line in stream:
   if line.startswith('@'):continue
   f=line.rstrip().split('\t');name,flag,role,pos,mapq,cigar=f[:6];flag=int(flag)
   if flag&(4|256|2048) or int(mapq)<30:continue
   bc=name.split('.BC_')[1].split('.UMI_')[0]
   if bc not in barcodes:continue
   umi=name.split('.UMI_')[1];off=refs[role]['ORF_start_nt0'];length=len(refs[role]['nucleotides'])
   for p,base,q in aligned_bases(int(pos),cigar,f[9],f[10]):
    index=p-off
    if 0<=index<length and q>=30 and base in 'ACGT':evidence[bc][role][index][base].add(umi)
 rows=[]
 for bc in sorted(barcodes):
  row=dict(barcode=bc,roles={})
  for role,ref in refs.items():
   counts=[]
   for i,base in enumerate(ref['nucleotides']):
    e=evidence[bc][role][i];other=set().union(*(v for b,v in e.items() if b!=base)) if e else set();counts.append(len(e.get(base,set())-other))
   row['roles'][role]=dict(n_bases=len(counts),minimum_unambiguous_UMI_support=min(counts),n_bases_without_Q30_support=sum(n==0 for n in counts),per_base_unambiguous_UMI_counts=counts)
  row['every_base_in_both_roles_has_Q30_support']=all(x['n_bases_without_Q30_support']==0 for x in row['roles'].values());rows.append(row)
 result=dict(n_independently_assembled_pairs=len(rows),n_pairs_with_Q30_support_at_every_base=sum(x['every_base_in_both_roles_has_Q30_support'] for x in rows),quality_rule='same Q30 bases, MAPQ30 primary alignments and conflicting-UMI exclusion as previous validation; no base substitution',cells=rows,identity_admitted=False)
 args.out.write_text(json.dumps(result,indent=2)+'\n');print({k:result[k] for k in ['n_independently_assembled_pairs','n_pairs_with_Q30_support_at_every_base']})
if __name__=='__main__':main()
