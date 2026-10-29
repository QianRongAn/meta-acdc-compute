"""Trace whole-run assembly variants to observed bases and CDR3-linked UMI support; no repair."""
from pathlib import Path
import argparse,json,collections,hashlib,sys
sys.path[:0]=['.','src']
from scripts.validate_np105_raw import aligned_bases
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'docs/fingerprint_np105_raw_20261029'
CDR3={'TRA':'TGTGCAGTCATGAATTATAACCAGGGGAAGCTTATCTTT','TRB':'TGTGCCAGCAGCCAAGATAGACGGTCTTCCTATAATTCGCCCCTCTACTTT'}
def haplotype_umis(calls,positions,haplotype):
 if not positions or len(positions)!=len(haplotype):return set()
 return {umi for umi,observed in calls.items() if all(observed.get(i,set())=={base} for i,base in zip(positions,haplotype))}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--numbering',type=Path,required=True);ap.add_argument('--sam',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();numbered=json.loads(a.numbering.read_text());refs=json.loads((P/'raw_sequence_validation.json').read_text())['domains'];targets={}
 for v in numbered['variants']:
  ref=refs[v['role']]['nucleotides'];nt=v['nucleotides'];key=(v['barcode'],v['role'],nt)
  if key in targets:continue
  differences=[i for i,(x,y) in enumerate(zip(ref,nt)) if x!=y] if len(ref)==len(nt) else []
  targets[key]=dict(barcode=v['barcode'],role=v['role'],candidate_nt_sha256=hashlib.sha256(nt.encode()).hexdigest(),candidate_nt=nt,reference_nt=ref,positions=differences,scope='whole-run de novo assembly candidate',equal_length_snv_audit_supported=bool(differences) and len(ref)==len(nt))
 old=json.loads((ROOT/'docs/fingerprint_np105_validation_20261029/variant_raw_support.json').read_text())['targets'][0];ref=refs['TRA']['nucleotides'];nt=ref[:200]+'C'+ref[201:];targets[(old['barcode'],'TRA',nt)]=dict(barcode=old['barcode'],role='TRA',candidate_nt_sha256=hashlib.sha256(nt.encode()).hexdigest(),candidate_nt=nt,reference_nt=ref,positions=[200],scope='predeclared prefix K82N candidate; check against whole run',equal_length_snv_audit_supported=True)
 for t in targets.values():t['_calls']=collections.defaultdict(lambda:collections.defaultdict(set));t['_cdr3']=set();t['_phased']=collections.defaultdict(set);t['_all_q_reads']=collections.Counter()
 bybar=collections.defaultdict(list)
 for t in targets.values():bybar[(t['barcode'],t['role'])].append(t)
 with a.sam.open() as stream:
  for line in stream:
   if line.startswith('@'):continue
   f=line.rstrip().split('\t');name,flag,role,pos,mapq,cigar=f[:6];flag=int(flag)
   if flag&(4|256|2048) or int(mapq)<30:continue
   bc=name.split('.BC_')[1].split('.UMI_')[0];ts=bybar.get((bc,role))
   if not ts:continue
   umi=name.split('.UMI_')[1];assert len(umi)==10;bases={p:(base,q) for p,base,q in aligned_bases(int(pos),cigar,f[9],f[10])};off=refs[role]['ORF_start_nt0'];cdr3start=off+refs[role]['nucleotides'].index(CDR3[role]);linked=all(cdr3start+i in bases and bases[cdr3start+i][0]==base and bases[cdr3start+i][1]>=30 for i,base in enumerate(CDR3[role]))
   for t in ts:
    if linked:t['_cdr3'].add(umi)
    positions=t['positions']
    for i in positions:
     if off+i in bases and bases[off+i][1]>=30:t['_calls'][umi][i].add(bases[off+i][0])
    if positions and all(off+i in bases for i in positions):
     hap=''.join(bases[off+i][0] for i in positions);t['_all_q_reads'][hap]+=1
     if all(bases[off+i][1]>=30 for i in positions):t['_phased'][hap].add(umi)
 rows=[]
 for t in targets.values():
  positions=t['positions'];variant=''.join(t['candidate_nt'][i] for i in positions);canonical=''.join(t['reference_nt'][i] for i in positions);support={}
  for label,hap in [('candidate',variant),('canonical',canonical)]:
   UMIs=haplotype_umis(t['_calls'],positions,hap);support[label]=dict(haplotype=hap,Q30_unambiguous_UMIs_covering_all_differences=len(UMIs),Q30_exact_CDR3_linked_unambiguous_UMIs=len(UMIs&t['_cdr3']),Q30_same_read_phased_UMIs=len(t['_phased'].get(hap,set())))
  row={k:v for k,v in t.items() if not k.startswith('_') and k not in ['candidate_nt','reference_nt']};row.update(positions_nt1=[i+1 for i in positions],support=support,all_quality_same_read_haplotype_read_counts=dict(t['_all_q_reads']),Q30_same_read_haplotype_UMI_counts={h:len(v) for h,v in t['_phased'].items()},interpretation='raw support audit only; does not establish an independent native receptor, genomic phase, or mutation mechanism');rows.append(row)
 a.out.write_text(json.dumps(dict(scope='all selected reads from checksum-verified complete SRR34543498',targets=rows,reference_is_observed_prefix_candidate_only=True,biological_reference_used_for_assembly=False,bases_substituted=0,identity_admitted=False),indent=2)+'\n')
 for x in rows:print(x['barcode'],x['role'],x['positions_nt1'],x['support'])
if __name__=='__main__':main()
