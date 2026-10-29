"""Validate assembled sequence with observed Q30 bases; alignment never fills candidates."""
import os,sys,json,re,collections,subprocess,hashlib,argparse
from pathlib import Path
sys.path[:0]=['.','src','/tmp/np-classic-deps-20261029']
from Bio.Seq import Seq
from scripts.audit_np_adjudication import v_framework_references
from meta_acdc.structure.paired_variable_identity import fasta_records
from scripts.audit_np_boundaries import raw_j_references
B=Path('docs/fingerprint_np105_raw_20261029');RAW=Path('/home/administrator/docs/np105-raw-20261029')

def parse_cigar(cigar):
    ops=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if ''.join(n+op for n,op in ops)!=cigar:raise ValueError('invalid CIGAR')
    return [(int(n),op) for n,op in ops]

def aligned_bases(pos,cigar,seq,qual):
    q=0;r=pos-1
    for n,op in parse_cigar(cigar):
        if op in 'M=X':
            for i in range(n):yield r+i,seq[q+i],ord(qual[q+i])-33
            q+=n;r+=n
        elif op in 'IS':q+=n
        elif op in 'DN':r+=n
    if q!=len(seq):raise ValueError('CIGAR query length differs')

def main():
 global RAW
 ap=argparse.ArgumentParser();ap.add_argument('--raw-root',type=Path,default=RAW);ap.add_argument('--out',type=Path,default=B/'raw_sequence_validation.json');args=ap.parse_args();RAW=args.raw_root
 c=json.loads((B/'raw_classic_numbering.json').read_text());neural=json.loads((B/'raw_neural_numbering.json').read_text())['results'];byrole=collections.defaultdict(list)
 for x,n in zip(c['inputs'],c['numbered']):
  mapping,start,end=n[0];nn=neural[x['name']];assert (start,end)==(nn['query_start'],nn['query_end'])
  aa=x['ORF'][start:end+1];nt=x['ORF_nt'][3*start:3*(end+1)];assert str(Seq(nt).translate())==aa and set(nt)<=set('ACGT')
  assert ''.join(a for p,a in mapping if a!='-')==aa
  byrole[x['role']].append(dict(name=x['name'],barcode=x['name'].split('_')[0],aa=aa,nt=nt,start=start,end=end,ORF=x['ORF'],ORF_nt=x['ORF_nt'],classic=mapping,neural=nn['numbering']))
 assert set(byrole)=={'TRA','TRB'} and all(len(v)==3 for v in byrole.values())
 refs={};domains={}
 for role,rows in byrole.items():
  assert len({x['nt'] for x in rows})==1,'independent cell variable sequences differ'
  x=rows[0];refs[role]=x['ORF_nt'];domains[role]=dict(amino_acids=x['aa'],nucleotides=x['nt'],length_aa=len(x['aa']),ORF_start_nt0=3*x['start'],independent_barcodes=[x['barcode'] for x in rows],same_nt_in_three_cells=True,classic_mapping=x['classic'],neural_mapping=x['neural'],numbering_identical=x['classic']==x['neural'],first_position=x['classic'][0][0],last_non_gap_position=[p for p,a in x['classic'] if a!='-'][-1])
 out=RAW/'pilot32_validation';out.mkdir(exist_ok=True);f=out/'assembled_ORFs.fasta';f.write_text(''.join('>'+role+'\n'+nt+'\n' for role,nt in refs.items()));toolroot=Path(os.environ.get('NP105_SPADES_ROOT','/home/administrator/docs/np105-raw-20261029/spades-tools'));env=dict(os.environ,LD_LIBRARY_PATH=str(toolroot/'usr/lib/x86_64-linux-gnu'));bwa=str(toolroot/'usr/bin/bwa')
 with (out/'bwa.log').open('w') as log:
  subprocess.run([bwa,'index',str(f)],env=env,stderr=log,check=True)
  sam=out/'observed_reads.sam'
  cmd=[bwa,'mem','-t','2',str(f),str(RAW/'pilot32_selected/pooled_R1.fastq'),str(RAW/'pilot32_selected/pooled_R2.fastq')]
  with sam.open('w') as dest:subprocess.run(cmd,env=env,stdout=dest,stderr=log,check=True)
 evidence={role:[collections.defaultdict(set) for _ in nt] for role,nt in refs.items()};alignment=collections.Counter()
 with sam.open() as stream:
  for line in stream:
   if line.startswith('@'):continue
   fields=line.rstrip().split('\t');name,flag,role,pos,mapq,cigar=fields[:6];flag=int(flag)
   if flag&(4|256|2048):continue
   alignment['primary_mapped']+=1
   if int(mapq)<30:alignment['MAPQ_below_30']+=1;continue
   bc=name.split('.BC_')[1].split('.UMI_')[0];umi=name.split('.UMI_')[1];molecule=bc+':'+umi
   for p,base,q in aligned_bases(int(pos),cigar,fields[9],fields[10]):
    if q>=30 and base in 'ACGT':evidence[role][p][base].add(molecule)
 report=[]
 for role,x in domains.items():
  off=x['ORF_start_nt0'];rows=[]
  for i,base in enumerate(x['nucleotides']):
   e=evidence[role][off+i];support=e.get(base,set());other=set().union(*(m for b,m in e.items() if b!=base)) if e else set()
   # A molecule with conflicting calls is explicitly excluded from unambiguous support.
   clean=support-other;cells={m.split(':')[0] for m in clean}
   rows.append(dict(position_nt1=i+1,base=base,Q30_unambiguous_molecules=len(clean),Q30_cells=len(cells),Q30_allele_molecules={b:len(m) for b,m in e.items()},top3_cells_supported=all(bc in cells for bc in x['independent_barcodes'])))
  x['minimum_Q30_molecules']=min(y['Q30_unambiguous_molecules'] for y in rows);x['minimum_Q30_cells']=min(y['Q30_cells'] for y in rows);x['all_bases_supported_in_top3_cells']=all(y['top3_cells_supported'] for y in rows);x['base_evidence']=rows
 facts=json.loads(Path('docs/fingerprint_np_boundaries_20261029/final_primary_cells.json').read_text());np=next(x for x in facts['rows'] if x['name']=='NP105');vrefs=v_framework_references(list(Path('docs/fingerprint_np_boundaries_20261029/sources/imgt').glob('*TR?V.fasta')));jrefs=raw_j_references(list(Path('docs/fingerprint_np_boundaries_20261029/sources/imgt').glob('*TR?J.fasta')))
 for role,x in domains.items():
  primary=np['chains']['alpha' if role=='TRA' else 'beta']['primary_labeled_variable'];x['primary_construct_variable']=primary;x['primary_construct_is_observed_ORF_substring']=primary in str(Seq(refs[role]).translate());
  intervals=[]
  for assembled in byrole[role]:
   begin=assembled['ORF'].find(primary);assert begin>=0 and assembled['ORF'].count(primary)==1
   intervals.append(assembled['ORF_nt'][3*begin:3*(begin+len(primary))])
  assert len(set(intervals))==1
  x['primary_labeled_variable_observed_nt']=intervals[0];x['primary_labeled_variable_identical_nt_in_three_cells']=True
  x['post_numbered_source_residues']=primary[len(x['amino_acids']):]
  x['post_numbered_source_base_evidence']=[]
  for i,base in enumerate(intervals[0][len(x['nucleotides']):],len(x['nucleotides'])):
   e=evidence[role][x['ORF_start_nt0']+i];support=e.get(base,set());other=set().union(*(m for b,m in e.items() if b!=base)) if e else set();clean=support-other;cells={m.split(':')[0] for m in clean}
   x['post_numbered_source_base_evidence'].append(dict(position_relative_to_source_variable_nt1=i+1,base=base,Q30_unambiguous_molecules=len(clean),Q30_cells=len(cells),top3_cells_supported=all(bc in cells for bc in x['independent_barcodes'])))
  x['exact_V_framework_protein_matches']=sorted(v['gene'] for v in vrefs if v['gene'].startswith('TRAV' if role=='TRA' else 'TRBV') and x['amino_acids'].startswith(v['framework_through_conserved_Cys']))
  scores=[]
  for j in jrefs:
   if not j['gene'].startswith('TRAJ' if role=='TRA' else 'TRBJ'):continue
   n=0
   for a,b in zip(reversed(x['amino_acids']),reversed(j['sequence'])):
    if a!=b:break
    n+=1
   scores.append(dict(gene=j['gene'],function=j['function'],exact_suffix_aa=n,reference=j['sequence']))
  mx=max(s['exact_suffix_aa'] for s in scores);x['longest_exact_J_suffix_matches']=[s for s in scores if s['exact_suffix_aa']==mx]
  x['exact_V_framework_nucleotide_matches']=[];x['observed_V_prefix_matches']=[]
  vlen=next((len(v['framework_through_conserved_Cys']) for v in vrefs if v['gene'] in x['exact_V_framework_protein_matches']),None)
  for path in Path('docs/fingerprint_np_boundaries_20261029/sources/imgt').glob('*TR?V.fasta'):
   for header,nt in fasta_records(path.read_text()):
    h=header.split('|')
    if len(h)<8 or h[1] not in x['exact_V_framework_protein_matches'] or h[7] not in ['1','2','3']:continue
    coding=nt.upper()[int(h[7])-1:].replace('.','')
    if vlen and coding[:vlen*3]==x['nucleotides'][:vlen*3]:
     x['exact_V_framework_nucleotide_matches'].append(h[1]);n=0
     for observed,reference in zip(x['nucleotides'],coding):
      if observed!=reference:break
      n+=1
     x['observed_V_prefix_matches'].append(dict(gene=h[1],exact_prefix_nt=n,reference_coding_nt_length=len(coding),domain_interval_nt0=[0,n],biological_recombination_cut_not_resolved=True))
  x['V_framework_length_aa']=vlen
  x['J_reference_nt_exact_substrings']=[]
  for path in Path('docs/fingerprint_np_boundaries_20261029/sources/imgt').glob('*TR?J.fasta'):
   for header,nt in fasta_records(path.read_text()):
    h=header.split('|');nt=nt.upper().replace('.','')
    if len(h)>3 and h[1] in {j['gene'] for j in x['longest_exact_J_suffix_matches']}:
     x['J_reference_nt_exact_substrings'].append(dict(gene=h[1],function=h[3],nt_length=len(nt),full_reference_nt_found_in_assembled_ORF=nt in refs[role],reference_nt=nt,unique_occurrences=refs[role].count(nt),ORF_interval_nt0=([refs[role].find(nt),refs[role].find(nt)+len(nt)] if nt in refs[role] else None),domain_interval_nt0=([refs[role].find(nt)-x['ORF_start_nt0'],refs[role].find(nt)+len(nt)-x['ORF_start_nt0']] if nt in refs[role] else None)))
  cm={tuple(p):a for p,a in x['classic_mapping'] if a!='-'};nm={tuple(p):a for p,a in x['neural_mapping'] if a!='-'}
  x['non_gap_numbering_identical']=cm==nm
  x['numbering_disputes']=[dict(position=list(p),classic=cm.get(p),neural=nm.get(p)) for p in sorted(cm.keys()|nm.keys()) if cm.get(p)!=nm.get(p)]
 result=dict(domains=domains,alignment_counts=dict(alignment),alignment_command=cmd,candidate_reference_is_de_novo_output_only=True,biological_reference_used_for_assembly=False,bases_substituted=0,complete_run_download_verified=False,identity_admitted=False,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
 args.out.write_text(json.dumps(result,indent=2)+'\n')
 print({k:{z:v[z] for z in ['length_aa','same_nt_in_three_cells','minimum_Q30_molecules','minimum_Q30_cells','all_bases_supported_in_top3_cells','exact_V_framework_protein_matches','longest_exact_J_suffix_matches','numbering_identical','primary_construct_is_observed_ORF_substring']} for k,v in domains.items()})
if __name__=='__main__':main()
