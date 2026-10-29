import sys,json,collections
from pathlib import Path
sys.path.insert(0,'scripts')
from annotate_np105_raw_assembly import annotate,fasta
from Bio.Seq import Seq
root=Path('/home/administrator/docs/np105-raw-20261029');base=Path('docs/fingerprint_np105_raw_20261029')
for out in sorted(root.glob('pilot32_cell_*')):
 if not out.is_dir():continue
 r=annotate(out/'transcripts.fasta');(base/(out.name+'_annotation.json')).write_text(json.dumps(r,indent=2)+'\n')
 print(out.name,[(x['role'],len(x['ORF']),x['CDR3_offset']) for x in r['matches']])
counts=collections.Counter();motifs={'TRA':'TGTGCAGTCATGAATTATAACCAGGGGAAGCTTATCTTT','TRB':'TGTGCCAGCAGCCAAGATAGACGGTCTTCCTATAATTCGCCCCTCTACTTT'}
for mate in [1,2]:
 with (root/'pilot32_selected'/('pooled_R'+str(mate)+'.fastq')).open() as f:
  for i,line in enumerate(f):
   if i%4!=1:continue
   s=line.strip();rev=str(Seq(s).reverse_complement())
   for role,m in motifs.items():
    if m in s or m in rev:counts[(role,mate)]+=1
print('whole_CDR3_nt_matches',dict((str(k),v) for k,v in counts.items()))
partial=[]
for name,nt in fasta(root/'pilot32_pool_assembly_v2/transcripts.fasta'):
 for s in [nt,str(Seq(nt).reverse_complement())]:
  for frame in range(3):
   aa=str(Seq(s[frame:len(s)-(len(s)-frame)%3]).translate())
   if 'NYNQGKL' in aa or 'CAVM' in aa:partial.append(dict(name=name,aa=aa))
(base/'partial_alpha_search.json').write_text(json.dumps(dict(raw_CDR3_read_counts={str(k):v for k,v in counts.items()},assembly_partial_hits=partial),indent=2)+'\n')
print('partial_transcript_hits',len(partial))
