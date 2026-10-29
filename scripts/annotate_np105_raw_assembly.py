"""Post-assembly CDR3 annotation only; never feeds sequences into assembly."""
import sys,json,argparse,hashlib
from pathlib import Path
sys.path.insert(0,"/tmp/np-classic-deps-20261029")
from Bio.Seq import Seq

def fasta(path):
    entries=[]
    for line in path.read_text().splitlines():
        if line.startswith('>'):entries.append([line[1:],''])
        else:entries[-1][1]+=line.strip()
    return entries

def annotate(path):
    rows=[]
    for name,nt in fasta(path):
        for strand,dna in [(1,nt),(-1,str(Seq(nt).reverse_complement()))]:
            for frame in range(3):
                aa=str(Seq(dna[frame:len(dna)-(len(dna)-frame)%3]).translate())
                for role,motif in [('TRA','CAVMNYNQGKLIF'),('TRB','CASSQDRRSSYNSPLYF')]:
                    at=aa.find(motif)
                    if at<0:continue
                    start=aa.rfind('*',0,at)+1;end=aa.find('*',at)
                    if end<0:end=len(aa)
                    rows.append(dict(transcript=name,nt_length=len(nt),strand=strand,frame=frame,role=role,ORF=aa[start:end],ORF_nt=dna[frame+3*start:frame+3*end],CDR3_offset=at-start,ORF_5prime_stop_bounded=start>0,ORF_3prime_stop_bounded=end<len(aa),N_count=nt.count('N')))
    return dict(input=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),transcripts=len(fasta(path)),annotation_only=True,biological_reference_used_for_assembly=False,matches=rows)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('fasta',type=Path);ap.add_argument('out',type=Path);a=ap.parse_args();r=annotate(a.fasta);a.out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(transcripts=r['transcripts'],matches=[dict(role=x['role'],transcript=x['transcript'],aa_len=len(x['ORF']),before_CDR3=x['CDR3_offset'],N_count=x['N_count']) for x in r['matches']]),indent=2))
