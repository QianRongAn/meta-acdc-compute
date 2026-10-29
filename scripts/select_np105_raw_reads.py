"""Select raw NP105 cell reads without any biological reference sequence or base substitution."""
from pathlib import Path
import argparse,collections,hashlib,json,itertools
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np105_raw_20261029'
def records(path):
    with path.open('rt',encoding='ascii') as f:
        while True:
            h=f.readline()
            if not h:return
            seq=f.readline();plus=f.readline();qual=f.readline()
            if not all([seq,plus,qual]) or not qual.endswith('\n'):return
            h,seq,plus,qual=[x.rstrip('\r\n') for x in [h,seq,plus,qual]]
            if not h.startswith('@') or not plus.startswith('+') or len(seq)!=len(qual):raise ValueError('malformed complete FASTQ record')
            yield h,seq,qual

def trim_ends(seq,qual,q=20):
    left=0;right=len(seq)
    while left<right and ord(qual[left])-33<q:left+=1
    while right>left and ord(qual[right-1])-33<q:right-=1
    return seq[left:right],qual[left:right]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--r1',type=Path,required=True);ap.add_argument('--r2',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    wl=json.loads((BASE/'barcode_whitelist.json').read_text());barcodes={b[:-2] for b in wl['barcodes']};assert len(barcodes)==92
    handles={};counts=collections.Counter();umis=collections.defaultdict(set);pooled=[(a.out/'pooled_R1.fastq').open('xb'),(a.out/'pooled_R2.fastq').open('xb')];side=(a.out/'read_provenance.tsv').open('w');side.write('spot\tbarcode\tUMI\n');n_pairs=0;selected=0;reasons=collections.Counter()
    source_counts=[0,0];unpaired_suffix=0
    for one,two in itertools.zip_longest(records(a.r1),records(a.r2)):
        source_counts[0]+=one is not None;source_counts[1]+=two is not None
        if one is None or two is None:unpaired_suffix+=1;continue
        n_pairs+=1;h1,s1,q1=one;h2,s2,q2=two;spot=h1.split()[0]
        if spot!=h2.split()[0]:raise ValueError('mate identifiers disagree')
        bc=s1[:16]
        if bc not in barcodes:continue
        if len(s1)<86 or len(s2)<60:reasons['short_input']+=1;continue
        if min(ord(c)-33 for c in q1[:26])<20:reasons['barcode_or_UMI_quality_below_20']+=1;continue
        umi=s1[16:26]
        if set(umi)-set('ACGT'):reasons['noncanonical_UMI']+=1;continue
        # Only the observed technical barcode/UMI prefix is removed; all biological bases remain observed bases.
        s1,q1=trim_ends(s1[26:],q1[26:]);s2,q2=trim_ends(s2,q2)
        if min(len(s1),len(s2))<60:reasons['short_after_Q20_end_trim']+=1;continue
        if set(s1+s2)-set('ACGTN'):raise ValueError('invalid raw DNA base')
        if bc not in handles:handles[bc]=[(a.out/(bc+'_R1.fastq')).open('xb'),(a.out/(bc+'_R2.fastq')).open('xb')]
        for mate,(seq,qual) in enumerate([(s1,q1),(s2,q2)],1):
            record=(spot+'.BC_'+bc+'.UMI_'+umi+'/'+str(mate)+'\n'+seq+'\n+\n'+qual+'\n').encode();pooled[mate-1].write(record);handles[bc][mate-1].write(record)
        side.write(spot[1:]+'\t'+bc+'\t'+umi+'\n');selected+=1;counts[bc]+=1;umis[bc].add(umi)
    for f in pooled+[f for hs in handles.values() for f in hs]:f.close()
    side.close();top=sorted(counts,key=lambda b:(-counts[b],b))[:3]
    result=dict(n_complete_records_per_source=source_counts,n_unpaired_prefix_suffix_records=unpaired_suffix,input_paths=[str(a.r1),str(a.r2)],n_complete_synchronized_prefix_pairs=n_pairs,n_selected_pairs=selected,n_selected_barcodes=len(counts),read_pairs_by_barcode=dict(counts),UMIs_by_barcode={b:len(v) for b,v in umis.items()},top_three_by_preassembly_read_count=top,excluded_reasons=dict(reasons),technical_bases_removed_from_R1=26,bases_substituted=0,reference_sequences_used=False,source_scope='complete paired FASTQ records from a bounded gzip prefix; not a complete run or whole-file MD5 verification',outputs={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in a.out.glob('*.fastq')})
    (a.out/'selection_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['n_complete_synchronized_prefix_pairs','n_selected_pairs','n_selected_barcodes','top_three_by_preassembly_read_count','excluded_reasons']},indent=2))
if __name__=='__main__':main()
