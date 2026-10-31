"""Unfiltered germline J evidence for audit only; not an admission policy."""
from .paired_variable_identity import fasta_records,_CODONS,AA

def translate_reference(header,sequence):
    fields=header.split('|')
    if len(fields)<8 or fields[7] not in ('1','2','3') or not fields[1].startswith(('TRAJ','TRBJ')):
        raise ValueError('expected framed IMGT TRAJ/TRBJ reference')
    nt=sequence.replace('.','').upper()
    if not nt or set(nt)-set('ACGTNRYWSKMBDHV'):
        raise ValueError('invalid reference nucleotides')
    offset=int(fields[7])-1
    protein=''.join(_CODONS.get(nt[i:i+3],'X') for i in range(offset,len(nt)-2,3))
    return dict(header=header,gene=fields[1],species=fields[2],functionality=fields[3],frame=int(fields[7]),nucleotides=nt,protein=protein,stop_positions0=[i for i,x in enumerate(protein) if x=='*'],ambiguous_positions0=[i for i,x in enumerate(protein) if x=='X'])

def read_j_evidence(paths):
    return [translate_reference(h,s) for p in paths for h,s in fasta_records(p.read_text())]

def suffix_evidence(variable,records,species,prefix):
    if set(variable)-AA:raise ValueError('observed protein must be canonical and stop-free')
    if len(variable)<10:return []
    result=[]
    for r in records:
        protein=r['protein']
        if r['species']!=species or not r['gene'].startswith(prefix) or len(protein)<10 or variable[-10:]!=protein[-10:]:continue
        n=10
        while n<min(len(variable),len(protein)) and variable[-n-1]==protein[-n-1]:n+=1
        result.append(dict(gene=r['gene'],functionality=r['functionality'],header=r['header'],reference_protein=protein,matched_last10=protein[-10:],longest_exact_suffix_length=n,reference_suffix_start0=len(protein)-n,stop_positions0=r['stop_positions0'],ambiguous_positions0=r['ambiguous_positions0'],old_stopfree_filter_would_drop=bool(r['stop_positions0'] or r['ambiguous_positions0']),all_stops_before_matching_suffix=bool(r['stop_positions0']) and all(i<len(protein)-n for i in r['stop_positions0']),identity_admitted=False))
    return result
