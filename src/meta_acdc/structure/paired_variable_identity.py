"""Conservative deposited paired-variable identity; never impute missing residues."""
import hashlib
import html
import itertools
import re
import shlex
from collections import defaultdict
from .graph import AA3TO1

AA=set('ACDEFGHIKLMNPQRSTVWY')
_CODONS={''.join(c):a for c,a in zip(itertools.product('TCAG',repeat=3),'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}

def fasta_records(text):
    rows=[]
    for block in text.split('>')[1:]:
        lines=block.splitlines()
        rows.append((lines[0].strip(),''.join(x.strip() for x in lines[1:])))
    return rows

def germline_references(paths):
    result=defaultdict(list)
    for path in paths:
        for header,seq in fasta_records(path.read_text()):
            p=header.split('|')
            if len(p)<8 or p[3] not in ('F','(F)') or p[7] not in ('1','2','3'):continue
            nt=seq.replace('.','').upper()[int(p[7])-1:]
            aa=''.join(_CODONS.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))
            if '*' in aa or 'X' in aa:continue
            species='human' if p[2]=='Homo sapiens' else 'mouse' if p[2].startswith('Mus musculus') else None
            if species:result[(species,p[1])].append(dict(amino_acids=aa,header=header))
    return result

def numbered_domains(text):
    out={}
    # Each outer chain row contains one inner numbering table. Fail on conflicts.
    pattern=r'<td>(VA|VB)</td>\s*<td>([^<]+)</td>.*?<table[^>]*>(.*?)</table>'
    for kind,chain,body in re.findall(pattern,text,re.S):
        rows=re.findall(r'<tr>(.*?)</tr>',body,re.S);positions=[];letters=[]
        for i in range(0,len(rows),2):
            if i+1==len(rows):raise ValueError('odd numbered rows')
            strip=lambda row:[html.unescape(re.sub('<[^>]*>','',x)).strip() for x in re.findall(r'<td[^>]*>(.*?)</td>',row,re.S)]
            labels,aas=strip(rows[i]),strip(rows[i+1])
            if len(labels)!=len(aas):raise ValueError('numbered row length mismatch')
            pairs=[(p,a) for p,a in zip(labels,aas) if p or a]
            labels=[p for p,a in pairs];aas=[a for p,a in pairs]
            if len(labels)!=len(aas) or not all(re.fullmatch('[AB][0-9]+[A-Z]?',x) for x in labels) or not all(x in AA for x in aas):raise ValueError('invalid numbered domain')
            positions.extend(labels);letters.extend(aas)
        if len(positions)!=len(set(positions)):raise ValueError('repeated IMGT position')
        domain=dict(kind=kind,positions=positions,sequence=''.join(letters))
        if chain in out and out[chain]!=domain:raise ValueError('conflicting domain annotation')
        out[chain]=domain
    return out

def validate_domain(seqres,annotation,fasta_text,html_text,refs):
    rows=fasta_records(fasta_text)
    if not rows:raise ValueError('missing annotation FASTA')
    header,variable=rows[0];parts=header.split('|');kind=annotation['kind']
    if len(parts)!=3 or parts[1]!=kind or variable!=annotation['sequence']:raise ValueError('FASTA/numbering mismatch')
    if len(variable)<80 or len(variable)>140 or set(variable)-AA:raise ValueError('ambiguous variable sequence')
    if seqres.count(variable)!=1:raise ValueError('variable not unique contiguous deposited sequence')
    gene=parts[2];species_matches=re.findall(r'<b>Species</b></td><td[^>]*>([^<]+)</td>',html_text)
    species=set(x.strip() for x in species_matches)
    if len(species)!=1:raise ValueError('species annotation absent or ambiguous')
    species=next(iter(species));vrefs=refs.get((species,gene),[])
    if not any(variable.startswith(x['amino_acids'][:10]) for x in vrefs):raise ValueError('full germline N terminus unverified; no imputation or tag trimming')
    # J allele need not be called: complete FR4 must exactly match a functional J reference.
    prefix='TRAJ' if kind=='VA' else 'TRBJ'
    jmatches=[(g,x) for (sp,g),values in refs.items() if sp==species and g.startswith(prefix) for x in values if variable.endswith(x['amino_acids'][-10:])]
    if not jmatches:raise ValueError('complete J/FR4 C terminus unverified')
    positions=annotation['positions'];numbers=[int(re.search('[0-9]+',p).group()) for p in positions]
    if min(numbers)>3 or max(numbers)<(128 if kind=='VA' else 127):raise ValueError('terminal numbering incomplete')
    for _,loop in rows[1:]:
        if not loop or loop not in variable:raise ValueError('CDR annotation conflict')
    start=seqres.index(variable)
    return dict(variable=variable,v_allele=gene,species=species,variable_start0=start,variable_end0=start+len(variable),n_terminal_extension=seqres[:start],c_terminal_extension=seqres[start+len(variable):],numbered_positions=positions,
                v_reference_headers=sorted({x['header'] for x in vrefs if variable.startswith(x['amino_acids'][:10])}),j_fr4_reference_headers=sorted({x['header'] for _,x in jmatches}),boundary_basis='exact STCRDab numbered SEQRES substring; IMGT functional V first10 and J final10 exact; engineered internal substitutions retained')

def pair_key(alpha,beta):
    if not alpha or not beta or (set(alpha)|set(beta))-AA:raise ValueError('noncanonical/empty paired sequence')
    return hashlib.sha256((alpha+'|'+beta).encode()).hexdigest()

def cif_declared_sequences(text):
    lines=text.splitlines();entities=defaultdict(dict);chains={}
    for i,line in enumerate(lines):
        if line.strip()!='loop_':continue
        j=i+1;header=[]
        while j<len(lines) and lines[j].startswith('_'):header.append(lines[j].strip());j+=1
        if not header or not any(h.startswith(('_entity_poly.','_entity_poly_seq.')) for h in header):continue
        while j<len(lines) and not lines[j].startswith(('#','_','loop_','data_')):
            if lines[j].strip():
                row=shlex.split(lines[j])
                if len(row)!=len(header):raise ValueError('unsupported multline CIF declaration')
                fields=dict(zip(header,row))
                if '_entity_poly.entity_id' in fields:
                    if fields.get('_entity_poly.type')!='polypeptide(L)':raise ValueError('unsupported polymer type')
                    for chain in fields['_entity_poly.pdbx_strand_id'].split(','):
                        if chain in chains:raise ValueError('duplicate declared chain')
                        chains[chain]=fields['_entity_poly.entity_id']
                if '_entity_poly_seq.entity_id' in fields:
                    e=fields['_entity_poly_seq.entity_id'];n=int(fields['_entity_poly_seq.num']);a=AA3TO1.get(fields['_entity_poly_seq.mon_id'],'X')
                    if n in entities[e] or a not in AA:raise ValueError('ambiguous polymer sequence')
                    entities[e][n]=a
            j+=1
    out={}
    for chain,entity in chains.items():
        residues=entities[entity]
        if not residues or sorted(residues)!=list(range(1,len(residues)+1)):raise ValueError('incomplete declared polymer sequence')
        out[chain]=''.join(residues[i] for i in range(1,len(residues)+1))
    if not out:raise ValueError('no declared polymer sequence')
    return out
