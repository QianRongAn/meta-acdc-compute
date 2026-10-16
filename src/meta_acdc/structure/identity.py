"""Depositor COMPND annotations + contact-checked class-I alpha/beta units.
Unknown/fusion/class-II records are explicitly excluded from this release scope.
No chain letter or length is used to guess an alpha/beta receptor identity.
"""
from __future__ import annotations
import re
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree


def pdb_entities(path: Path):
    text=[]
    with path.open() as f:
        for line in f:
            if line.startswith('COMPND'):text.append(line[10:].strip())
            if line.startswith(('ATOM  ','HETATM')):break
    entities=[];current={}
    for token in ' '.join(text).split(';'):
        if ':' not in token:continue
        k,v=token.split(':',1);k=k.strip();v=v.strip()
        if k=='MOL_ID':
            if current:entities.append(current)
            current={'mol_id':v}
        else:current[k.lower()]=v
    if current:entities.append(current)
    result={}
    for ent in entities:
        desc=(ent.get('molecule','')+' '+ent.get('synonym','')).upper().replace('Α','ALPHA').replace('Β','BETA')
        kind='unknown';role='other'
        tcr=any(x in desc for x in ('TCR','T-CELL','T CELL','T-LYMPHOCYTE','T LYMPHOCYTE','ANTIGEN RECEPTOR'))
        if 'BETA-2' in desc or 'BETA 2' in desc or 'MICROGLOBULIN' in desc:
            kind='b2m';role='b2m'
        elif tcr:
            a='ALPHA' in desc or 'LIGHT CHAIN' in desc or 'TRAC' in desc
            b='BETA' in desc or 'HEAVY CHAIN' in desc or 'TRBC' in desc
            if a and not b:kind='tra';role='tcr'
            elif b and not a:kind='trb';role='tcr'
        elif 'CLASS II' in desc or re.search(r'HLA[- ]?D[PRQ]',desc):kind='mhc_ii'
        elif 'CLASS I ' in desc or re.search(r'HLA[- ]?[ABCEFG](?:\b|[ *0-9])',desc):
            kind='mhc_i';role='mhc'
        for chain in ent.get('chain','').split(','):
            if chain.strip():result[chain.strip()]={'kind':kind,'role':role,'description':ent.get('molecule',''),'mol_id':ent['mol_id'],'basis':'PDB_COMPND'}
    return result


def annotate(residues, path):
    annotations=pdb_entities(path);chains=defaultdict(list)
    for res in residues:chains[res.chain].append(res)
    for chain,rs in chains.items():
        ent=annotations.get(chain,{'kind':'unknown','role':'other','basis':'missing_COMPND','description':''})
        # A short ATOM polymer may be a peptide, never an alpha/beta guess.
        if ent['kind']=='unknown' and 5<=len(rs)<=20:
            ent=dict(ent,kind='peptide',role='peptide',basis=ent['basis']+'+short_polymer')
        for res in rs:
            res.role=ent['role'];res.chain_kind=ent['kind']
    return annotations


def chain_sequence(rs):
    from .graph import AA3TO1
    return ''.join(AA3TO1.get(r.resname.upper(),'X') for r in sorted(rs,key=lambda x:(x.resid,x.icode)))


def contact_count(a,b,cutoff=6.0):
    def atoms(rs):
        xyz=[];indices=[]
        for i,res in enumerate(rs):
            for name,x,y,z in res.atoms:
                if not name.startswith('H'):xyz.append((x,y,z));indices.append(i)
        return np.asarray(xyz),indices
    aa,ai=atoms(a);bb,bi=atoms(b)
    if len(aa)==0 or len(bb)==0:return 0
    hits=cKDTree(aa).query_ball_tree(cKDTree(bb),cutoff)
    left=set();right=set()
    for i,js in enumerate(hits):
        if js:left.add(ai[i]);right.update(bi[j] for j in js)
    return min(len(left),len(right))


def class_i_units(residues):
    chains=defaultdict(list)
    for r in residues:chains[r.chain].append(r)
    kinds={c:rs[0].chain_kind for c,rs in chains.items()}
    al=sorted(c for c,k in kinds.items() if k=='tra');be=sorted(c for c,k in kinds.items() if k=='trb')
    pe=sorted(c for c,k in kinds.items() if k=='peptide');mh=sorted(c for c,k in kinds.items() if k=='mhc_i');b2=sorted(c for c,k in kinds.items() if k=='b2m')
    excluded=[];units=[]
    if not all([al,be,pe,mh,b2]):return [],[{'reason':'outside_complete_class_I_alpha_beta_scope','chain_kinds':kinds}]
    cache={}
    def score(a,b):
        key=tuple(sorted((a,b)))
        if key not in cache:cache[key]=contact_count(chains[a],chains[b])
        return cache[key]
    matrix=np.array([[score(a,b) for b in be] for a in al])
    ia,ib=linear_sum_assignment(-matrix);used=set()
    for i,j in zip(ia,ib):
        a,b=al[i],be[j];n=int(matrix[i,j])
        if n<3 or any(matrix[i,k]>=n for k in range(len(be)) if k!=j) or any(matrix[k,j]>=n for k in range(len(al)) if k!=i):
            excluded.append({'reason':'ambiguous_or_absent_alpha_beta_contacts','alpha':a,'beta':b});continue
        def choose(candidates,metric):
            ranked=sorted(((metric(c),c) for c in candidates),reverse=True)
            if not ranked or ranked[0][0]<=0 or (len(ranked)>1 and ranked[0][0]==ranked[1][0]):return None
            return ranked[0][1]
        p=choose(pe,lambda c:score(a,c)+score(b,c))
        h=choose(mh,lambda c:score(p,c)) if p else None
        z=choose(b2,lambda c:score(h,c)) if h else None
        if not p or not h or not z or p in used or not(score(a,h)+score(b,h)):
            excluded.append({'reason':'ambiguous_or_missing_pmhc_assignment','alpha':a,'beta':b});continue
        used.add(p);keep={a,b,p,h,z};rs=[r for r in residues if r.chain in keep]
        units.append((rs,{'alpha_chain':a,'beta_chain':b,'peptide_chain':p,'mhc_chain':h,'b2m_chain':z,'alpha':chain_sequence(chains[a]),'beta':chain_sequence(chains[b]),'peptide':chain_sequence(chains[p]),'identity_basis':'PDB_COMPND+unique_contact_assignment','mhc_class':'I'}))
    return units,excluded


# Canonical heavy atoms required for peptide side-chain graft controls.
_SIDECHAINS = dict(A='CB',R='CB CG CD NE CZ NH1 NH2',N='CB CG OD1 ND2',D='CB CG OD1 OD2',C='CB SG',Q='CB CG CD OE1 NE2',E='CB CG CD OE1 OE2',G='',H='CB CG ND1 CD2 CE1 NE2',I='CB CG1 CG2 CD1',L='CB CG CD1 CD2',K='CB CG CD CE NZ',M='CB CG SD CE',F='CB CG CD1 CD2 CE1 CE2 CZ',P='CB CG CD',S='CB OG',T='CB OG1 CG2',W='CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',Y='CB CG CD1 CD2 CE1 CE2 CZ OH',V='CB CG1 CG2')

def validate_peptide_atoms(residues):
    from .graph import AA3TO1
    peptide=[r for r in residues if r.role=='peptide']
    if not peptide:raise ValueError('missing explicit peptide')
    for r in peptide:
        aa=AA3TO1.get(r.resname.upper())
        if aa not in _SIDECHAINS:raise ValueError('noncanonical peptide residue')
        required={'N','CA','C','O'} | set(_SIDECHAINS[aa].split())
        missing=required-{a[0] for a in r.atoms}
        if missing:raise ValueError(f'incomplete peptide heavy atoms {r.chain}:{r.resid}{r.icode}: {sorted(missing)}')
