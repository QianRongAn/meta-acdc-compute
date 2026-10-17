"""Build class-I alpha/beta structural-decoy data, with fixed family folds.
Binding/activation labels are NOT supplied by this structural task.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,random,re,subprocess
from difflib import SequenceMatcher
from pathlib import Path
import numpy as np
from .release import STRUCTURAL_TASK,sha256
from ..structure.artifact_qc import require_construction_screen
from ..structure.feature_schema import FEATURE_SCHEMA
from ..structure.graph import parse_pdb,build_interface_graph_from_residues,AA_ORDER
from ..structure.identity import class_i_units,chain_sequence,validate_peptide_atoms
from ..models.dataset import make_decoy,make_displaced_peptide_decoy,_peptide_chains


def domain(sequence):
    m=re.search(r'C([A-Z]{4,35})F[GAS][A-Z]G',sequence)
    return (sequence[:m.end()+4],m.group(1)) if m else (sequence,'')


def similar(a,b):
    av,ac=domain(a);bv,bc=domain(b)
    if min(len(av),len(bv))/max(len(av),len(bv))<.9:return False
    if SequenceMatcher(None,av,bv,autojunk=False).ratio()<.9:return False
    return not(ac and bc) or SequenceMatcher(None,ac,bc,autojunk=False).ratio()>=.9


def assign_folds(records,seed=42):
    n=len(records);parent=list(range(n))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    def union(i,j):
        a,b=find(i),find(j)
        if a!=b:parent[max(a,b)]=min(a,b)
    for i,a in enumerate(records):
        for j in range(i):
            b=records[j]
            if a['pdb']==b['pdb'] or (similar(a['alpha'],b['alpha']) and similar(a['beta'],b['beta'])):union(i,j)
    grouped=collections.defaultdict(list)
    for i in range(n):grouped[find(i)].append(i)
    families=[]
    for members in grouped.values():
        names=sorted(records[i]['complex_id'] for i in members)
        fid='F'+hashlib.sha256('|'.join(names).encode()).hexdigest()[:16];families.append(fid)
        for i in members:records[i]['family_id']=fid
    families.sort();random.Random(seed).shuffle(families)
    if len(families)<10:raise ValueError('too few isolated receptor families for four independent folds')
    nt=max(1,round(.15*len(families)));nc=max(1,round(.10*len(families)));nv=max(1,round(.15*len(families)))
    folds={f:'test' for f in families[:nt]};folds.update({f:'calibration' for f in families[nt:nt+nc]});folds.update({f:'validation' for f in families[nt+nc:nt+nc+nv]})
    for x in records:x['fold']=folds.get(x['family_id'],'train')
    return len(families)


def inherit_folds(records, previous):
    old={r['complex_id']:r for r in previous if r['domain']=='crystal'}
    retained=[];excluded=[]
    for row in records:
        before=old.get(row['complex_id'])
        if before is None:
            excluded.append(dict(complex_id=row['complex_id'],reason='not_in_frozen_reference_release'));continue
        if any(row[k]!=before[k] for k in ('alpha','beta','peptide','pdb')):
            raise ValueError('frozen complex identity changed: '+row['complex_id'])
        row['family_id']=before['family_id'];row['fold']=before['fold'];retained.append(row['complex_id'])
    return set(retained),excluded


def build(*args,**kwargs):
    raise ValueError('legacy per-residue graft and fixed-axis builders are retired; rebuild from --previous-release-dir')


def main():
    parser=argparse.ArgumentParser(description='Rebuild paired rigid-peptide structural controls with inherited families')
    parser.add_argument('--previous-release-dir',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from .repair_release import rebuild
    rebuild(args.previous_release_dir,args.out)


if __name__=='__main__':main()
