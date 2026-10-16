"""Build class-I alpha/beta structural-decoy data, with fixed family folds.
Binding/activation labels are NOT supplied by this structural task.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,random,re,subprocess
from difflib import SequenceMatcher
from pathlib import Path
import numpy as np
from .release import STRUCTURAL_TASK,sha256
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


def build(source,output,with_electrostatics=False,predicted_root=None):
    import torch
    if (output/'graphs.pt').exists():raise ValueError('immutable release graph already exists; choose a new version')
    output.mkdir(parents=True,exist_ok=True);source=source.resolve();units=[];excluded=[];raw_hashes={}
    with (source/'complexes.tsv').open() as f:pdbs=sorted({x['pdb'].lower() for x in csv.DictReader(f,delimiter='\t')})
    for pdb in pdbs:
        path=source/f'{pdb}.pdb';raw_hashes[str(path)]=sha256(path)
        try:
            residues=parse_pdb(path);parts,bad=class_i_units(residues)
            excluded.extend(dict(pdb=pdb,**x) for x in bad)
            for i,(rs,meta) in enumerate(parts):
                identifier=pdb if len(parts)==1 else f'{pdb}_{i}'
                meta.update(domain='crystal',mhc_sequence=chain_sequence([x for x in rs if x.chain==meta['mhc_chain']]),b2m_sequence=chain_sequence([x for x in rs if x.chain==meta['b2m_chain']]),complex_id=identifier,pdb=pdb,source_path=str(path),source_sha256=raw_hashes[str(path)])
                if 'X' in meta['alpha']+meta['beta']+meta['peptide'] or any(a not in AA_ORDER for a in meta['peptide']):raise ValueError('noncanonical sequence; manual review')
                try:
                    validate_peptide_atoms(rs)
                    graph=build_interface_graph_from_residues(rs,electrostatics=with_electrostatics)
                    if graph.n_nodes<20:raise ValueError('too few interface nodes')
                    units.append((rs,meta,graph))
                except ValueError as e:excluded.append(dict(pdb=pdb,complex_id=identifier,reason=str(e)))
        except Exception as e:excluded.append(dict(pdb=pdb,reason=f'{type(e).__name__}: {e}'))
    records=[x[1] for x in units];nf=assign_folds(records)
    if predicted_root is not None:
        from ..structure.quality import build_qc,annotate_prediction
        from ..structure.cif import parse_cif
        qc=build_qc([predicted_root],records,output/'native_qc.tsv')
        by_id={x['complex_id']:x for x in records}
        predicted=[]
        for q in qc:
            if q['flag'] or q['native']!='yes' or not q['model_id'].endswith('_model_0'):
                continue
            ref=by_id[q['reference_id']];path=Path(q['source_path']);rs=annotate_prediction(parse_cif(path),ref,ref['peptide'])
            meta=dict(ref,domain='predicted',complex_id='af3-'+q['model_id'],source_path=str(path),source_sha256=q['file_sha256'],qc_model_id=q['model_id'],confidence_path=q['confidence_path'],confidence_sha256=q['confidence_sha256'])
            for field,kind in [('alpha_chain','tra'),('beta_chain','trb'),('mhc_chain','mhc_i'),('b2m_chain','b2m'),('peptide_chain','peptide')]:
                meta[field]=next(x.chain for x in rs if x.chain_kind==kind)
            meta['reference_id']=ref['complex_id']
            meta['reference_alpha']=ref['alpha'];meta['reference_beta']=ref['beta']
            meta['alpha']=chain_sequence([x for x in rs if x.chain==meta['alpha_chain']])
            meta['beta']=chain_sequence([x for x in rs if x.chain==meta['beta_chain']])
            meta['identity_basis']='reference_observed_variable_sequence_match+per_model_QC'
            g=build_interface_graph_from_residues(rs,electrostatics=with_electrostatics)
            predicted.append((rs,meta,g));raw_hashes[str(path)]=q['file_sha256']
        units.extend(predicted);records=[x[1] for x in units]
        # Predicted copies inherit the already frozen reference family/fold.
        print(f'accepted predicted natives={len(predicted)}',flush=True)
    graphs=[];labels=[];names=[];masks=[];nodes=[];provenance=[];seen={};rng=random.Random(43)
    def add(graph,meta,kind,donor=None):
        arrays=(graph.node_features,graph.node_coords,graph.edge_index,graph.edge_features);h=hashlib.sha256()
        for a in arrays:
            x=np.asarray(a);h.update(str(x.shape).encode());h.update(x.tobytes())
        digest=h.hexdigest();label=int(kind=='native')
        if digest in seen:
            excluded.append({'complex_id':meta['complex_id'],'kind':kind,'donor':donor,'reason':'identical_model_input','existing':seen[digest]});return False
        if not graph.peptide_mask or not any(graph.peptide_mask):raise ValueError('no peptide nodes')
        seen[digest]=meta['complex_id']+':'+kind
        graphs.append(arrays);labels.append(label);names.append(meta['complex_id'] if label else meta['complex_id']+('-disp' if kind=='displaced' else '<-'+donor))
        masks.append(graph.peptide_mask);nodes.append(graph.node_keys)
        provenance.append(dict(meta,graph_id=names[-1],kind=kind,label=label,label_scope=STRUCTURAL_TASK,donor=donor,input_sha256=digest))
        return True
    native=[]
    for rs,meta,g in units:
        if add(g,meta,'native'):native.append((rs,meta,g))
    for k,(rs,meta,g) in enumerate(native):
        forced=set(g.node_keys)
        try:add(build_interface_graph_from_residues(make_displaced_peptide_decoy(rs),forced_nodes=forced,electrostatics=with_electrostatics),meta,'displaced')
        except ValueError as e:excluded.append({'complex_id':meta['complex_id'],'kind':'displaced','reason':str(e)})
        pool=[x for x in native if x[1]['fold']==meta['fold'] and len(x[1]['peptide'])==len(meta['peptide']) and x[1]['peptide']!=meta['peptide']]
        rng.shuffle(pool)
        for dr,dm,dg in pool[:20]:
            try:
                foreign=next(iter(_peptide_chains(dr).values()));mutated=make_decoy(rs,foreign)
                gg=build_interface_graph_from_residues(mutated,forced_nodes=forced,electrostatics=with_electrostatics)
                if add(gg,meta,'graft',dm['complex_id']):break
            except ValueError:continue
        else:excluded.append({'complex_id':meta['complex_id'],'kind':'graft','reason':'no_valid_distinct_same_length_same_fold_donor'})
        if k%40==0:print(f'build {k}/{len(native)} natives; graphs={len(graphs)}',flush=True)
    splits={fold:[i for i,x in enumerate(provenance) if x['fold']==fold] for fold in ['train','validation','calibration','test']}
    for fold,idx in splits.items():
        if not idx or {labels[i] for i in idx}!={0,1}:raise ValueError('unusable frozen fold '+fold)
    data={'graphs':graphs,'labels':labels,'pdbs':names,'peptide_masks':masks,'node_keys':nodes,'records':provenance,'provenance':{'schema_version':1,'task':STRUCTURAL_TASK,'identity':'depositor_metadata+unique_contacts','family_rule':'same_PDB_or_alpha_and_beta_variable_domain_90pct_and_CDR3_90pct_when_available','seed':42,'raw_sha256':raw_hashes,'feature_dims':[len(graphs[0][0][0]),len(graphs[0][3][0])]}}
    torch.save(data,output/'graphs.pt')
    (output/'structure_records.json').write_text(json.dumps(records,indent=2))
    (output/'graph_records.json').write_text(json.dumps(provenance,indent=2))
    (output/'exclusions.json').write_text(json.dumps(excluded,indent=2))
    (output/'splits.json').write_text(json.dumps(splits,indent=2))
    with (output/'family_map.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['complex_id','pdb','family_id','fold','alpha','beta','peptide'],delimiter='\t');w.writeheader()
        for meta in records:w.writerow({k:meta[k] for k in w.fieldnames})
    stats={'n_raw_PDB':len(pdbs),'n_native':sum(labels),'n_graphs':len(graphs),'n_families':nf,'n_exclusions':len(excluded),'fold_graphs':{k:len(v) for k,v in splits.items()},'task':STRUCTURAL_TASK,'electrostatics':with_electrostatics}
    (output/'build_stats.json').write_text(json.dumps(stats,indent=2));print(json.dumps(stats),flush=True)
    return stats


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--structures',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--with-electrostatics',action='store_true');p.add_argument('--predicted-root',type=Path);a=p.parse_args();build(a.structures,a.out,a.with_electrostatics,a.predicted_root)
if __name__=='__main__':main()
