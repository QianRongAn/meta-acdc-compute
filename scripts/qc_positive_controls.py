"""Rebuild candidate native graphs with existing atomic policies; no training labels."""
from pathlib import Path
import json,hashlib,argparse
from dataclasses import replace
import numpy as np
import torch
from meta_acdc.structure.graph import parse_pdb,build_interface_graph_from_residues
from meta_acdc.structure.identity import annotate,class_i_units,chain_sequence
from meta_acdc.structure.receptor_identity_audit import seqres_sequences,is_subsequence
from meta_acdc.structure.rigid_controls import atomic_screen,policy_sha256
ROOT=Path(__file__).resolve().parents[1]
def local_path(value):
    value=value.replace('\\','/')
    return Path('/mnt/'+value[0].lower()+value[2:]) if len(value)>2 and value[1]==':' else Path(value)
def main():
    base=ROOT/'docs/fingerprint_positive_controls_20261027';ap=argparse.ArgumentParser();ap.add_argument('--out-dir',type=Path,default=base);args=ap.parse_args();out=args.out_dir;out.mkdir(parents=True,exist_ok=True);data=json.loads((base/'source_balanced_design.json').read_text());manifest=json.loads((base/'coordinate_download_manifest.json').read_text());by_pdb={}
    for x in data['rows']:
        if x['status']=='identity_admitted_structure_QC_pending':by_pdb.setdefault(x['pdb'],[]).append(x)
    rows=[];graphs=[];records=[]
    for source in manifest:
        path=(ROOT/source['repository_snapshot']) if source.get('repository_snapshot') else local_path(source['local_cache']);assert hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256'];pdb=source['pdb'];
        if pdb not in by_pdb:
            rejected=[x for x in data['rows'] if x['pdb']==pdb]
            rows.append(dict(pdb=pdb,status='quarantined_before_atomic_QC',stage='identity_role_gate',reasons=sorted({reason for x in rejected for reason in x.get('reasons',[])}),source_sha256=source['sha256']));continue
        residues=parse_pdb(path);annotate(residues,path);seqres=seqres_sequences(path.read_text())
        for identity in by_pdb.get(pdb,[]):
            for role,kind in [('alpha','tra'),('beta','trb')]:
                for residue in residues:
                    if residue.chain==identity[role+'_chain']:residue.role='tcr';residue.chain_kind=kind
        units,excluded=class_i_units(residues)
        rows.append(dict(pdb=pdb,stage='unit_extraction',n_units=len(units),excluded=excluded,source_sha256=source['sha256']))
        # Compact-MHC constructs lacking beta2m remain outside the existing graph scope.
        for unit,meta in units:
            row=dict(pdb=pdb,alpha_chain=meta['alpha_chain'],beta_chain=meta['beta_chain'],status='quarantined',reasons=[],source_sha256=source['sha256'])
            try:
                candidates=[x for x in by_pdb.get(pdb,[]) if x['alpha_chain']==meta['alpha_chain'] and x['beta_chain']==meta['beta_chain']]
                if len(candidates)!=1:raise ValueError('no unique complete-variable identity for contact-selected unit')
                identity=candidates[0]
                for role in ['alpha','beta']:
                    if seqres.get(meta[role+'_chain'])!=identity[role]['full_construct']:raise ValueError('coordinate-file SEQRES differs from primary identity header')
                    if not is_subsequence(meta[role],identity[role]['full_construct']):raise ValueError('observed atoms disagree with complete construct')
                graph=build_interface_graph_from_residues([replace(r,plddt=None) for r in unit]);screen=atomic_screen(unit);row['atomic_qc']=screen
                if not screen['passed']:raise ValueError('existing fixed native atomic screen rejected unit')
                row['unit_metadata']=dict(meta,pmhc_construct_key=hashlib.sha256(('|'.join(seqres[meta[k]] for k in ['mhc_chain','b2m_chain','peptide_chain'])).encode()).hexdigest())
                row.update(status='native_graph_QC_pass',receptor_id=identity['receptor_id'],V_context=identity['V_context'],source_study_key=identity['source_study_key'],n_nodes=graph.n_nodes)
                arrays=(graph.node_features,graph.node_coords,graph.edge_index,graph.edge_features);h=hashlib.sha256()
                for a in arrays:h.update(str(np.asarray(a).shape).encode());h.update(np.asarray(a).tobytes())
                row['graph_sha256']=h.hexdigest();records.append(dict(row,peptide_mask=graph.peptide_mask,node_keys=graph.node_keys));graphs.append(arrays)
            except ValueError as error:row['reasons'].append(str(error))
            rows.append(row)
    torch.save(dict(task='cross-source native candidates only',trainable=False,binding_labels_available=False,graphs=graphs,records=records,atomic_policy_sha256=policy_sha256()),out/'native_graph_candidates.pt')
    (out/'native_atomic_qc.json').write_text(json.dumps(dict(rows=rows,n_graphs=len(graphs),n_pdbs_with_pass=len({x['pdb'] for x in records}),scope='existing complete class-I alpha/beta graph scope and fixed atomic policy; no imputations, no relaxed thresholds',model_training_started=False),indent=2)+'\n')
    print('QC_PASS_GRAPHS',len(graphs),'PDBS',len({x['pdb'] for x in records}))
if __name__=='__main__':main()
