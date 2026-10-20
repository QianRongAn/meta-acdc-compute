"""Rebuild atom-screened, paired structural controls from immutable source records.

No new train/test assignments, atom imputation, old local grafts or fixed-axis
translations. Native inventory and paired training arrays have separate counts.
"""
from __future__ import annotations
import collections, copy, csv, hashlib, json, os
from pathlib import Path
import numpy as np
from dataclasses import replace
from .release import sha256,STRUCTURAL_TASK
from ..structure.artifact_qc import load_record_residues
from ..structure.rigid_controls import atomic_screen,transfer_peptide,policy,policy_sha256,POLICY_PATH
from ..structure.identity import chain_sequence
from ..structure.graph import build_interface_graph_from_residues
from ..structure.feature_schema import FEATURE_SCHEMA


def without_confidence(residues):
    return [replace(r,plddt=None) for r in residues]


def graph_hash(arrays):
    digest=hashlib.sha256()
    for a in arrays:
        x=np.asarray(a);digest.update(str(x.shape).encode());digest.update(x.tobytes())
    return digest.hexdigest()


def rebuild(previous,output):
    import torch
    previous=previous.resolve();output=output.resolve()
    initial_policy_sha=policy_sha256()
    if output.exists():raise ValueError('new release directory must not exist')
    previous_manifest=json.loads((previous/'release_manifest.json').read_text())
    # Forensic source version may be revoked, but its accepted bytes must still match.
    for name,record in previous_manifest['files'].items():
        if sha256(previous/name)!=record['sha256']:raise ValueError('changed historical artifact: '+name)
    old=json.loads((previous/'structure_records.json').read_text())
    exclusions=[];inventory=[];loaded={};native_graphs={};qc_rows=[];corrections=[]
    from ..structure.quality import qc_file
    by_pdb=collections.defaultdict(list)
    for item in old:
        if item['domain']=='crystal':by_pdb[item['pdb']].append(item)
    # Crystal references must pass first; a predicted copy cannot revive a quarantined reference.
    for original in sorted(old,key=lambda r:(r['domain']!='crystal',r['complex_id'])):
        meta=copy.deepcopy(original);identifier=meta['complex_id']
        try:
            if sha256(meta['source_path'])!=meta['source_sha256']:raise ValueError('source hash changed')
            if meta['domain']=='predicted':
                clean_refs=collections.defaultdict(list)
                for r in inventory:
                    if r['domain']=='crystal':clean_refs[r['pdb']].append(r)
                q=qc_file(Path(meta['source_path']),clean_refs);qc_rows.append(q)
                if q['flag'] or q['native']!='yes':raise ValueError('predicted identity/confidence/QC: '+q['flag'])
                if q['reference_id']!=meta['reference_id']:raise ValueError('reference reassignment forbidden')
                if sha256(meta['confidence_path'])!=meta['confidence_sha256']:raise ValueError('confidence hash changed')
                rs=load_record_residues(meta,strict_identity=False)
                for key,field in [('alpha','alpha_chain'),('beta','beta_chain'),('peptide','peptide_chain')]:
                    actual=chain_sequence([r for r in rs if r.chain==meta[field]])
                    if actual!=meta[key]:raise ValueError('observed predicted receptor/peptide identity changed')
                for key,field in [('mhc_sequence','mhc_chain'),('b2m_sequence','b2m_chain')]:
                    actual=chain_sequence([r for r in rs if r.chain==meta[field]])
                    if actual!=meta[key]:
                        corrections.append(dict(complex_id=identifier,field=key,old=meta[key],observed=actual,
                                                basis='same hashed source atoms; reference identity revalidated'))
                        meta['historical_declared_'+key]=meta[key];meta[key]=actual
            else:rs=load_record_residues(meta)
            g=build_interface_graph_from_residues(without_confidence(rs))
            if g.n_nodes<20:raise ValueError('small interface')
            atom_qc=atomic_screen(rs)
            if not atom_qc['passed']:
                exclusions.append(dict(complex_id=identifier,domain=meta['domain'],fold=meta['fold'],
                                       stage='native_atoms',reason='preset_heavy_atom_screen',qc=atom_qc));continue
            meta['construction_qc']=atom_qc['gross_screen'];meta['atomic_qc']=atom_qc
            inventory.append(meta);loaded[identifier]=rs;native_graphs[identifier]=g
        except ValueError as exc:
            exclusions.append(dict(complex_id=identifier,domain=meta['domain'],fold=meta['fold'],
                                   stage='native_identity_or_completeness',reason=str(exc)))
        if len(inventory)%40==0:print('source records checked; accepted inventory',len(inventory),flush=True)
    from .experimental_evidence import assert_reserved_absent
    assert_reserved_absent(inventory)
    output.mkdir(parents=True)
    graphs=[];labels=[];names=[];masks=[];nodes=[];records=[];seen=set();pairs=[]
    def arrays(g):return (g.node_features,g.node_coords,g.edge_index,g.edge_features)
    def append(g,meta,kind,donor=None,transfer=None,atom_qc=None):
        name=meta['complex_id'] if kind=='native' else meta['complex_id']+'<-rigid-'+donor
        digest=graph_hash(arrays(g))
        if digest in seen:raise ValueError('duplicate model input')
        seen.add(digest);graphs.append(arrays(g));labels.append(int(kind=='native'));names.append(name)
        masks.append(g.peptide_mask);nodes.append(g.node_keys)
        item=dict(meta,graph_id=name,kind=kind,label=int(kind=='native'),label_scope=STRUCTURAL_TASK,
                  donor=donor,input_sha256=digest,pair_id=meta['complex_id'],policy_id=policy()['policy_id'])
        if atom_qc is not None:item.update(atomic_qc=atom_qc,construction_qc=atom_qc['gross_screen'])
        if transfer is not None:item['rigid_transfer']=transfer
        records.append(item)
    for receiver in sorted(inventory,key=lambda r:r['complex_id']):
        identifier=receiver['complex_id'];rs=loaded[identifier];ng=native_graphs[identifier]
        donors=[d for d in inventory if d['fold']==receiver['fold'] and d['domain']==receiver['domain']
                and d['mhc_sequence']==receiver['mhc_sequence'] and d['b2m_sequence']==receiver['b2m_sequence']
                and len(d['peptide'])==len(receiver['peptide']) and d['peptide']!=receiver['peptide']]
        accepted=None;attempts=[]
        for donor in sorted(donors,key=lambda r:r['complex_id']):
            try:
                changed,transfer=transfer_peptide(rs,loaded[donor['complex_id']])
                cq=atomic_screen(changed)
                if not cq['passed']:
                    attempts.append(dict(donor=donor['complex_id'],reason='preset_heavy_atom_screen',qc=cq));continue
                # Same UNION of current-atom interface nodes; no missingness or membership label shortcut.
                cg=build_interface_graph_from_residues(without_confidence(changed))
                keys=set(ng.node_keys)|set(cg.node_keys)
                same_native=build_interface_graph_from_residues(without_confidence(rs),forced_nodes=keys)
                same_control=build_interface_graph_from_residues(without_confidence(changed),forced_nodes=keys)
                if same_native.node_keys!=same_control.node_keys:raise ValueError('paired node identities differ')
                nh,ch=graph_hash(arrays(same_native)),graph_hash(arrays(same_control))
                if nh==ch or nh in seen or ch in seen:raise ValueError('identical model inputs')
                accepted=(same_native,same_control,donor,transfer,cq);break
            except ValueError as exc:attempts.append(dict(donor=donor['complex_id'],reason=str(exc)))
        if accepted is None:
            exclusions.append(dict(complex_id=identifier,domain=receiver['domain'],fold=receiver['fold'],
                stage='paired_control',reason='no_acceptable_same_fold_domain_MHC_B2M_donor',attempts=attempts));continue
        a,b,d,t,q=accepted
        append(a,receiver,'native');append(b,receiver,'rigid_transfer',d['complex_id'],t,q)
        pairs.append(dict(receiver=identifier,donor=d['complex_id'],fold=receiver['fold'],family_id=receiver['family_id'],
                          domain=receiver['domain'],n_nodes=a.n_nodes,native_atomic_qc=receiver['atomic_qc'],
                          control_atomic_qc=q,rigid_transfer=t,rejected_donor_attempts=attempts))
        if len(pairs)%20==0:print('retained pairs',len(pairs),flush=True)
    splits={fold:[i for i,r in enumerate(records) if r['fold']==fold] for fold in ['train','validation','calibration','test']}
    usable=all(idx and {labels[i] for i in idx}=={0,1} for idx in splits.values())
    if policy_sha256()!=initial_policy_sha:raise ValueError('atomic policy changed during reconstruction')
    provenance=dict(schema_version=1,feature_schema=FEATURE_SCHEMA,task=STRUCTURAL_TASK,
        complete_interface_heavy_atoms=True,construction_screen='construction-screen-v1',
        atomic_policy=policy()['policy_id'],atomic_policy_sha256=policy_sha256(),confidence_masked=True,
        inherited_records_path=str(previous/'structure_records.json'),inherited_records_sha256=sha256(previous/'structure_records.json'),
        inherited_manifest_path=str(previous/'release_manifest.json'),inherited_manifest_sha256=sha256(previous/'release_manifest.json'),
        feature_dims=[26,16],fixed_family_folds=True)
    torch.save(dict(graphs=graphs,labels=labels,pdbs=names,peptide_masks=masks,node_keys=nodes,records=records,provenance=provenance),output/'graphs.pt')
    # Unlabelled native assets are useful independently of the tiny synthetic-control pilot.
    native_assets=[]
    for meta in inventory:
        g=native_graphs[meta['complex_id']]
        native_assets.append(dict(record=meta,graph=arrays(g),peptide_mask=g.peptide_mask,node_keys=g.node_keys,
                                  input_sha256=graph_hash(arrays(g))))
    torch.save(dict(task='unlabeled_native_inventory',assets=native_assets,labels_are_available=False,
                    provenance=provenance),output/'native_inventory_graphs.pt')
    for name,data in [('structure_records.json',inventory),('graph_records.json',records),('exclusions.json',exclusions),
                      ('splits.json',splits),('paired_atomic_qc.json',pairs),('identity_corrections.json',corrections)]:
        (output/name).write_text(json.dumps(data,indent=2)+'\n')
    (output/'atomic_policy.json').write_bytes(POLICY_PATH.read_bytes())
    if qc_rows:
        with (output/'native_qc.tsv').open('w',newline='') as f:
            fields=list(dict.fromkeys(key for row in qc_rows for key in row))
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(qc_rows)
    for name in ['vdjdb.associations.tsv','vdjdb.reserved_quarantine.tsv','iedb.assays.tsv','iedb.assays.tsv.rejected.tsv',
                 'clinical_index_cases.tsv','reservation_report.json','input_sources.json']:
        if (previous/name).exists():os.link(previous/name,output/name)
    stats=dict(n_previous_native=len(old),n_clean_inventory=len(inventory),n_pairs=len(pairs),n_graphs=len(graphs),
        graph_kinds=dict(collections.Counter(r['kind'] for r in records)),
        inventory_domains=dict(collections.Counter(r['domain'] for r in inventory)),
        paired_domains=dict(collections.Counter(r['domain'] for r in pairs)),
        n_families=len({r['family_id'] for r in records}),fold_graph_counts={k:len(v) for k,v in splits.items()},
        usable_all_four_folds=usable,policy_id=policy()['policy_id'],thresholds_modified_after_build=False,
        training_allowed=False,binding_training_allowed=False,clinical_training_allowed=False,
        labels_are_experimental_nonbinding=False,n_identity_corrections=len(corrections),
        excluded_by_stage=dict(collections.Counter(r['stage'] for r in exclusions)))
    (output/'build_stats.json').write_text(json.dumps(stats,indent=2)+'\n')
    print(json.dumps(stats,indent=2),flush=True)
    return stats


def validate_native_inventory(folder):
    """Independently reproduce every unlabelled native graph from hashed source atoms."""
    import torch
    from ..structure.identity import missing_heavy_atoms
    data=torch.load(folder/'native_inventory_graphs.pt',map_location='cpu',weights_only=False)
    if not data.get('assets'):raise ValueError('empty native asset inventory')
    if data.get('task')!='unlabeled_native_inventory' or data.get('labels_are_available') is not False or 'labels' in data:
        raise ValueError('native asset inventory must not fabricate response labels')
    old={r['complex_id']:r for r in json.loads(Path(data['provenance']['inherited_records_path']).read_text())}
    seen=set();errors=[]
    for asset in data['assets']:
        row=asset['record'];identifier=row['complex_id']
        if identifier in seen:errors.append(identifier+': duplicate asset identity')
        seen.add(identifier)
        if sha256(row['source_path'])!=row['source_sha256']:errors.append(identifier+': source hash changed')
        if identifier not in old or any(row[k]!=old[identifier][k] for k in ['fold','family_id','alpha','beta','peptide','domain']):
            errors.append(identifier+': frozen identity or family changed')
        atoms=load_record_residues(row);screen=atomic_screen(atoms)
        if not screen['passed'] or screen!=row['atomic_qc']:errors.append(identifier+': atom QC does not reproduce')
        g=build_interface_graph_from_residues(without_confidence(atoms))
        expected=(g.node_features,g.node_coords,g.edge_index,g.edge_features)
        if asset['input_sha256']!=graph_hash(asset['graph']) or graph_hash(expected)!=asset['input_sha256']:
            errors.append(identifier+': graph arrays do not reproduce from source')
        if g.node_keys!=asset['node_keys'] or g.peptide_mask!=asset['peptide_mask']:errors.append(identifier+': node/mask differs')
        if g.missing_atom_nodes:errors.append(identifier+': retained interface has missing atoms')
        if any(r[3]!=0.5 for r in asset['graph'][0]):errors.append(identifier+': confidence is not masked')
    return dict(passed=not errors,errors=errors,n_assets=len(data['assets']),
                domains=dict(collections.Counter(a['record']['domain'] for a in data['assets'])),
                labels_are_available=False,missing_retained_interface_atoms=0 if not errors else None,
                training_allowed=False,binding_training_allowed=False,clinical_training_allowed=False)
