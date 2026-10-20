"""Independent array, identity, source and fold validation for a draft release."""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from .release import STRUCTURAL_TASK,sha256,frozen_indices


def validate(folder):
    import torch
    data=torch.load(folder/'graphs.pt',map_location='cpu',weights_only=False)
    errors=[];seen={};folds=defaultdict(set);pdb_folds=defaultdict(set);counts=Counter();source_cache={}
    if len({len(data[k]) for k in ['graphs','labels','records','peptide_masks','node_keys']})!=1:raise ValueError('unaligned dataset arrays')
    from ..structure.feature_schema import FEATURE_SCHEMA
    policy=data.get('provenance',{})
    if policy.get('feature_schema')!=FEATURE_SCHEMA or not policy.get('complete_interface_heavy_atoms') or policy.get('construction_screen')!='construction-screen-v1':
        errors.append('new atomic completeness/construction QC policy requires data regeneration')
    from ..structure.artifact_qc import load_record_residues, construction_screen
    from ..structure.identity import missing_heavy_atoms
    from ..models.dataset import make_decoy, make_displaced_peptide_decoy
    from ..structure.rigid_controls import atomic_screen,transfer_peptide,policy as atomic_policy,policy_sha256
    from ..structure.graph import build_interface_graph_from_residues
    from .repair_release import without_confidence
    atomic_release=policy.get('atomic_policy')==atomic_policy()['policy_id']
    if not atomic_release or policy.get('atomic_policy_sha256')!=policy_sha256():
        errors.append('rigid atomic repair policy missing or changed; legacy graft/translation release refused')
    if not policy.get('confidence_masked'):errors.append('confidence masking must be identical for both classes')
    reconstructed={}
    def original(identifier):
        if identifier not in reconstructed:
            reconstructed[identifier]=load_record_residues(meta[identifier])
        return reconstructed[identifier]
    inherited_path=policy.get('inherited_records_path')
    inherited_digest=policy.get('inherited_records_sha256')
    if not inherited_path or not inherited_digest or sha256(inherited_path)!=inherited_digest:
        errors.append('frozen reference inheritance missing or changed')
    else:
        previous={r['complex_id']:r for r in json.loads(Path(inherited_path).read_text())}
        for row in json.loads((folder/'structure_records.json').read_text()):
            before=previous.get(row['complex_id'])
            if before is None or any(row[k]!=before[k] for k in ('family_id','fold','alpha','beta','peptide','pdb','domain','source_path','source_sha256')):
                errors.append('frozen receptor identity/fold changed: '+row['complex_id'])
    records=data['records'];meta={x['complex_id']:x for x in json.loads((folder/'structure_records.json').read_text())}
    splits=json.loads((folder/'splits.json').read_text())
    try:frozen_indices({'splits':splits},data)
    except ValueError as exc:errors.append('training folds unusable: '+str(exc))
    qc={x['model_id']:x for x in csv.DictReader((folder/'native_qc.tsv').open(),delimiter='\t')} if (folder/'native_qc.tsv').exists() else {}
    for i,(graph,y,row,mask,keys) in enumerate(zip(data['graphs'],data['labels'],records,data['peptide_masks'],data['node_keys'])):
        nf,x,ei,ef=map(np.asarray,graph);name=row['graph_id'];counts[row['kind']]+=1
        def fail(msg):errors.append(name+': '+msg)
        if len(nf)!=len(x) or len(keys)!=len(nf) or len(mask)!=len(nf):fail('unaligned arrays');continue
        if not all(np.isfinite(a).all() for a in (nf,x,ei,ef)):fail('nonfinite')
        if len(nf)<20 or not np.any(mask) or not np.array_equal(np.asarray(mask),nf[:,24]==1):fail('peptide mask invalid')
        if not np.allclose(nf[:,4:24].sum(1),1):fail('unknown amino-acid identity')
        if ei.ndim!=2 or ei.shape[1]!=2 or len(ei)!=len(ef) or not len(ei) or ei.min()<0 or ei.max()>=len(nf):fail('edge indices invalid');continue
        distances=np.linalg.norm(x[:,None]-x[None,:],axis=-1)
        expected={(a,b) for a,b in zip(*np.where((distances<=8)&(~np.eye(len(x),dtype=bool))))}
        if set(map(tuple,ei))!=expected:fail('edge topology disagrees with coordinates')
        ed=distances[ei[:,0],ei[:,1]]
        if not np.allclose(ef[:,:12],np.exp(-((ed[:,None]-np.arange(12)*8/12)**2)/2),atol=1e-7):fail('RBF distance mismatch')
        digest=hashlib.sha256()
        for a in (nf,x,ei,ef):digest.update(str(a.shape).encode());digest.update(a.tobytes())
        h=digest.hexdigest()
        if h!=row['input_sha256']:fail('stored input hash mismatch')
        if h in seen:fail('duplicate graph input '+seen[h])
        seen[h]=name
        try:
            atoms=original(row['complex_id'])
            if row['kind']=='rigid_transfer':
                donor=meta[row['donor']]
                for field in ('fold','domain','mhc_sequence','b2m_sequence'):
                    if donor[field]!=row[field]:fail('invalid rigid donor matching: '+field)
                if donor['peptide']==row['peptide'] or len(donor['peptide'])!=len(row['peptide']):fail('invalid donor peptide')
                atoms,transform=transfer_peptide(atoms,original(row['donor']))
                if transform!=row.get('rigid_transfer'):fail('rigid geometry evidence mismatch')
            elif row['kind']=='graft':atoms=make_decoy(atoms,[r for r in original(row['donor']) if r.role=='peptide'])
            elif row['kind']=='displaced':atoms=make_displaced_peptide_decoy(atoms)
            expected_keys=set(map(tuple,keys))
            selected=[r for r in atoms if (r.chain,r.resid,r.icode) in expected_keys]
            if len(selected)!=len(keys) or any(missing_heavy_atoms(r) for r in selected):fail('source has incomplete or missing interface residues')
            if construction_screen(atoms)!=row.get('construction_qc'):fail('construction QC does not reproduce from source atoms')
            if atomic_release:
                if row['kind'] not in {'native','rigid_transfer'}:fail('retired control type')
                result=atomic_screen(atoms)
                if not result['passed'] or result!=row.get('atomic_qc'):fail('atomic QC does not reproduce')
                rebuilt=build_interface_graph_from_residues(without_confidence(atoms),forced_nodes=expected_keys)
                expected_arrays=(rebuilt.node_features,rebuilt.node_coords,rebuilt.edge_index,rebuilt.edge_features)
                for observed,expected_array in zip(graph,expected_arrays):
                    if np.asarray(observed).shape!=np.asarray(expected_array).shape or not np.allclose(observed,expected_array,rtol=0,atol=1e-10):
                        fail('model inputs do not reproduce from original atoms')
                if not np.all(nf[:,3]==0.5):fail('native/control confidence shortcut')
        except (ValueError,KeyError) as exc:fail('source atom reconstruction: '+str(exc))
        if not row.get('construction_qc',{}).get('passed') or row.get('construction_qc',{}).get('policy')!='construction-screen-v1':fail('missing construction QC')
        if y!=int(row['kind']=='native') or row['label_scope']!=STRUCTURAL_TASK:fail('label semantics')
        if row['donor'] and meta[row['donor']]['fold']!=row['fold']:fail('cross-fold graft donor')
        folds[row['family_id']].add(row['fold']);pdb_folds[row['pdb']].add(row['fold'])
        source=Path(row['source_path'])
        if str(source) not in source_cache:source_cache[str(source)]=sha256(source)
        if source_cache[str(source)]!=row['source_sha256']:fail('changed source structure')
        if row['domain']=='predicted':
            q=qc.get(row['qc_model_id'])
            if not q or q['flag'] or q['native']!='yes' or q['file_sha256']!=row['source_sha256']:fail('prediction was not accepted by model-level QC')
            if sha256(row['confidence_path'])!=row['confidence_sha256']:fail('changed confidence file')
    if atomic_release:
        paired=defaultdict(list)
        for i,row in enumerate(records):paired[row.get('pair_id')].append(i)
        for pair,indices in paired.items():
            if pair is None or len(indices)!=2 or {records[i]['kind'] for i in indices}!={'native','rigid_transfer'}:
                errors.append('unmatched native/control pair: '+str(pair));continue
            a,b=indices
            if data['node_keys'][a]!=data['node_keys'][b] or data['peptide_masks'][a]!=data['peptide_masks'][b]:
                errors.append('pair node-set/mask shortcut: '+str(pair))
            control=next(records[i] for i in indices if records[i]['kind']=='rigid_transfer')
            changed,_=transfer_peptide(original(pair),original(control['donor']))
            native_keys=set(build_interface_graph_from_residues(without_confidence(original(pair))).node_keys)
            control_keys=set(build_interface_graph_from_residues(without_confidence(changed)).node_keys)
            if set(map(tuple,data['node_keys'][a]))!=native_keys|control_keys:
                errors.append('pair does not use independently reconstructed union interface: '+str(pair))
    if len(records)!=len(data['graphs']) or len(seen)!=len(data['graphs']):errors.append('data/record cardinality mismatch')
    if any(len(v)!=1 for v in folds.values()) or any(len(v)!=1 for v in pdb_folds.values()):errors.append('family/PDB leakage across frozen folds')
    # Independent exact variable-chain overlap check across all four folds.
    from .build_release import domain,similar
    items=list(meta.values())
    from .experimental_evidence import assert_reserved_absent
    try:assert_reserved_absent(items)
    except ValueError as exc:errors.append(str(exc))
    for i,a in enumerate(items):
        for b in items[:i]:
            if a['fold']!=b['fold'] and similar(a['alpha'],b['alpha']) and similar(a['beta'],b['beta']):
                errors.append('near-receptor cross-fold overlap: '+a['complex_id']+'/'+b['complex_id'])
    with (folder/'vdjdb.associations.tsv').open() as f:
        vdj=list(csv.DictReader(f,delimiter='\t'))
    if not vdj or any(x['label'] or x['label_scope']!='unlabeled' for x in vdj):errors.append('VDJdb carries fabricated binding labels')
    try:assert_reserved_absent(vdj)
    except ValueError as exc:errors.append(str(exc))
    assay=Counter()
    with (folder/'iedb.assays.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            assay[row['endpoint']]+=1
            if row['label_scope'] not in {'mhc_assay','mhc_presentation'} or not row['assay_id']:errors.append('IEDB scope/provenance missing');break
    if not assay:errors.append('empty IEDB evidence ledger')
    with (folder/'clinical_index_cases.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if ('A118T' in row['tcr'] or 'unverified' in row['evidence'].lower()) and row['pdb']:errors.append('unverified clinical receptor linked to structure')
    engineering_errors=[e for e in errors if not e.startswith('training folds unusable: ')]
    return {'atomic_engineering_passed':not engineering_errors,'training_ready':not errors,'passed':not errors,'errors':errors,'n_graphs':len(records),'graph_kinds':dict(counts),'n_native':counts['native'],'n_families':len(folds),'n_vdjdb':len(vdj),'n_iedb_assays':sum(assay.values()),'iedb_endpoints':dict(assay),'domains':sorted({x['domain'] for x in records}),'fold_counts':{k:len(v) for k,v in splits.items()},'folds':splits,'task':STRUCTURAL_TASK,'binding_training_allowed':False}


def publish(folder,code_root):
    import subprocess
    report=validate(folder)
    (folder/'validation.json').write_text(json.dumps(report,indent=2))
    if not report['passed']:raise ValueError('data validation failed: '+str(report['errors'][:8]))
    manifest=folder/'release_manifest.json'
    if manifest.exists():raise ValueError('release is immutable; refusing to replace a manifest')
    code_files={}
    for p in (code_root/'src/meta_acdc').rglob('*.py'):code_files[str(p.relative_to(code_root))]=sha256(p)
    code_files['scripts/pretrain_audit.py']=sha256(code_root/'scripts/pretrain_audit.py')
    code_files['docs/structure_repair_20261017/policy_v2.json']=sha256(code_root/'docs/structure_repair_20261017/policy_v2.json')
    code_files['docs/evidence/holdout_registry.json']=sha256(code_root/'docs/evidence/holdout_registry.json')
    files={p.name:{'sha256':sha256(p),'bytes':p.stat().st_size} for p in sorted(folder.iterdir()) if p.is_file() and not p.name.endswith('.tmp')}
    sources={}
    d=json.loads((folder/'structure_records.json').read_text())
    for x in d:sources[x['source_path']]=x['source_sha256']
    import torch
    provenance=torch.load(folder/'graphs.pt',map_location='cpu',weights_only=False)['provenance']
    sources[provenance['inherited_records_path']]=provenance['inherited_records_sha256']
    m={'schema_version':1,'status':'validated','task':STRUCTURAL_TASK,'training_allowed':False,'acceptance_status':'machine_validation_only','acceptance_blocker':'A separate explicit training release decision is required; publishing data does not start or authorize training','binding_training_allowed':False,'atomic_policy':'rigid-observed-peptide-transfer-v2','mask_plddt_required':True,'clinical_training_allowed':False,'domains':report['domains'],'dataset':'graphs.pt','code_commit':subprocess.run(['git','-C',str(code_root),'rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip(),'code_files':code_files,'files':files,'source_structures':sources,'splits':report['folds'],'scope':'class-I alpha-beta crystal/predicted native versus synthetic decoys; not experimentally validated binding probabilities'}
    manifest.write_text(json.dumps(m,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='folds'}),flush=True)
    return m


def main():
    p=argparse.ArgumentParser();p.add_argument('--release-dir',type=Path,required=True);p.add_argument('--publish',action='store_true');a=p.parse_args()
    if a.publish:publish(a.release_dir,Path(__file__).resolve().parents[3])
    else:
        report=validate(a.release_dir);print(json.dumps({k:v for k,v in report.items() if k!='folds'}),flush=True)
        raise SystemExit(0 if report['passed'] else 2)
if __name__=='__main__':main()
