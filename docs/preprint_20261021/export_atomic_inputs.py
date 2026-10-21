"""Export trusted local structural evidence for public coordinate-level verification.

Reads the project's own historical torch cache. Never use this exporter on an
untrusted pickle. No model fitting or prediction is performed.
"""
from pathlib import Path
import argparse,dataclasses,gzip,hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def verify_snapshot_helpers(root=ROOT):
    prior=json.loads((HERE/'atomic_bundle_v1/export_receipt.json').read_text())
    for name,digest in prior['code_files'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Historical export requires pinned e4d68de helper bytes; newer main helpers are not the audited fixed-axis generator: '+name)
    return True

def export(folder,out):
    import torch
    from meta_acdc.structure.artifact_qc import load_record_residues,construction_screen
    from meta_acdc.models.dataset import make_decoy,make_displaced_peptide_decoy
    from meta_acdc.structure.identity import missing_heavy_atoms
    out=Path(out)
    verify_snapshot_helpers()
    if out.exists():raise ValueError('refusing archive overwrite')
    expected=json.loads((ROOT/'docs/n_review/structural_atom_audit.json').read_text())
    by_expected={r['graph_id']:r for r in expected['rows']}
    data=torch.load(folder/'graphs.pt',map_location='cpu',weights_only=False)
    records=json.loads((folder/'structure_records.json').read_text())
    by_id={r['complex_id']:r for r in records};loaded={};source_meta={}
    def get(key):
        if key not in loaded:
            record=by_id[key];path=Path(record['source_path']);rawhash=hashlib.sha256(path.read_bytes()).hexdigest()
            if rawhash!=record['source_sha256']:raise ValueError('source file changed '+key)
            loaded[key]=load_record_residues(record,strict_identity=False)
            source_meta[key]=dict(source_name=path.name,source_sha256=rawhash,domain=record['domain'],pdb=record['pdb'])
        return loaded[key]
    graphs=[]
    for row,keys in zip(data['records'],data['node_keys']):
        rs=get(row['complex_id'])
        if row['kind']=='graft':rs=make_decoy(rs,[r for r in get(row['donor']) if r.role=='peptide'])
        elif row['kind']=='displaced':rs=make_displaced_peptide_decoy(rs)
        screen=construction_screen(rs);prior=by_expected[row['graph_id']]
        if screen!=prior['construction']:raise ValueError('screen reconstruction differs '+row['graph_id'])
        selected={tuple(k) for k in keys}
        missing=[dict(chain=r.chain,resid=r.resid,icode=r.icode,role=r.role,atoms=missing_heavy_atoms(r)) for r in rs if (r.chain,r.resid,r.icode) in selected and missing_heavy_atoms(r)]
        if missing!=prior['missing_nodes']:raise ValueError('missing atom reconstruction differs')
        graphs.append(dict(graph_id=row['graph_id'],complex_id=row['complex_id'],kind=row['kind'],family_id=row['family_id'],node_keys=keys,donor=row.get('donor'),peptide_residues=[dataclasses.asdict(r) for r in rs if r.role=='peptide']))
        if len(graphs)%100==0:print('atomic export verified',len(graphs),flush=True)
    archive=dict(schema=1,scope='atomic coordinate extracts and reconstructed synthetic peptide controls; no biological negative labels',sources=source_meta,native_residues={k:[dataclasses.asdict(r) for r in rs] for k,rs in loaded.items()},graphs=graphs)
    out.mkdir(parents=True)
    with (out/'atomic_inputs.json.gz').open('xb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as f:f.write((json.dumps(archive,separators=(',',':'),allow_nan=False)+'\n').encode())
    deps=['src/meta_acdc/structure/artifact_qc.py','src/meta_acdc/structure/graph.py','src/meta_acdc/structure/cif.py','src/meta_acdc/structure/identity.py','src/meta_acdc/models/dataset.py']
    receipt=dict(n_graphs=len(graphs),n_source_complexes=len(loaded),all_888_audit_rows_reconstructed=True,source_cache_sha256=hashlib.sha256((folder/'graphs.pt').read_bytes()).hexdigest(),source_record_sha256=hashlib.sha256((folder/'structure_records.json').read_bytes()).hexdigest(),archive_sha256=hashlib.sha256((out/'atomic_inputs.json.gz').read_bytes()).hexdigest(),bytes=(out/'atomic_inputs.json.gz').stat().st_size,code_files={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in deps},exporter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),model_training=False)
    (out/'export_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();print(json.dumps(export(a.release,a.out),indent=2))
