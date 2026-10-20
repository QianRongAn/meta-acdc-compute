"""Read-only reconstruction of historical release atom defects; no training."""
from __future__ import annotations
import argparse, collections, json
from pathlib import Path
import torch
from meta_acdc.structure.artifact_qc import load_record_residues, construction_screen
from meta_acdc.structure.identity import missing_heavy_atoms
from meta_acdc.models.dataset import make_decoy, make_displaced_peptide_decoy
from meta_acdc.data.release import sha256


def audit(folder):
    data=torch.load(folder/'graphs.pt',map_location='cpu',weights_only=False)
    meta={r['complex_id']:r for r in json.loads((folder/'structure_records.json').read_text())}
    loaded={};rows=[];by_kind=collections.defaultdict(collections.Counter)
    def get(identifier):
        if identifier not in loaded: loaded[identifier]=load_record_residues(meta[identifier],strict_identity=False)
        return loaded[identifier]
    for row,keys in zip(data['records'],data['node_keys']):
        rs=get(row['complex_id'])
        if row['kind']=='graft':
            donor=[r for r in get(row['donor']) if r.role=='peptide']
            rs=make_decoy(rs,donor)
        elif row['kind']=='displaced':rs=make_displaced_peptide_decoy(rs)
        selected={(c,n,i) for c,n,i in keys}
        missing=[dict(chain=r.chain,resid=r.resid,icode=r.icode,role=r.role,atoms=missing_heavy_atoms(r))
                 for r in rs if (r.chain,r.resid,r.icode) in selected and missing_heavy_atoms(r)]
        screen=construction_screen(rs)
        from meta_acdc.structure.identity import chain_sequence
        base=get(row['complex_id']); declared=meta[row['complex_id']]
        mismatches=[field for field,seqkey in [('alpha_chain','alpha'),('beta_chain','beta'),('mhc_chain','mhc_sequence'),('b2m_chain','b2m_sequence'),('peptide_chain','peptide')]
                    if chain_sequence([r for r in base if r.chain==declared[field]])!=declared[seqkey]]
        result=dict(graph_id=row['graph_id'],kind=row['kind'],domain=row['domain'],
                    family_id=row['family_id'],fold=row['fold'],declared_source_sequence_mismatches=mismatches,missing_nodes=missing,
                    construction=screen)
        rows.append(result);stats=by_kind[row['kind']];stats['graphs']+=1
        stats['declared_sequence_mismatch_graphs']+=bool(mismatches)
        stats['graphs_with_incomplete_interface']+=bool(missing)
        stats['incomplete_interface_nodes']+=len(missing)
        stats['graphs_failing_construction_screen']+=not screen['passed']
        for m in missing:stats['missing_'+m['role']]+=1
        if len(rows)%100==0:print('audited',len(rows),flush=True)
    return dict(source_manifest_sha256=sha256(folder/'release_manifest.json'),
                read_only=True, threshold_is_screening_not_physical_validation=True,
                summary={k:dict(v) for k,v in by_kind.items()}, rows=rows)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--release-dir',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    if args.out.exists():raise ValueError('refusing to overwrite audit')
    report=audit(args.release_dir);args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['summary'],indent=2))


if __name__=='__main__':main()
