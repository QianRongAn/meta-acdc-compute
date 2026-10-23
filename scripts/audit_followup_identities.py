"""Overlay newly recovered research chains without changing any frozen cohort."""
from pathlib import Path
import argparse, hashlib, json, sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT))
from meta_acdc.data.experimental_evidence import reserved_hits
from scripts.audit_measured_intake_connections import components

def audit(root):
    d=root/'docs/admission_followup_20261022'
    old=json.loads((root/'docs/measured_intake_20261021/connections.json').read_text())
    new=json.loads((d/'derived/sequence_evidence.json').read_text())
    overrides={(x['pmid'],x['receptor_name']):x for x in new['receptors']}
    nodes=[];changed=[]
    for previous in old['nodes']:
        x=dict(previous);key=(x['pmid'],x['receptor_name'])
        if key in overrides:
            s=overrides[key]
            x.update(alpha_variable=s['alpha_published'],beta_variable=s['beta_published'],
                     cdr3_alpha=s['cdr3_alpha'],cdr3_beta=s['cdr3_beta'],
                     sequence_scope=s['sequence_scope'],sequence_pair_review_pass=False,
                     training_allowed=False,independent_acceptance_allowed=False)
            changed.append(x['receptor_id'])
        nodes.append(x)
    if len(changed)!=6:raise ValueError('Recovered receptor mapping not one-to-one')
    edges=[]
    for e in old['edges']:
        reasons=[r for r in e['reasons'] if r!='corresponding_chain_or_CDR3_near_match']
        if reasons:edges.append(dict(e,reasons=reasons))
    near=[]
    for i,a in enumerate(nodes):
        for b in nodes[:i]:
            if reserved_hits(a,[b]):
                e=dict(a=a['receptor_id'],b=b['receptor_id'],reasons=['corresponding_chain_or_CDR3_near_match'])
                near.append(e);edges.append(e)
    paths=['docs/clinical_primary_20261018/clinical_reservations.json','docs/database_binding/combined_reservation_registry.json']
    registries={}
    for path in paths:
        raw=json.loads((root/path).read_text());registries[path]=raw if isinstance(raw,list) else raw['receptors']
    exposures=[dict(receptor_id=x['receptor_id'],receptor_name=x['receptor_name'],
                    hits={p:reserved_hits(x,v) for p,v in registries.items()},
                    training_allowed=False,independent_acceptance_allowed=False) for x in nodes if x['receptor_id'] in changed]
    blocks=components([x['receptor_id'] for x in nodes],edges)
    families=components([x['receptor_id'] for x in nodes],near)
    inputs=paths+['docs/measured_intake_20261021/connections.json','docs/admission_followup_20261022/derived/sequence_evidence.json']
    return dict(changed_receptor_ids=changed,nodes=nodes,edges=edges,exposure=exposures,
                provisional_joint_blocks=blocks,provisional_sequence_components=families,
                summary=dict(nodes=len(nodes),updated_receptors=len(changed),provisional_joint_blocks=len(blocks),
                             provisional_sequence_components=len(families),training_release_created=False,new_test_created=False),
                rule='Prior conservative source edges retained; sequence edges recomputed using published research chains with leaders. No-hit does not establish family/source independence. Clinical/combined reservation matches are exclusions, not clinical identity equivalence.',
                source_sha256={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in inputs},
                partition_assignment=None,test_labels_used=False,model_scored=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    result=audit(ROOT);a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['summary']))
