"""Expand exact paired-chain discovery to all 44 complete identities; unchanged admission gates."""
from pathlib import Path
import argparse,hashlib,json
from urllib.parse import quote
ROOT=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    paths=[ROOT/'docs/fingerprint_positive_controls_20261027/source_balanced_design.json',ROOT/'docs/fingerprint_construct_followup_20261027/construct_adjudication_final.json']
    rows={}
    for p in paths:
        for r in json.loads(p.read_text())['rows']: rows[(r['pdb'],r['alpha_chain'],r['beta_chain'])]=r
    reps={}
    for r in rows.values():
        if r['status']=='identity_admitted_structure_QC_pending':reps.setdefault(r['receptor_id'],r)
    assert len(reps)==44
    targets=[dict(target_id=rid,anchor=r['pdb'],anchor_study=r['source_study_key'],full_boundary_validated=True,sequences={k:r[k]['variable'] for k in ['alpha','beta']}) for rid,r in sorted(reps.items())]
    unique=sorted({(role,seq) for t in targets for role,seq in t['sequences'].items()})
    queries=[];requests=[]
    for i,(role,sequence) in enumerate(unique):
        name=f'q{i:03d}.json';q=dict(query=dict(type='terminal',service='sequence',parameters=dict(evalue_cutoff=.1,identity_cutoff=1.,sequence_type='protein',value=sequence)),return_type='polymer_entity',request_options=dict(return_all_hits=True))
        queries.append(dict(file=name,role=role,sequence=sequence,query=q));requests.append([name,'https://search.rcsb.org/rcsbsearch/v2/query?json='+quote(json.dumps(q,separators=(',',':')))])
    plan=dict(targets=targets,queries=queries,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},scope='44 complete paired variable identities including anchors outside fixed atomic QC; discovery does not admit candidates',protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    (args.out/'search_plan.json').write_text(json.dumps(plan,indent=2)+'\n');(args.out/'requests.json').write_text(json.dumps(requests,indent=2)+'\n');print(len(targets),len(queries))
if __name__=='__main__':main()
