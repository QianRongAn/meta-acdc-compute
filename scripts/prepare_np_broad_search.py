"""Replay broadened paired-chain discovery plan without reading independent acceptance data."""
from pathlib import Path
import argparse,hashlib,json
from urllib.parse import quote
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_boundaries_20261029'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    p=ROOT/'docs/fingerprint_cross_source_extension_20261024/discovery_leads.json';d=json.loads(p.read_text());reps={}
    for r in d['rows']:reps.setdefault(r['receptor_id'],r)
    strict=json.loads((ROOT/'docs/fingerprint_expanded_controls_20261029/complete_search/search_plan.json').read_text())
    targets={r['target_id']:dict(r,scope='previous complete fixed identity anchor') for r in strict['targets']}
    for rid,r in sorted(reps.items()):
        if rid not in targets:targets[rid]=dict(target_id=rid,anchor=r['pdb'],anchor_study=('doi:'+r['study'] if r['study'].startswith('10.') else ('pmid:'+r['study'] if r['study'] else None)),full_boundary_validated=False,sequences={k:r[k]['variable'] for k in ['alpha','beta']},scope='functional V first10/J last10 discovery from old public archive; never complete identity admission')
    raw=json.loads((BASE/'numbering/numbering_raw.json').read_text());rs=raw['results']
    for r in json.loads((BASE/'final_primary_cells.json').read_text())['rows']:
        seqs={}
        for role in ['alpha','beta']:
            q=rs[r['name']+'_'+role]
            if q.get('error') or q.get('chain_type')!={'alpha':'A','beta':'B'}[role]:break
            seqs[role]=''.join(aa for _,aa in q['numbering'] if aa!='-')
        if len(seqs)==2:
            rid=hashlib.sha256((seqs['alpha']+'|'+seqs['beta']).encode()).hexdigest();targets.setdefault(rid,dict(target_id=rid,anchor=r['pdb'],anchor_study='doi:10.1016/j.celrep.2026.117908',full_boundary_validated=False,sequences=seqs,scope='ANARCII numbered boundaries; independent annotation only; no admission'))
    queries=[];req=[]
    for i,(role,seq) in enumerate(sorted({(k,s) for r in targets.values() for k,s in r['sequences'].items()})):
        q=dict(query=dict(type='terminal',service='sequence',parameters=dict(evalue_cutoff=.1,identity_cutoff=1.,sequence_type='protein',value=seq)),return_type='polymer_entity',request_options=dict(return_all_hits=True));name=f'q{i:03d}.json';queries.append(dict(file=name,role=role,sequence=seq,query=q));req.append([name,'https://search.rcsb.org/rcsbsearch/v2/query?json='+quote(json.dumps(q,separators=(',',':')))])
    (a.out/'search_plan.json').write_text(json.dumps(dict(targets=list(targets.values()),queries=queries,n_previous_discovery_IDs=len(reps),source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()},scope='broad public development discovery only; functional endpoints/numbering do not grant admission',protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False),indent=2)+'\n')
    (a.out/'requests.json').write_text(json.dumps(req,indent=2)+'\n');print(len(targets),len(queries))
if __name__=='__main__':main()
