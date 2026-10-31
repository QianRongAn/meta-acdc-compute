"""Find complete-construct rescue leads using unchanged observed fragments only."""
from pathlib import Path
import hashlib,json
from urllib.parse import urlencode
from meta_acdc.structure.paired_variable_identity import numbered_domains
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'docs/fingerprint_control_closure_20261101'
def main():
 out=BASE/'fragment_rescue';out.mkdir(exist_ok=False)
 old=ROOT/'docs/fingerprint_positive_controls_20261027/source_balanced_design.json'
 rows=json.loads(old.read_text())['rows'];oldplan=json.loads((ROOT/'docs/fingerprint_np_boundaries_20261029/broad_search/search_plan.json').read_text())
 prior={(role,seq) for t in oldplan['targets'] for role,seq in t['sequences'].items()};targets=[]
 for pdb in ['1ao7','3mv7','1bd2','1mi5']:
  row=next(x for x in rows if x['pdb']==pdb)
  path=ROOT/'docs/fingerprint_cross_source_extension_20261024/sources/stcrdab_all'/(pdb+'.html');ann=numbered_domains(path.read_text());seqs={role:ann[row[role+'_chain']]['sequence'] for role in ['alpha','beta']}
  for role,seq in seqs.items():assert row[role]['full_construct'].count(seq)==1
  targets.append(dict(anchor=pdb,sequences=seqs,annotation_file=str(path.relative_to(ROOT)),annotation_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),unchanged_quarantine=row['reasons'],anchor_study=row['source_study_key'],complete_identity_claimed=False,previous_exact_queries={role:(role,seq) in prior for role,seq in seqs.items()}))
 queries=[];requests=[]
 for i,(role,seq) in enumerate(sorted({(role,seq) for t in targets for role,seq in t['sequences'].items()})):
  query=dict(query=dict(type='terminal',service='sequence',parameters=dict(evalue_cutoff=.1,identity_cutoff=1.0,sequence_type='protein',value=seq)),return_type='polymer_entity',request_options=dict(return_all_hits=True));name=f'fragment_{i:02d}.json';queries.append(dict(file=name,role=role,sequence=seq,query=query));requests.append([name,'https://search.rcsb.org/rcsbsearch/v2/query?'+urlencode({'json':json.dumps(query,separators=(',',':'))})])
 plan=dict(date='2026-11-01',targets=targets,queries=queries,interpretation='Existing exact deposited fragments only; sequence search identity 1 is local alignment, not full variable identity. New full-construct ends, explicit pairing, source and both-end atomic gates remain mandatory.',independent_acceptance_accessed=False,protocol_frozen=False)
 (out/'plan.json').write_text(json.dumps(plan,indent=2)+'\n');(out/'requests.json').write_text(json.dumps(requests,indent=2)+'\n');print([(t['anchor'],t['previous_exact_queries']) for t in targets])
if __name__=='__main__':main()
