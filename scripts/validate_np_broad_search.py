"""Validate broad public sequence leads; preserve incomplete identity and both-end QC requirements."""
import json,sys
from pathlib import Path
from scripts import validate_native_search
from scripts.validate_expanded_controls import complete_pair_ready
ROOT=Path(__file__).resolve().parents[1]
def main():
    validate_native_search.BASE=ROOT/'docs/fingerprint_np_boundaries_20261029/broad_search'
    validate_native_search.main()
    p=Path(sys.argv[sys.argv.index('--out')+1]);r=json.loads(p.read_text())
    anchors={x['target_id']:any(y['existing_fixed_QC_pass'] for y in x['chain_pair_checks']) for x in r['rows'] if x['pdb']==x['anchor']}
    for x in r['rows']:
        x['anchor_existing_fixed_QC_pass']=anchors.get(x['target_id'],False)
        x['both_ends_fixed_QC_pass']=complete_pair_ready(x['anchor_existing_fixed_QC_pass'],any(y['existing_fixed_QC_pass'] for y in x['chain_pair_checks']),x['query_is_complete_identity'])
        if x['decision']=='requires_new_protocol_readiness_review' and not x['both_ends_fixed_QC_pass']:x['decision']='candidate_QC_pass_but_anchor_fails_fixed_QC'
    r['n_new_qualified_cross_study_leads']=sum(x['cross_study'] and x['both_ends_fixed_QC_pass'] for x in r['rows'])
    r['complete_cross_study_leads']=[x for x in r['rows'] if x['query_is_complete_identity'] and x['cross_study'] and x['full_query_pair_on_distinct_entities']]
    r['discovery_only_cross_study_leads']=[x for x in r['rows'] if not x['query_is_complete_identity'] and x['cross_study'] and x['full_query_pair_on_distinct_entities']]
    r['n_discovery_only_cross_study_leads']=len(r['discovery_only_cross_study_leads'])
    p.write_text(json.dumps(r,indent=2)+'\n');print('FINAL',r['n_new_qualified_cross_study_leads'],'DISCOVERY_ONLY',r['n_discovery_only_cross_study_leads'])
if __name__=='__main__':main()
