"""Require both ends of an expanded cross-study lead to pass unchanged fixed QC."""
import json
from pathlib import Path
import sys
from scripts import validate_native_search

def complete_pair_ready(anchor_pass,candidate_pass,complete):
    return bool(anchor_pass and candidate_pass and complete)

def main():
    validate_native_search.BASE=Path(__file__).resolve().parents[1]/'docs/fingerprint_expanded_controls_20261029/complete_search'
    validate_native_search.main()
    out=Path(sys.argv[sys.argv.index('--out')+1]);report=json.loads(out.read_text())
    anchors={r['target_id']:any(p['existing_fixed_QC_pass'] for p in r['chain_pair_checks']) for r in report['rows'] if r['pdb']==r['anchor']}
    for r in report['rows']:
        r['anchor_existing_fixed_QC_pass']=anchors.get(r['target_id'],False)
        r['both_ends_fixed_QC_pass']=complete_pair_ready(r['anchor_existing_fixed_QC_pass'],any(p['existing_fixed_QC_pass'] for p in r['chain_pair_checks']),r['query_is_complete_identity'])
        if r['decision']=='requires_new_protocol_readiness_review' and not r['both_ends_fixed_QC_pass']:r['decision']='candidate_QC_pass_but_anchor_fails_fixed_QC'
    report['n_new_qualified_cross_study_leads']=sum(r['cross_study'] and r['both_ends_fixed_QC_pass'] for r in report['rows'])
    report['complete_cross_study_leads']=[r for r in report['rows'] if r['query_is_complete_identity'] and r['cross_study'] and r['full_query_pair_on_distinct_entities']]
    report['interpretation']+=' Both anchor and partner must pass: a partner alone never supplies a qualified cross-study pair.'
    out.write_text(json.dumps(report,indent=2)+'\n');print('BOTH_ENDS_QUALIFIED',report['n_new_qualified_cross_study_leads'])
if __name__=='__main__':main()
