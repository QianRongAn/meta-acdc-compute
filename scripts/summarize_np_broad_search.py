"""Summarize all broadened cross-study discovery clusters without upgrading their labels."""
from pathlib import Path
import argparse,collections,gzip,json,hashlib
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_boundaries_20261029'
REASONS={
'3huj':'CD1d/NKT complex; outside complete classical class-I peptide scope',
'3tfk':'42F3 compact MHC without beta2m; outside complete class-I scope',
'6pvc':'MR1/MAIT ligand context; outside classical class-I peptide scope',
'2ckb':'existing full numbered alpha N-terminus quarantine; endpoint-only substring variants cannot repair it',
'7n1d':'existing full numbered beta N-terminus quarantine; 7RTR also lacks interface heavy atoms',
'9bbc':'TCR-CD3 detergent/nanodisc complexes without required pMHC; no qualifying paired class-I unit',
'1kb5':'existing alpha N-terminus quarantine and single-chain fusion context; no paired-chain admission',
'5fka':'existing alpha complete J/FR4 terminal quarantine; endpoint-only alternative boundary is not verified numbering',
'2uwe':'existing full numbered alpha N-terminus quarantine; endpoint-only substring variants cannot repair it',
'5l2k':'GEM42 CD1b/apo context; outside complete classical class-I peptide scope',
'5vcj':'CD1d/iNKT context; outside complete classical class-I peptide scope'}
def report_bytes():
    path=BASE/'broad_search/primary_validation_final.json'
    return path.read_bytes() if path.exists() else gzip.decompress(path.with_suffix('.json.gz').read_bytes())
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    r=json.loads(report_bytes());g=collections.defaultdict(list)
    for row in r['discovery_only_cross_study_leads']:g[row['anchor']].append(row)
    assert set(g)==set(REASONS)
    old=json.loads((ROOT/'docs/fingerprint_cross_source_extension_20261024/sequence_curation_all_class_I.json').read_text());oldrows=collections.defaultdict(list)
    for x in old['rows']:oldrows[x['pdb']].append(x)
    groups=[]
    for anchor,rows in sorted(g.items()):groups.append(dict(anchor=anchor,n_discovery_variable_pair_hypotheses=len({x['target_id'] for x in rows}),target_ids=sorted({x['target_id'] for x in rows}),partner_pdbs=sorted({x['pdb'] for x in rows}),n_target_partner_leads=len(rows),scope_or_unresolved_identity_reason=REASONS[anchor],existing_annotation_quarantine_reasons=sorted({s for x in oldrows.get(anchor,[]) for s in x['reasons']}),admitted=False,atomic_QC_newly_run=False))
    oldids={e['rcsb_id'].lower() for p in (ROOT/'docs/fingerprint_cross_source_20261024/sources').glob('rcsb_batch_*.json') for e in json.loads(p.read_text())['data']['entries'] if e}
    hits=json.loads((BASE/'broad_search/search_results.json').read_text());plan=json.loads((BASE/'broad_search/search_plan.json').read_text())
    out=dict(n_targets=len(plan['targets']),n_previously_complete_identity_targets=sum(x['full_boundary_validated'] for x in plan['targets']),n_discovery_only_targets=sum(not x['full_boundary_validated'] for x in plan['targets']),n_chain_queries=hits['n_queries'],n_successful_queries=hits['n_successful_queries'],n_fresh_metadata_PDBs=r['n_fresh_metadata_pdbs'],PDBs_outside_old_799_archive=sorted(set(hits['discovery_pdbs'])-oldids),n_discovery_only_cross_study_leads=r['n_discovery_only_cross_study_leads'],n_discovery_only_variable_pair_hypotheses=len({x['target_id'] for x in r['discovery_only_cross_study_leads']}),n_anchor_PDB_clusters=len(groups),discovery_only_clusters=groups,n_previously_complete_cross_study_leads=len(r['complete_cross_study_leads']),n_new_qualified_cross_study_pairs=r['n_new_qualified_cross_study_leads'],interpretation='Hypotheses can differ only by candidate boundaries and can share an anchor. 197 leads/27 hypotheses/11 anchor PDB clusters are not independent receptors, accepted labels or positives.',protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    a.out.write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k.startswith('n_')})
if __name__=='__main__':main()
