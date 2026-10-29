"""Verify weekly public-source discovery without treating failures as negative hits."""
import json,hashlib
from pathlib import Path
B=Path('docs/fingerprint_np105_validation_20261029')
def main():
 result=json.loads((B/'sources/weekly/weekly_TCR.json').read_text());metadata=json.loads((B/'sources/weekly_metadata/weekly_entities.json').read_text());assert not metadata.get('errors');entries=metadata['data']['entries'];assert all(entries);es={e['rcsb_id']:e for e in entries};assert set(es)=={x['identifier'] for x in result['result_set']} and len(es)==41
 plan=json.loads(Path('docs/fingerprint_np_boundaries_20261029/broad_search/search_plan.json').read_text());matches=[]
 for t in plan['targets']:
  for e in entries:
   pairs={role:[x['rcsb_id'] for x in e['polymer_entities'] if seq in x['entity_poly']['pdbx_seq_one_letter_code_can']] for role,seq in t['sequences'].items()}
   if all(pairs.values()) and any(a!=b for a in pairs['alpha'] for b in pairs['beta']):matches.append((t['target_id'],e['rcsb_id']))
 assert not matches
 def chains(pdb):return {x['rcsb_polymer_entity']['pdbx_description'].rsplit(' ',1)[-1]:x['entity_poly']['pdbx_seq_one_letter_code_can'] for x in es[pdb]['polymer_entities'] if x['rcsb_polymer_entity']['pdbx_description'].startswith('TCR1414')}
 assert chains('9XC3')==chains('9X87') and set(chains('9XC3'))=={'alpha','beta'};assert len(es['9X87']['polymer_entities'])==2
 for pdb in ['9XC3','9X9I','9X87']:
  citation=es[pdb]['rcsb_primary_citation'];assert citation['pdbx_database_id_DOI'] is None and citation['pdbx_database_id_PubMed'] is None
 assert len({es[p]['rcsb_primary_citation']['title'].rstrip('.') for p in ['9XC3','9X9I','9X87']})==1
 summary=json.loads((B/'weekly_summary.json').read_text());n=json.loads((B/'weekly_numbering.json').read_text());byname={x['name']:x for x in n['rows']};passed=0;failed=0
 for row in json.loads((B/'sources/weekly_chain_retry/manifest.json').read_text()):
  payload=json.loads((B/'sources/weekly_chain_retry'/row['file']).read_text());assert hashlib.sha256((B/'sources/weekly_chain_retry'/row['file']).read_bytes()).hexdigest()==row['sha256']
  if row['status']==200:
   passed+=1;seq=byname[Path(row['file']).stem]['boundary']['sequence']
   for hit in payload['result_set']:
    entity=next(x for e in entries for x in e['polymer_entities'] if x['rcsb_id']==hit['identifier']);assert entity['entity_poly']['pdbx_seq_one_letter_code_can'].count(seq)==1
  else:failed+=1;assert row['status']==500 and payload['status']==500
 assert (passed,failed)==(1,3);assert summary['n_chain_queries_unresolved']==3 and summary['new_qualified_cross_study_pairs']==0
 assert json.loads((B/'sources/publication/europepmc_title.json').read_text())['hitCount']==0
 receipt=dict(status='pass',weekly_metadata_entries=41,old_query_targets_compared=len(plan['targets']),old_target_exact_pair_matches=0,same_deposited_complete_pair_verified=['9XC3','9X87'],missing_primary_citation_identifiers_verified=True,unbound_two_entity_9X87_verified=True,successful_sequence_queries=passed,unresolved_sequence_queries=failed,failed_queries_not_counted_as_zero_hits=True,successful_hits_exact_substrings=True,new_complete_identity_admissions=0,new_qualified_cross_study_pairs=0,protocol_frozen=False,independent_acceptance_accessed=False)
 (B/'weekly_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
