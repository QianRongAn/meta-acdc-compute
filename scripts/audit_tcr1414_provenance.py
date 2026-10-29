"""Separate receptor-generation provenance from structural study identity; no admission changes."""
from pathlib import Path
import json,hashlib,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'docs/fingerprint_np105_recovery_20261029';OLD=ROOT/'docs/fingerprint_np105_validation_20261029'
FIGURE_CHAINS={'1A':dict(role='alpha',V='TRAV21',J='TRAJ28',CDR3='CAVKSGAGSYQLTF'),'2A':dict(role='alpha',V='TRAV22',J='TRAJ29',CDR3='CAAPSGNTPLVF'),'1B':dict(role='beta',V='TRBV18',J='TRBJ2-5',CDR3='CASSPFQETQYF'),'3B':dict(role='beta',V='TRBV28',J='TRBJ2-6',CDR3='CASRVSGANVLTF')}
def main():
 path=B/'sources/primary/tcr1414_primary.xml';root=ET.fromstring(path.read_bytes());doi=next(x.text for x in root.iter('article-id') if x.get('pub-id-type')=='doi');assert doi=='10.1136/jitc-2022-006001';paragraphs=[''.join(p.itertext()) for p in root.iter('p')]
 anchors={'separate_chain_RACE':'clones per TCR chain','combinatorial_recovery':'15 different','pair_naming':'called TCR1414_1 and TCR1414_2','humanized_mouse_origin':'ABabDR4 mice'}
 evidence={}
 for name,anchor in anchors.items():
  texts=[text for text in paragraphs if anchor in text];assert texts,name;evidence[name]=dict(matching_paragraph_count=len(texts),paragraph_sha256=[hashlib.sha256(x.encode()).hexdigest() for x in texts])
 entries={e['rcsb_id']:e for e in json.loads((OLD/'sources/weekly_metadata/weekly_entities.json').read_text())['data']['entries']};rows=[]
 for pdb,pair in [('9XC3',['1A','1B']),('9X9I',['2A','3B'])]:
  chains=[]
  for label in pair:
   fact=FIGURE_CHAINS[label];entities=[x for x in entries[pdb]['polymer_entities'] if x['rcsb_polymer_entity']['pdbx_description'].endswith(fact['role']) and fact['CDR3'] in x['entity_poly']['pdbx_seq_one_letter_code_can']];assert len(entities)==1;chains.append(dict(figure_chain=label,**fact,entity=entities[0]['rcsb_id'],CDR3_exact_unique_entity_substring=entities[0]['entity_poly']['pdbx_seq_one_letter_code_can'].count(fact['CDR3'])==1))
  citation=entries[pdb]['rcsb_primary_citation'];assert citation['pdbx_database_id_DOI'] is None;rows.append(dict(pdb=pdb,receptor_generation_DOI=doi,structure_publication_DOI=None,chains=chains,native_single_cell_pairing_established_by_reported_method=False,decision='not admitted as a verified independent native paired receptor; functional combinatorial recovery is not same-cell pairing evidence',full_variable_sequence_identity_across_papers_established=False))
 numbered={x['name']:x for x in json.loads((OLD/'weekly_numbering.json').read_text())['rows']};results={};checked_queries=[]
 for name in numbered:
  source=OLD/'sources/weekly_chain_retry'/(name+'.json') if name=='9XC3_5' else B/'sources/retry'/(name+'.json');response=json.loads(source.read_text());assert 'result_set' in response
  hits=response['result_set'];seq=numbered[name]['boundary']['sequence'];entities=[]
  for hit in hits:
   entry=entries[hit['identifier'].split('_')[0]];entity=next(x for x in entry['polymer_entities'] if x['rcsb_id']==hit['identifier']);assert entity['entity_poly']['pdbx_seq_one_letter_code_can'].count(seq)==1;entities.append(hit['identifier'])
  results[name]={x.split('_')[0] for x in entities};checked_queries.append(dict(query=name,entities=entities,full_interval_exact_substrings=True))
 paired={'TCR1414_1':sorted(results['9XC3_4']&results['9XC3_5']),'TCR1414_2':sorted(results['9X9I_4']&results['9X9I_5'])};assert paired=={'TCR1414_1':['9X87','9XC3'],'TCR1414_2':['9X9I']}
 record=dict(primary_article=dict(DOI=doi,PMCID='PMC9950979',PMID='36822673',source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),license='CC BY-NC 4.0; original retained locally'),methods_evidence=evidence,reported_chain_candidates=dict(alpha=3,beta=5,tested_combinations=15),figure3_facts=FIGURE_CHAINS,figure_transcription_manually_visually_reviewed=True,figure_sha256=hashlib.sha256((B/'sources/figure/figure3.jpg').read_bytes()).hexdigest(),receptors=rows,all_four_chain_queries_successful=True,checked_queries=checked_queries,paired_PDB_hits=paired,receptor_generation_DOI_not_imputed_as_structure_DOI=True,n_new_complete_identity_admissions=0,n_new_qualified_cross_study_pairs=0,n_new_verified_same_source_same_V_native_groups=0,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
 (B/'tcr1414_source_audit.json').write_text(json.dumps(record,indent=2)+'\n');print(paired,'NATIVE_ADMISSIONS',0)
if __name__=='__main__':main()
