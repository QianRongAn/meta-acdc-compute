from pathlib import Path
import json,csv,collections
r=Path(__file__).resolve().parents[1];base=r/'docs/fingerprint_cross_source_20261024';src=base/'sources';entries={}
for p in sorted(src.glob('rcsb_batch_*.json')):
 for x in json.loads(p.read_text())['data']['entries']:entries[x['rcsb_id'].lower()]=x
reg=json.loads((r/'docs/fingerprint_paired_identity_20261024/curation_v3/identity_registry.json').read_text());targets={x['receptor_id']:x for x in reg['rows'] if x['status']=='admitted'}
def study(e):
 c=e.get('rcsb_primary_citation') or {};doi=c.get('pdbx_database_id_DOI');pmid=c.get('pdbx_database_id_PubMed');return ('doi:'+doi.lower()) if doi else ('pmid:'+str(pmid)) if pmid else None
known=[]
for rid,t in targets.items():
 hits=[]
 for pdb,e in entries.items():
  entities=e.get('polymer_entities') or []
  a=[x for x in entities if x.get('entity_poly') and t['alpha']['variable'] in (x['entity_poly'].get('pdbx_seq_one_letter_code_can') or '')]
  b=[x for x in entities if x.get('entity_poly') and t['beta']['variable'] in (x['entity_poly'].get('pdbx_seq_one_letter_code_can') or '')]
  if any(x['rcsb_id']!=y['rcsb_id'] for x in a for y in b):hits.append(dict(pdb=pdb,study=study(e),alpha_entities=[x['rcsb_id'] for x in a],beta_entities=[x['rcsb_id'] for x in b]))
 known.append(dict(receptor_id=rid,anchor=t['complex_id'],hits=hits,n_pdb=len(hits),n_studies=len({x['study'] for x in hits if x['study']})))
print('KNOWN EXACT PAIRS')
for x in known:
 if x['n_pdb']>1:print(x['anchor'],x['n_pdb'],x['n_studies'],[x['pdb'] for x in x['hits']])
(base/'known_pair_search_leads.json').write_bytes((json.dumps(known,indent=2)+'\n').encode())
rows=list(csv.DictReader((src/'tcr_chain_data.tsv').read_text().splitlines(),delimiter='\t'));byname=collections.defaultdict(set);bycontext=collections.defaultdict(set)
for x in rows:
 if x['complex_type']=='CLASSI' and x['tcr_name']:byname[x['tcr_name']].add(x['pdb_id'])
leadnames={name:sorted(pdbs) for name,pdbs in byname.items() if len({study(entries[p]) for p in pdbs if p in entries and study(entries[p])})>1}
# All class-I entries sharing both V genes with one already admitted pair are search candidates, never labels.
by_pdb=collections.defaultdict(list)
for x in rows:by_pdb[x['pdb_id']].append(x)
context_targets={(x['alpha']['v_allele'].split('*')[0],x['beta']['v_allele'].split('*')[0]) for x in targets.values()}
context_leads=[]
for pdb,rr in by_pdb.items():
 for a in rr:
  for b in rr:
   if a['tcr_type']=='Alpha' and b['tcr_type']=='Beta' and a['complex_type']=='CLASSI' and b['complex_type']=='CLASSI' and (a['V_gene'],b['V_gene']) in context_targets:
    context_leads.append(pdb)
selected=sorted(set(pdb for values in leadnames.values() for pdb in values)|set(context_leads)|{h['pdb'] for x in known for h in x['hits']})
(base/'candidate_selection.json').write_bytes((json.dumps(dict(scope='search leads only; names and genes are never identity truth',multi_study_name_leads=leadnames,known_V_pair_leads=sorted(set(context_leads)),selected_pdbs=selected),indent=2)+'\n').encode())
print('MULTI_STUDY_NAME_LEADS',len(leadnames),'PDBS',len(set(p for x in leadnames.values() for p in x)),'KNOWN_CONTEXT_PDBS',len(set(context_leads)),'SELECTED',len(selected));print(list(leadnames.items())[:25])
