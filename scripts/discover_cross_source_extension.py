from pathlib import Path
import json,hashlib,csv,re,collections,subprocess
from meta_acdc.structure.paired_variable_identity import germline_references,pair_key
r=Path(__file__).resolve().parents[1];old=r/'docs/fingerprint_cross_source_20261024';d=r/'docs/fingerprint_cross_source_extension_20261024';d.mkdir(exist_ok=True)
refs=germline_references(list((r/'docs/fingerprint_paired_identity_20261024/sources').glob('Homo_sapiens_TR*.fasta'))+list((r/'docs/fingerprint_paired_identity_20261024/sources').glob('Mus_musculus_TR*.fasta')))
entries={e['rcsb_id'].lower():e for p in (old/'sources').glob('rcsb_batch_*.json') for e in json.loads(p.read_text())['data']['entries'] if e}
rows=list(csv.DictReader((old/'sources/tcr_chain_data.tsv').read_text().splitlines(),delimiter='\t'));by_pdb=collections.defaultdict(list)
for x in rows:by_pdb[x['pdb_id'].lower()].append(x)
def candidates(seq,role):
 out=[];pref='TRAV' if role=='Alpha' else 'TRBV';jp='TRAJ' if role=='Alpha' else 'TRBJ'
 for sp in ['human','mouse']:
  vs=[(g,x['amino_acids'][:10]) for (s,g),values in refs.items() if s==sp and g.startswith(pref) for x in values]
  js=set(x['amino_acids'][-10:] for (s,g),values in refs.items() if s==sp and g.startswith(jp) for x in values)
  for g,v in vs:
   start=seq.find(v)
   if start<0 or start>30:continue
   for j in js:
    stop=seq.find(j,start+70)
    if stop<0:continue
    end=stop+10
    if 80<=end-start<=140:out.append((sp,g.split('*')[0],seq[start:end]))
 return sorted(set(out))
leads=[]
for pdb,e in entries.items():
 entities={c:z for z in e.get('polymer_entities') or [] for c in (z.get('rcsb_polymer_entity_container_identifiers') or {}).get('auth_asym_ids',[])}
 variants=collections.defaultdict(lambda:collections.defaultdict(list))
 for x in by_pdb[pdb]:
  if x['tcr_type'] not in ['Alpha','Beta']:continue
  for chain in x['cdr3_chains'].split(','):
   chain=chain.strip()
   if chain not in entities:continue
   z=entities[chain];seq=z['entity_poly']['pdbx_seq_one_letter_code_can'];cs=candidates(seq,x['tcr_type'])
   for sp,v,variable in cs:
    if x['cdr3_sequences'] and x['cdr3_sequences'] not in variable:continue
    variants[x['tcr_name']][x['tcr_type']].append(dict(chain=chain,species=sp,V=v,variable=variable))
 citation=e.get('rcsb_primary_citation') or {};study=(citation.get('pdbx_database_id_DOI') or '').lower() or str(citation.get('pdbx_database_id_PubMed') or '')
 for name,roles in variants.items():
  for a in roles['Alpha']:
   for b in roles['Beta']:
    if a['species']!=b['species']:continue
    leads.append(dict(pdb=pdb,name=name,alpha=a,beta=b,receptor_id=pair_key(a['variable'],b['variable']),V_context='|'.join([a['species'],a['V'],b['V']]),study=study,discovery_only=True))
by_id=collections.defaultdict(list)
for x in leads:by_id[x['receptor_id']].append(x)
cross={k:v for k,v in by_id.items() if len({x['study'] for x in v if x['study']})>1}
contexts={x['V_context'] for group in cross.values() for x in group};studies={x['study'] for group in cross.values() for x in group}
selected=sorted({x['pdb'] for x in leads if x['receptor_id'] in cross or x['V_context'] in contexts or x['study'] in studies})
oldsel=set(json.loads((old/'candidate_selection_with_quarantine_followup.json').read_text())['selected_pdbs'])
(d/'discovery_leads.json').write_text(json.dumps(dict(scope='candidate discovery only, never identity admission; explicit STCRDab pairing still required',n_archive_entries=len(entries),rows=leads,cross_study_groups=cross),indent=2)+'\n')
(d/'candidate_selection.json').write_text(json.dumps(dict(selected_pdbs=selected,new_pdbs=sorted(set(selected)-oldsel),n_old_overlap=len(set(selected)&oldsel),independent_acceptance_accessed=False),indent=2)+'\n')
print('LEADS',len(leads),'CROSS_GROUPS',len(cross),'SELECTED',len(selected),'NEW',len(set(selected)-oldsel))
for group in cross.values():print(group[0]['name'],group[0]['V_context'],sorted({x['pdb'] for x in group}),sorted({x['study'] for x in group}))
