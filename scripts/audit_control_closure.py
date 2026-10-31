"""Recheck fragment discoveries and source gaps without relaxing complete identity gates."""
from pathlib import Path
import argparse,hashlib,json
from scripts.validate_native_search import exact_entity_pairs
from scripts.curate_cross_source_receptors import chain_pairs
from scripts.curate_cross_source_receptors import study_key
from scripts.audit_np_boundaries import checked_numbering,raw_j_references
from meta_acdc.structure.paired_variable_identity import numbered_domains,germline_references,validate_domain,pair_key
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'docs/fingerprint_control_closure_20261101'
def sha(data):return hashlib.sha256(data).hexdigest()
def checked_sources(folder):
 rows=json.loads((folder/'manifest.json').read_text())
 for row in rows:
  p=folder/row['file'];raw=p.read_bytes() if p.exists() else b''
  assert len(raw)==row['bytes'] and sha(raw)==row['sha256'],str(p)
  assert row['returncode']==0 and row['status'] in ['200','204'],row['file']
 return rows
def absent_numbering(path):
 text=path.read_text();return dict(domains=len(numbered_domains(text)),explicit_pairs=chain_pairs(text),explicit_missing_entry='does not exist in our database' in text)
def numbering_diagnosis():
 folder=BASE/'4PRH_numbering';raw=json.loads((folder/'raw.json').read_text());summary=json.loads((folder/'summary.json').read_text());plan=json.loads((BASE/'4PRH_numbering_plan.json').read_text());assert raw['inputs']==plan['inputs'] and sha((BASE/'4PRH_numbering_plan.json').read_bytes())==raw['input_plan_sha256']
 refs_dir=ROOT/'docs/fingerprint_paired_identity_20261024/sources';refs=germline_references(list(refs_dir.glob('Homo_sapiens_TR*.fasta')));allj=raw_j_references(list(refs_dir.glob('*_TRAJ.fasta'))+list(refs_dir.glob('*_TRBJ.fasta')));rows=[]
 for i,x in enumerate(raw['inputs']):
  expected='A' if x['role']=='alpha' else 'B';n=raw['neural'][x['name']];mapping,start,end=raw['classic'][i][0];a=checked_numbering(x['sequence'],mapping,start,end,raw['classic_details'][i][0]['chain_type'],expected);b=checked_numbering(x['sequence'],n['numbering'],n['query_start'],n['query_end'],n['chain_type'],expected);assert a==b
  prior=next(y for y in summary['rows'] if y['name']==x['name']);assert a==json.loads(json.dumps(prior['methods']['ANARCI'])) and prior['exact_two_method_position_agreement']
  seq=a['sequence'];vp='TRAV' if expected=='A' else 'TRBV';jp='TRAJ' if expected=='A' else 'TRBJ';v=sorted(g for (sp,g),rr in refs.items() if sp=='human' and g.startswith(vp) and any(seq.startswith(y['amino_acids'][:10]) for y in rr));j=sorted(g for (sp,g),rr in refs.items() if sp=='human' and g.startswith(jp) and any(seq.endswith(y['amino_acids'][-10:]) for y in rr));assert v==prior['functional_V_first10_matches'] and j==prior['functional_J_last10_matches']
  rows.append(dict(name=x['name'],role=x['role'],sequence_sha256=sha(x['sequence'].encode()),numbered_interval=[a['start0'],a['end0_exclusive']],observed_variable=seq,first10=seq[:10],last10=seq[-10:],last_IMGT=a['last_IMGT'],two_method_positions_exact=True,functional_human_V_matches=v,functional_human_J_matches=j,all_function_status_human_or_mouse_J_exact_matches=[dict(gene=y['gene'],function=y['function'],reference_header=y['header']) for y in allj if y['gene'].startswith(jp) and seq.endswith(y['sequence'][-10:])],identity_admitted=False))
 return dict(rows=rows,inputs_trimmed=False,imputation=False,raw_output_sha256=sha((folder/'raw.json').read_bytes()),inference_not_replayed_by_this_audit=True,external_STCRDab_still_missing=True,identity_admitted=False)
def audit():
 f=BASE/'fragment_rescue';plan=json.loads((f/'plan.json').read_text());sources=checked_sources(f/'sources');byfile={x['file']:x for x in sources};hits={}
 for q in plan['queries']:
  assert q['query']['query']['parameters']['value']==q['sequence']
  data=json.loads((f/'sources'/q['file']).read_text());assert data['result_type']=='polymer_entity'
  hits[(q['role'],q['sequence'])]={x['identifier'].split('_')[0] for x in data['result_set']}
 metadata={}
 for s in checked_sources(f/'metadata'):
  data=json.loads((f/'metadata'/s['file']).read_text());assert not data.get('errors')
  metadata.update({e['rcsb_id']:e for e in data['data']['entries']})
 checked_sources(f/'annotation');checked_sources(f/'fasta');refs_dir=ROOT/'docs/fingerprint_paired_identity_20261024/sources';refs=germline_references(list(refs_dir.glob('Homo_sapiens_TR*.fasta'))+list(refs_dir.glob('Mus_musculus_TR*.fasta')));oldrows=json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/source_balanced_design.json').read_text())['rows'];known={x['pdb'].upper() for x in oldrows};records=[];discovery=[]
 for t in plan['targets']:
  source=ROOT/t['annotation_file'];assert sha(source.read_bytes())==t['annotation_sha256']
  seqs=t['sequences'];pdbs=sorted(hits[('alpha',seqs['alpha'])]&hits[('beta',seqs['beta'])]);discovery.extend(pdbs)
  for pdb in pdbs:
   e=metadata[pdb];pairs=exact_entity_pairs(e,seqs['alpha'],seqs['beta']);stored=[x for x in oldrows if x['pdb']==pdb.lower()]
   row=dict(anchor=t['anchor'],pdb=pdb,exact_distinct_entity_fragment_pair=bool(pairs),local_alignment_only=not bool(pairs),study=study_key(e.get('rcsb_primary_citation') or {}),different_paper=bool(study_key(e.get('rcsb_primary_citation') or {}) and study_key(e.get('rcsb_primary_citation') or {})!=t['anchor_study']),title=e['struct']['title'],molecules=[dict(entity=x['rcsb_id'],description=x['rcsb_polymer_entity']['pdbx_description']) for x in e['polymer_entities']],previous_identity_quarantine=sorted({reason for x in stored for reason in x['reasons']}),not_in_existing_candidate_identity_rows=pdb not in known,complete_identity_claimed=False,identity_admitted=False)
   descriptions=[x['rcsb_polymer_entity']['pdbx_description'].upper() for x in e['polymer_entities']];row['component_inventory']=dict(explicit_class_I=any('CLASS I ' in d for d in descriptions),explicit_class_II=any('CLASS II ' in d for d in descriptions),beta2m=any('MICROGLOBULIN' in d for d in descriptions),Fab=any('FAB' in d for d in descriptions),n_entities=len(descriptions))
   row['fragment_entity_pairs']=[dict(alpha=a['rcsb_id'],beta=b['rcsb_id'],prefix_before_alpha_fragment=a['entity_poly']['pdbx_seq_one_letter_code_can'].split(seqs['alpha'])[0],prefix_before_beta_fragment=b['entity_poly']['pdbx_seq_one_letter_code_can'].split(seqs['beta'])[0]) for a,b in pairs]
   if pdb not in known:
    page=f/'annotation'/(pdb.lower()+'.html');row['external_annotation']=absent_numbering(page);text=page.read_text();domains=numbered_domains(text);checks=[]
    entity_by_chain={chain:entity for entity in e['polymer_entities'] for chain in entity['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']}
    for ac,bc in chain_pairs(text):
     pair=dict(alpha_chain=ac,beta_chain=bc,roles={},complete_external_sequence_gate_pass=False,identity_admitted=False)
     for role,chain in [('alpha',ac),('beta',bc)]:
      try:
       entity=entity_by_chain[chain];ann=domains[chain];assert ann['kind']==('VA' if role=='alpha' else 'VB')
       domain=validate_domain(entity['entity_poly']['pdbx_seq_one_letter_code_can'],ann,(f/'fasta'/(pdb.lower()+'_'+chain+'.fasta')).read_text(),text,refs)
       pair['roles'][role]=dict(sequence_gate_pass=True,variable=domain['variable'],v_allele=domain['v_allele'],variable_start0=domain['variable_start0'],variable_end0=domain['variable_end0'],n_terminal_extension=domain['n_terminal_extension'],numbered_positions=domain['numbered_positions'])
      except (ValueError,KeyError,OSError,AssertionError) as error:pair['roles'][role]=dict(sequence_gate_pass=False,reason=str(error) or type(error).__name__)
     pair['complete_external_sequence_gate_pass']=all(x['sequence_gate_pass'] for x in pair['roles'].values())
     if pair['complete_external_sequence_gate_pass']:pair['paired_variable_sequence_key']=pair_key(pair['roles']['alpha']['variable'],pair['roles']['beta']['variable'])
     checks.append(pair)
    row['external_sequence_checks']=checks
   records.append(row)
 assert set(discovery)==set(metadata)
 checked_sources(BASE/'sources');np_status={p:absent_numbering(BASE/'sources'/(p+'_stcrdab.html')) for p in ['9nu7','9nv6','9nv2','9nve','9nw2','9nva']}
 weekly=json.loads((ROOT/'docs/fingerprint_np105_validation_20261029/sources/weekly_metadata/weekly_entities.json').read_text());tcr1414=[]
 for e in weekly['data']['entries']:
  if e['rcsb_id'] not in ['9XC3','9X9I','9X87','9X85']:continue
  descriptions=[x['rcsb_polymer_entity']['pdbx_description'] for x in e['polymer_entities']];nII=sum('HLA class II' in d for d in descriptions)
  tcr1414.append(dict(pdb=e['rcsb_id'],molecules=descriptions,n_explicit_class_II_entities=nII,scope='class_II_outside_current_class_I' if nII else 'unbound_TCR_outside_current_complete_pMHC_scope',qualifies_current_class_I=False))
 fresh=[]
 for pdb in ['7RTR','4JFE']:
  e=json.loads((BASE/'sources'/('rcsb_'+pdb+'.json')).read_text());fresh.append(dict(pdb=pdb,latest_revision=e['rcsb_accession_info']['revision_date'],changed_since_previous_coordinate_recheck=e['rcsb_accession_info']['revision_date'][:10]>'2026-10-29',coordinate_download_repeated=False))
 new_release=next(x for x in checked_sources(BASE/'sources') if x['file']=='new_release_search.json')
 result=dict(date='2026-11-01',n_fragment_targets=len(plan['targets']),n_queries=len(plan['queries']),n_novel_exact_role_sequences=sum(not x for t in plan['targets'] for x in t['previous_exact_queries'].values()),n_fresh_metadata_PDBs=len(metadata),n_exact_fragment_pair_leads=sum(x['exact_distinct_entity_fragment_pair'] for x in records),n_local_alignment_false_pair_leads=sum(x['local_alignment_only'] for x in records),unreviewed_PDBs=sorted(set(metadata)-known),new_external_sequence_gate_pass_PDBs=sorted({x['pdb'] for x in records if any(c['complete_external_sequence_gate_pass'] for c in x.get('external_sequence_checks',[]))}),records=records,np_external_numbering=np_status,tcr1414_scope=tcr1414,scope_entities_source='docs/fingerprint_np105_validation_20261029/sources/weekly_metadata/weekly_entities.json',scope_entities_source_sha256=sha((ROOT/'docs/fingerprint_np105_validation_20261029/sources/weekly_metadata/weekly_entities.json').read_bytes()),independent_4PRH_numbering=numbering_diagnosis(),coordinate_revision_recheck=fresh,new_release_query_status=new_release['status'],new_release_search_interpretation='No hits for recorded T-cell-receptor full-text plus initial_release>=2026-10-29 query; not a database-wide proof of absence.',new_identity_admissions=0,new_QC_cross_study_pairs=0,new_native_same_source_same_V_groups=0,protocol_frozen=False,training_started=False,independent_acceptance_accessed=False)
 return result
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 if a.out.exists():raise FileExistsError(a.out)
 result=audit();a.out.write_text(json.dumps(result,indent=2)+'\n');print({k:result[k] for k in ['n_fragment_targets','n_queries','n_novel_exact_role_sequences','n_fresh_metadata_PDBs','n_exact_fragment_pair_leads','n_local_alignment_false_pair_leads','unreviewed_PDBs']})
if __name__=='__main__':main()
