"""Exact complete variable-pair curation from frozen primary deposition/API snapshots."""
from pathlib import Path
import json,re,hashlib,argparse
from collections import defaultdict,Counter
from meta_acdc.structure.paired_variable_identity import germline_references,numbered_domains,validate_domain,pair_key
ROOT=Path(__file__).resolve().parents[1]
def study_key(citation):
    doi=citation.get('pdbx_database_id_DOI');pmid=citation.get('pdbx_database_id_PubMed')
    if doi:return 'doi:'+doi.strip().lower()
    if pmid:return 'pmid:'+str(pmid)
    return None

def chain_pairs(text):
    return sorted(set((a.strip(),b.strip()) for b,a in re.findall(r'<b>VB chain</b></td><td[^>]*>([^<]+)</td></tr>\s*<tr><td[^>]*><b>VA chain</b></td><td[^>]*>([^<]+)</td>',text)))

def main():
    base=ROOT/'docs/fingerprint_cross_source_20261024';src=base/'sources';annotations=src/'stcrdab';oldsrc=ROOT/'docs/fingerprint_paired_identity_20261024/sources'
    ap=argparse.ArgumentParser();ap.add_argument('--selection',type=Path,default=base/'candidate_selection_with_quarantine_followup.json');ap.add_argument('--out',type=Path,default=base/'sequence_curation_v2.json');args=ap.parse_args()
    refs=germline_references(list(oldsrc.glob('Homo_sapiens_TR*.fasta'))+list(oldsrc.glob('Mus_musculus_TR*.fasta')))
    entries={e['rcsb_id'].lower():e for path in src.glob('rcsb_batch_*.json') for e in json.loads(path.read_text())['data']['entries'] if e}
    if (src/'rcsb_extra.json').exists():entries.update({e['rcsb_id'].lower():e for e in json.loads((src/'rcsb_extra.json').read_text())['data']['entries'] if e})
    selection=json.loads(args.selection.read_text());rows=[];entry_failures=[]
    for pdb in selection['selected_pdbs']:
        entry=entries[pdb];citation=entry.get('rcsb_primary_citation') or {};textpath=annotations/(pdb+'.html')
        try:
            text=textpath.read_text();domains=numbered_domains(text);pairs=chain_pairs(text)
            if not pairs:raise ValueError('no explicit external alpha/beta pairing')
        except (ValueError,OSError) as error:
            entry_failures.append(dict(pdb=pdb,reason=str(error),study_key=study_key(citation)));continue
        entities={}
        for e in entry.get('polymer_entities') or []:
            for chain in (e.get('rcsb_polymer_entity_container_identifiers') or {}).get('auth_asym_ids',[]):
                if chain in entities:raise ValueError('ambiguous author chain in primary API')
                entities[chain]=e
        for a,b in pairs:
            row=dict(pdb=pdb,alpha_chain=a,beta_chain=b,source_study_key=study_key(citation),primary_citation=citation,status='quarantined',reasons=[],scope='complete deposited engineered variable sequences; not binding truth',structural_QC='pending')
            for role,chain,kind in [('alpha',a,'VA'),('beta',b,'VB')]:
                try:
                    entity=entities[chain];seq=entity['entity_poly']['pdbx_seq_one_letter_code_can']
                    desc=(entity.get('rcsb_polymer_entity') or {}).get('pdbx_description','')
                    row[role]=dict(entity_id=entity['rcsb_id'],full_construct=seq,entity_description=desc)
                    candidates=[chain]+sorted(c for c in domains if c!=chain and c in entities and entities[c]['rcsb_id']==entity['rcsb_id'])
                    results=[];errors=[]
                    for sourcechain in candidates:
                        try:
                            ann=domains[sourcechain]
                            if ann['kind']!=kind:raise ValueError('external role mismatch')
                            annotation=annotations/(pdb+'_'+sourcechain+'.fasta');item=validate_domain(seq,ann,annotation.read_text(),text,refs)
                            results.append(dict(item,domain_candidate=ann,annotation_chain=sourcechain,annotation_sha256=hashlib.sha256(annotation.read_bytes()).hexdigest(),annotation_inheritance='same RCSB polymer entity with identical complete construct' if sourcechain!=chain else 'direct chain'))
                        except (ValueError,KeyError,OSError) as error:errors.append(str(error))
                    if not results:raise ValueError('; '.join(sorted(set(errors))))
                    if len({x['variable'] for x in results})!=1:raise ValueError('same-entity variable annotation conflict')
                    row[role].update(results[0])
                except (ValueError,KeyError,TypeError,OSError) as error:row['reasons'].append(role+': '+str(error))
            if not row['reasons']:
                if row['alpha']['entity_id']==row['beta']['entity_id']:row['reasons'].append('single-chain fusion outside paired-chain scope')
                elif row['alpha']['species']!=row['beta']['species']:row['reasons'].append('species mismatch')
                else:
                    row.update(status='sequence_supported_pending_primary_header',receptor_id=pair_key(row['alpha']['variable'],row['beta']['variable']),construct_id=pair_key(row['alpha']['full_construct'],row['beta']['full_construct']),V_context='|'.join([row['alpha']['species'],row['alpha']['v_allele'].split('*')[0],row['beta']['v_allele'].split('*')[0]]))
            rows.append(row)
    accepted=[x for x in rows if x['status']=='sequence_supported_pending_primary_header'];groups=defaultdict(list)
    for row in accepted:groups[row['receptor_id']].append(row)
    multi=[]
    for rid,items in groups.items():
        pdbs=sorted({x['pdb'] for x in items});studies=sorted({x['source_study_key'] for x in items if x['source_study_key']})
        if len(pdbs)>1:multi.append(dict(receptor_id=rid,pdbs=pdbs,studies=studies,n_studies=len(studies),V_context=items[0]['V_context']))
    summary=dict(n_selected_pdbs=len(selection['selected_pdbs']),n_annotated_paired_units=len(rows),n_sequence_supported_units=len(accepted),n_supported_pdbs=len({x['pdb'] for x in accepted}),n_complete_paired_ids=len(groups),n_multi_pdb_ids=len(multi),n_multi_study_ids=sum(x['n_studies']>1 for x in multi),quarantine_reasons=dict(Counter(z for x in rows for z in x['reasons'])),n_entry_annotation_failures=len(entry_failures))
    out=dict(summary=summary,identity_semantics='same complete paired variable amino-acid sequences, engineering retained; no name-based merge',rows=rows,entry_failures=entry_failures,multi_pdb_groups=multi,independent_acceptance_accessed=False)
    target=args.out
    if target.exists():raise FileExistsError(target)
    target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(summary,indent=2));print('MULTI_STUDY',json.dumps([x for x in multi if x['n_studies']>1],indent=2))
if __name__=='__main__':main()
