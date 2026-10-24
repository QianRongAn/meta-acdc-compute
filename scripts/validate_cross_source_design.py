"""Validate raw deposit headers and source-balanced design before any model protocol."""
from pathlib import Path
import json,hashlib,itertools,re
from collections import Counter,defaultdict
from meta_acdc.structure.receptor_identity_audit import seqres_sequences
from meta_acdc.structure.identity import pdb_entities
from meta_acdc.structure.source_balanced_retrieval import author_keys,source_balanced_blocks,joint_components
ROOT=Path(__file__).resolve().parents[1]
def main():
    base=ROOT/'docs/fingerprint_cross_source_20261024';data=json.loads((base/'sequence_curation_v2.json').read_text());rows=[]
    for original in data['rows']:
        row=dict(original)
        if row['status']=='sequence_supported_pending_primary_header':
            path=base/'sources/headers'/(row['pdb']+'.pdb')
            try:
                text=path.read_text();seq=seqres_sequences(text);entities=pdb_entities(path)
                for role,kind in [('alpha','tra'),('beta','trb')]:
                    chain=row[role+'_chain']
                    declared=entities.get(chain,{})
                    if declared.get('kind')!=kind:
                        description=declared.get('description','').upper()
                        own='ALPHA' if role=='alpha' else 'BETA';opposite='BETA' if role=='alpha' else 'ALPHA'
                        if declared.get('kind') not in ['unknown',None] or not re.search(r'\b'+own+r'\b',description) or re.search(r'\b'+opposite+r'\b',description):raise ValueError(role+': primary COMPND role conflict or unspecified')
                        row[role]['primary_role_basis']='explicit primary ALPHA/BETA designation plus external TRAV/TRBV domain; legacy parser lacks TCR-name token'
                    if seq.get(chain)!=row[role]['full_construct']:raise ValueError(role+': primary complete SEQRES/API mismatch')
                row.update(status='identity_admitted_structure_QC_pending',primary_header_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),primary_header=str(path.relative_to(ROOT)))
            except (OSError,ValueError) as error:row.update(status='quarantined',reasons=row['reasons']+['primary header: '+str(error)])
        rows.append(row)
    admitted=[x for x in rows if x['status']=='identity_admitted_structure_QC_pending'];groups=defaultdict(list)
    for x in admitted:groups[x['receptor_id']].append(x)
    multi=[]
    for rid,items in groups.items():
        bypdb={}
        for x in sorted(items,key=lambda x:(x['pdb'],x['alpha_chain'],x['beta_chain'])):bypdb.setdefault(x['pdb'],x)
        distinct=[]
        for a,b in itertools.combinations(bypdb.values(),2):
            if not a['source_study_key'] or not b['source_study_key'] or a['source_study_key']==b['source_study_key']:continue
            ak,bk=author_keys(a),author_keys(b);shared=ak&bk
            distinct.append(dict(pdb_a=a['pdb'],pdb_b=b['pdb'],study_a=a['source_study_key'],study_b=b['source_study_key'],authors_known=bool(ak and bk),shared_author_keys=sorted(shared),author_disjoint=bool(ak and bk and not shared)))
        if distinct:multi.append(dict(receptor_id=rid,V_context=items[0]['V_context'],pdbs=sorted(bypdb),distinct_study_pairs=distinct,has_author_disjoint_pair=any(x['author_disjoint'] for x in distinct)))
    blocks=source_balanced_blocks(admitted);strict=source_balanced_blocks(admitted,True);components=joint_components(admitted)
    eligible_contexts={b['V_context'] for b in blocks};eligible_components=[g for g in components if any(x['V_context'] in eligible_contexts for x in g)]
    controls=[]
    for context in sorted({x['V_context'] for x in admitted}):
        members=[x for x in admitted if x['V_context']==context];rids=sorted({x['receptor_id'] for x in members})
        if len(rids)>1:controls.append(dict(V_context=context,receptor_ids=rids,pdbs=sorted({x['pdb'] for x in members}),studies=sorted({x['source_study_key'] for x in members if x['source_study_key']}),has_cross_study_positive=any(x['V_context']==context for x in multi)))
    result=dict(summary=dict(n_identity_admitted_units=len(admitted),n_identity_admitted_pdbs=len({x['pdb'] for x in admitted}),n_complete_paired_ids=len(groups),n_cross_study_positive_ids=len(multi),n_author_disjoint_positive_ids=sum(x['has_author_disjoint_pair'] for x in multi),n_same_V_contexts_with_distinct_ids=len(controls),n_source_balanced_query_blocks=len(blocks),n_author_disjoint_balanced_blocks=len(strict),n_joint_components=len(components),n_source_balanced_joint_components=len(eligible_components)),rows=rows,entry_annotation_failures=data['entry_failures'],cross_study_positive_groups=multi,same_V_distinct_receiver_controls=controls,source_balanced_query_blocks=blocks,author_disjoint_source_balanced_blocks=strict,joint_components=[dict(pdbs=sorted({x['pdb'] for x in g}),studies=sorted({x['source_study_key'] for x in g if x['source_study_key']}),V_contexts=sorted({x['V_context'] for x in g})) for g in components],training_release_created=False,model_training_started=False,independent_acceptance_accessed=False,minimum_design_gate=dict(source_balanced_query_blocks_required=True,minimum_disjoint_source_family_components=3,rationale='at least one fit and two development-validation components; computational minimum only, not statistical power or biological release',all_training_graphs_require_atomic_QC=True),status='blocked_before_protocol_freeze',blockers=[])
    if not blocks:result['blockers'].append('no query has both exact complete same-receptor positive and same-V distinct-receptor negative from one common different-source gallery study')
    if not strict:result['blockers'].append('no author-disjoint source-balanced query block')
    if len(eligible_components)<3:result['blockers'].append('fewer than three eligible disjoint joint-source/family components')
    result['blockers'].append('structural atom/chain QC pending outside the explicitly checked cross-study candidates; no training release')
    (base/'source_balanced_design.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['summary'],indent=2));print('BLOCKERS',result['blockers'])
if __name__=='__main__':main()
