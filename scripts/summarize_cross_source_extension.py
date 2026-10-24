"""Join immutable source identity, atomic QC and same-study control evidence."""
from pathlib import Path
from collections import defaultdict
import json,itertools,hashlib,xml.etree.ElementTree as ET
from meta_acdc.structure.source_balanced_retrieval import source_balanced_blocks
ROOT=Path(__file__).resolve().parents[1]
def qualified_rows(identity_rows,qc_rows):
 by={}
 for r in identity_rows:
  if r['status']!='identity_admitted_structure_QC_pending':continue
  k=(r['pdb'],r['alpha_chain'],r['beta_chain'])
  if k in by:raise ValueError('duplicate exact source-chain key')
  by[k]=r
 out=[]
 for q in qc_rows:
  if q.get('status')!='native_graph_QC_pass':continue
  k=(q['pdb'],q['alpha_chain'],q['beta_chain'])
  if k not in by or by[k]['receptor_id']!=q['receptor_id']:raise ValueError('QC/complete-variable identity mismatch')
  out.append(dict(by[k],native_graph_sha256=q['graph_sha256'],unit_metadata=q['unit_metadata']))
 return out

def control_groups(rows):
 grouped=defaultdict(list)
 for r in rows:
  if r.get('source_study_key'):grouped[(r['source_study_key'],r['V_context'])].append(r)
 out=[]
 for (study,context),members in sorted(grouped.items()):
  rids=sorted({r['receptor_id'] for r in members})
  if len(rids)<2:continue
  pairs=[]
  for a,b in itertools.combinations(members,2):
   if a['receptor_id']==b['receptor_id']:continue
   same=a['unit_metadata']['pmhc_construct_key']==b['unit_metadata']['pmhc_construct_key']
   pairs.append(dict(pdb_a=a['pdb'],pdb_b=b['pdb'],receptor_a=a['receptor_id'],receptor_b=b['receptor_id'],same_complete_pmhc_construct=same,peptide_a=a['unit_metadata']['peptide'],peptide_b=b['unit_metadata']['peptide']))
  out.append(dict(study=study,V_context=context,receptor_ids=rids,pdbs=sorted({r['pdb'] for r in members}),pairs=pairs,independent_clonotypes_proven=False,meaning='different complete engineered variable sequences; not experimental nonbinding'))
 return out

def main():
 base=ROOT/'docs/fingerprint_cross_source_extension_20261024';identity=json.loads((base/'source_balanced_design.json').read_text());qc=json.loads((base/'native_atomic_qc.json').read_text());rows=qualified_rows(identity['rows'],qc['rows']);groups=[]
 for g in identity['cross_study_positive_groups']:
  members=[r for r in rows if r['receptor_id']==g['receptor_id']];pairs=[]
  for a,b in itertools.combinations(members,2):
   if a['pdb']==b['pdb'] or not a['source_study_key'] or not b['source_study_key'] or a['source_study_key']==b['source_study_key']:continue
   pairs.append(dict(pdb_a=a['pdb'],pdb_b=b['pdb'],study_a=a['source_study_key'],study_b=b['source_study_key']))
  groups.append(dict(receptor_id=g['receptor_id'],V_context=g['V_context'],identity_supported_pdbs=g['pdbs'],identity_author_disjoint=g['has_author_disjoint_pair'],QC_pass_pdbs=sorted({r['pdb'] for r in members}),QC_cross_study_pairs=pairs))
 controls=control_groups(rows);blocks=source_balanced_blocks(rows);strict=source_balanced_blocks(rows,True)
 summary=dict(n_identity_admitted_pdbs=identity['summary']['n_identity_admitted_pdbs'],n_complete_paired_ids=identity['summary']['n_complete_paired_ids'],n_identity_cross_study_groups=len(groups),n_QC_graphs=len(rows),n_QC_pdbs=len({r['pdb'] for r in rows}),n_QC_complete_paired_ids=len({r['receptor_id'] for r in rows}),n_QC_cross_study_pairs=sum(len(g['QC_cross_study_pairs']) for g in groups),n_QC_same_study_same_V_control_groups=len(controls),n_QC_source_balanced_blocks=len(blocks),n_QC_author_disjoint_balanced_blocks=len(strict))
 result=dict(summary=summary,cross_study_groups=groups,same_study_same_V_controls=controls,source_balanced_blocks=blocks,author_disjoint_balanced_blocks=strict,status='blocked_before_protocol_freeze',blockers=['no atomic-QC-passing cross-study same-complete-variable pair','no source-balanced query block; same-V controls do not repair missing cross-study positives','fewer than three eligible disjoint source/family components; minimum computability is not statistical power'],protocol_frozen=False,training_release_created=False,model_training_started=False,independent_acceptance_accessed=False,scope='identity retrieval construction evidence; no experimental negative labels')
 (base/'joint_readiness.json').write_text(json.dumps(result,indent=2)+'\n')
 root=ET.parse(base/'sources/jdi_paper/article.xml').getroot();table=next(t for t in root.iter('table-wrap') if t.get('id')=='Tab4');tab=[[''.join(c.itertext()) for c in tr] for tr in table.findall('.//tbody/tr')];proof=[]
 for pdb,name in [('7pb2','JDI (parent TCR)'),('7ow5','JDIa41b1 (high affinity TCR)'),('7ow6','JDIa41b1 (high affinity TCR)')]:
  r=next(r for r in rows if r['pdb']==pdb);t=next(t for t in tab if t[0]==name);matches=[s in r['alpha']['variable'] for s in t[1:4]]+[s in r['beta']['variable'] for s in t[4:6]];assert all(matches)
  proof.append(dict(pdb=pdb,construct_name=name,receptor_id=r['receptor_id'],CDR_substrings_match_primary_table_4=matches,full_variable_sequences_from_primary_SEQRES_and_numbering=True))
 receipt=dict(doi='10.1038/s41467-022-32811-1',pmcid='PMC9464187',primary_xml_sha256=hashlib.sha256((base/'sources/jdi_paper/article.xml').read_bytes()).hexdigest(),proof=proof,finding='7PB2 parent JDI and 7OW5/7OW6 affinity-enhanced JDIa41b1 are one engineered lineage; 7OW6/7PB2 share deposited pMHC construct and G12D peptide; no independent-clonotype claim',source_label='same research article, not cross-study replication',binding_numbers_not_used_for_training=True,independent_acceptance_accessed=False)
 (base/'jdi_construct_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
