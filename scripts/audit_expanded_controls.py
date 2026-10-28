"""Audit expanded discovery evidence without granting identity or atomic admission."""
from pathlib import Path
import argparse,hashlib,json
from scripts.curate_cross_source_receptors import chain_pairs
from meta_acdc.structure.paired_variable_identity import numbered_domains,germline_references
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'docs/fingerprint_expanded_controls_20261029'
SCOPE_REASONS={
'1tce':'SHC plus CD3 zeta phosphopeptide; no paired alpha/beta variable domains',
'1ygr':'CD45 phosphatase plus CD3 zeta ITAM; no paired alpha/beta variable domains',
'2e5k':'Sts-1 signaling suppressor SH3 domain; not a paired TCR',
'2h0q':'Sts-1 signaling suppressor PGM domain; not a paired TCR',
'2ikq':'Sts-1 signaling suppressor PGM domain; not a paired TCR',
'3ri4':'Ets1 plus TCR alpha promoter DNA; not a paired receptor protein',
'6mf8':'TCR alpha transmembrane fragment; no complete paired variable domains',
'6u07':'engineered alpha/beta constant domains only; no paired variable domains',
'7k0z':'marsupial mu/gamma TCR; outside alpha/beta class-I scope',
'9ci8':'gamma/delta-CD3 complex; outside alpha/beta class-I scope',
'9cia':'gamma/delta-CD3 complex; outside alpha/beta class-I scope',
'9jy2':'gamma/delta-CD3-Fab complex; outside alpha/beta class-I scope',
'9jy3':'gamma/delta-CD3-Fab complex; outside alpha/beta class-I scope',
'9zzv':'engineered orthogonal alpha/beta constant domains only; no variable domains'}

def metadata():
    result={}
    for path in list((BASE/'sources/metadata_seed').glob('batch*.json'))+[BASE/'sources/followup/3ri4.json']+list((BASE/'complete_search/sources/metadata').glob('batch*.json')):
        data=json.loads(path.read_text());assert not data.get('errors')
        for e in data['data']['entries']:
            assert e is not None
            result[e['rcsb_id'].lower()]=e
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    if args.out.exists():raise FileExistsError(args.out)
    entries=metadata()
    old={e['rcsb_id'].lower() for p in (ROOT/'docs/fingerprint_cross_source_20261024/sources').glob('rcsb_batch_*.json') for e in json.loads(p.read_text())['data']['entries'] if e}
    inventory=set();query_counts={}
    for p in [BASE/'sources/expanded/rcsb_phrase_0.json',BASE/'sources/expanded/rcsb_phrase_1.json',BASE/'sources/tcr/rcsb_tcr.json']:
        d=json.loads(p.read_text());assert d['total_count']==len(d['result_set']);query_counts[str(p.relative_to(BASE))]=d['total_count'];inventory|={x['identifier'].lower() for x in d['result_set']}
    novel=sorted(inventory-old);assert set(novel)==set(SCOPE_REASONS)
    outside=[dict(pdb=p,title=entries[p]['struct']['title'],primary_citation=entries[p].get('rcsb_primary_citation'),entities=[dict(entity=e['rcsb_id'],description=e['rcsb_polymer_entity']['pdbx_description']) for e in entries[p]['polymer_entities']],reason=SCOPE_REASONS[p],admitted=False) for p in novel]
    facts=json.loads((BASE/'natural_primary_cells.json').read_text())
    refs=germline_references(list((ROOT/'docs/fingerprint_paired_identity_20261024/sources').glob('Homo_sapiens_TR*.fasta'))+list((ROOT/'docs/fingerprint_paired_identity_20261024/sources').glob('Mus_musculus_TR*.fasta')))
    natural=[]
    for item in facts['rows']:
        p=item['pdb'];entry=entries[p];checks={}
        for role,kind,prefix,jprefix in [('alpha','VA','TRAV','TRAJ'),('beta','VB','TRBV','TRBJ')]:
            seq=item['chains'][role]['primary_labeled_variable'];cdr3=item['CDR3_'+role]
            matches=[e for e in entry['polymer_entities'] if e['entity_poly']['pdbx_seq_one_letter_code_can'].count(seq)==1 and seq.count(cdr3)==1]
            assert len(matches)==1,(p,role)
            checks[role]=dict(primary_cell=item['chains'][role],entity=matches[0]['rcsb_id'],author_chains=matches[0]['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'],exact_primary_variable_in_deposition=True,CDR3=cdr3,CDR3_unique=True,
                functional_mouse_V_prefix_matches=sorted(g for (sp,g),rr in refs.items() if sp=='mouse' and g.startswith(prefix) and any(seq.startswith(r['amino_acids'][:10]) for r in rr)),
                functional_mouse_J_suffix_matches=sorted(g for (sp,g),rr in refs.items() if sp=='mouse' and g.startswith(jprefix) and any(seq.endswith(r['amino_acids'][-10:]) for r in rr)))
        assert checks['alpha']['entity']!=checks['beta']['entity']
        text=(BASE/'sources/annotation'/(p+'.html')).read_text()
        pairs=chain_pairs(text);domains=numbered_domains(text)
        assert not pairs and not domains
        natural.append(dict(name=item['name'],pdb=p,primary_TRAV_TRAJ=item['TRAV_TRAJ'],primary_TRBV_TRBJ=item['TRBV_TRBJ'],constant_domain=item['constant_domain'],chain_checks=checks,external_numbered_pairs=pairs,numbered_domain_count=len(domains),fixed_identity_admitted=False,atomic_QC_run=False,reason='primary sequence/role evidence restored, but no valid external numbered full-domain pairing; fixed identity gate remains unresolved'))
    quartet=[r for r in facts['rows'] if r['name'] in ['NP202','NP101','NP207','NP203']]
    assert len({r['TRAV_TRAJ'] for r in quartet})==len({r['TRBV_TRBJ'] for r in quartet})==1
    assert len({(r['CDR3_alpha'],r['CDR3_beta']) for r in quartet})==4
    cross=json.loads((BASE/'complete_search/primary_validation_both_ends.json').read_text())
    natural_search=json.loads((BASE/'natural_search/search_results.json').read_text())
    output=dict(description_search_query_counts=query_counts,n_description_union=len(inventory),n_previous_archive=len(old),n_outside_previous_archive=len(novel),new_description_leads=outside,
        n_complete_paired_targets=44,n_complete_chain_queries=87,n_complete_queries_successful=87,n_complete_search_pdbs=77,n_full_cross_study_discovery_leads=len(cross['complete_cross_study_leads']),
        n_new_qualified_cross_study_pairs=cross['n_new_qualified_cross_study_leads'],
        natural_source=dict(primary_preprint_doi=facts['article_doi'],publication_doi=facts['primary_publication_doi'],source_status=facts['article_status'],workbook_sha256=facts['workbook_sha256'],evidence='distinct naturally observed paired clonotypes in a pooled mouse repertoire; not independent donor assignments or independent studies',same_V_label_quartet=[r['name'] for r in quartet],V_label_caveat='primary TRAV16 label and TCR3d TRAV16D/DV11 nomenclature are not normalized into an allele identity without numbered/germline evidence',rows=natural),
        n_natural_discovery_queries=natural_search['n_queries'],n_natural_successful_queries=natural_search['n_successful_queries'],natural_search_pdbs=natural_search['discovery_pdbs'],
        current_qualified_inventory=dict(identity_PDBs=64,paired_variable_IDs=44,QC_graphs=43,QC_PDBs=31,QC_IDs=29,cross_study_positive_pairs=0,independent_natural_control_groups=0,source_balanced_blocks=0),
        protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False,
        limits='Finite description and exact-sequence discovery is not exhaustive public-database coverage. HTTP failures are unresolved access, never negative evidence. No atomic QC is claimed for identity-blocked candidates.')
    args.out.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k.startswith('n_') or k.endswith('pdbs')},indent=2))
if __name__=='__main__':main()
