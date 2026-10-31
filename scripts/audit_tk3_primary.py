"""Verify TK3 source conflict and germline annotation facts, not binding labels."""
from pathlib import Path
import argparse,hashlib,html,json,re,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_tk3_j_evidence_20261101'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def primary_audit():
    source=BASE/'sources/paper2010.xml';root=ET.fromstring(source.read_bytes());assert root.tag=='article'
    doi=next(x.text for x in root.findall('.//article-id') if x.attrib.get('pub-id-type')=='doi');assert doi=='10.1084/jem.20100603'
    paragraphs=[' '.join(''.join(p.itertext()).split()) for p in root.iter('p')]
    assert any('TRAV20/TRAJ58 combined with TRBV9/TRBJ2-2' in t for t in paragraphs)
    assert any('TK3, MW1, and CS1' in t for t in paragraphs)
    assert any('three unrelated' in t for t in paragraphs)
    assert any('cDNA was reverse transcribed from TK3 T cell RNA' in t for t in paragraphs)
    assert any('engineered disulfide linkage in the constant domains' in t for t in paragraphs)
    assert any('All SPR experiments were conducted at 25' in t for t in paragraphs)
    page=BASE/'sources/imgt_gene.html';raw=page.read_text();m=re.search(r"<li id='notes-11'>(.*?)</li>",raw,re.S);assert m
    note=html.unescape(re.sub('<[^>]+>',' ',m.group(1)));assert 'first codon' in note and 'STOP-CODON' in note and 'may disappear during rearrangements' in note
    flat=' '.join(html.unescape(re.sub('<[^>]+>',' ',raw)).split());assert 'TRAJ58 TRAJ58*01 ORF (11)' in flat
    iedb=BASE/'sources/iedb_TK3.html';page_text=iedb.read_text();m=re.search(r'var\s+assayReceptorData\s*=\s*',page_text);data,_=json.JSONDecoder().raw_decode(page_text[m.end():]);receptors=[r for r in data['data'] if r['ref_name']=='TK3'];assert len(receptors)==1;r=receptors[0];assert r['curated_chain1']['j_gene']=='TCRAJ57*01' and r['calculated_chain1']['j_gene']=='TRAJ57*01'
    variable=r['calculated_chain1']['v_dom_seq'];old=json.loads((ROOT/'docs/fingerprint_control_closure_20261101/audit_final_v3.json').read_text());case=next(x for x in old['independent_4PRH_numbering']['rows'] if x['role']=='alpha');assert case['observed_variable'].endswith(variable)
    clone_manifest=json.loads((BASE/'natural_clone_search_retry/manifest.json').read_text());assert next(x for x in clone_manifest if x['file']=='MW1_search.json')['status']=='204';hits=[]
    for item in clone_manifest:
        assert item['returncode']==0 and item['status'] in ['200','204']
        path=BASE/'natural_clone_search_retry'/item['file']
        if item['status']=='200':hits.extend(x['identifier'] for x in json.loads(path.read_text())['result_set'])
    metadata=json.loads((BASE/'clone_metadata/CS1_candidates.json').read_text())['data']['entries'];assert set(hits)=={x['rcsb_id'] for x in metadata}
    exclusions=[]
    for e in metadata:
        descriptions=[p['rcsb_polymer_entity']['pdbx_description'] for p in e['polymer_entities']]
        assert not any('T cell receptor' in x or 'T-cell receptor' in x for x in descriptions)
        exclusions.append(dict(pdb=e['rcsb_id'],title=e['struct']['title'],molecules=descriptions,citation=e['rcsb_primary_citation'],decision='name_homonym_without_paired_TCR_entities',identity_admitted=False))
    return dict(primary_paper=dict(DOI=doi,PMCID='PMC2901058',source_sha256=sha(source),locations=['RESULTS / Allelic variation in the TRBV9 gene','Figure 1a caption','METHODS / TCR cloning, protein expression, and purification'],reported_V_J=dict(alpha=['TRAV20','TRAJ58'],beta=['TRBV9','TRBJ2-2']),natural_clone_leads=['TK3','MW1','CS1'],three_unrelated_individuals_reported=True,TK3_cDNA_from_T_cell_RNA_reported=True,constant_domain_disulfide_engineering_reported=True,complete_engineered_DNA_not_recovered=True),IMGT=dict(source_URL='https://www.imgt.org/IMGTrepertoire/LocusGenes/genetable/autotable.php?group=TRAJ&species=Homo+sapiens',source_sha256=sha(page),gene='TRAJ58*01',functionality='ORF',note_locator='notes-11',meaning='The germline first-codon stop can disappear on rearrangement; germline ORF annotation is not a mature-receptor nonfunctionality verdict.'),IEDB_conflict=dict(source_URL='https://www.iedb.org/details_v3.php?id=2107002&type=assay',source_sha256=sha(iedb),receptor_id=r['curated_receptor_id'],receptor_name=r['ref_name'],curated_alpha_J=r['curated_chain1']['j_gene'],calculated_alpha_J=r['calculated_chain1']['j_gene'],observed_alpha_variable=variable,primary_paper_alpha_J='TRAJ58',sequence_suffix_supports_TRAJ58=True,annotation_status='source_conflict_quarantined',external_record_edited=False,assay_page_development_exposed=True,source_not_added_to_independent_acceptance=True,assay_endpoint_not_converted_to_binding_negative=True),natural_clone_search=dict(MW1_named_query_no_hits=True,CS1_name_hits=hits,excluded_homonyms=exclusions,scope_limit='Named clone-description plus TCR full-text queries only; not a sequence-wide exclusion of unnamed structures.'),full_primary_pages_republished=False,new_identity_admissions=0,new_QC_cross_study_pairs=0,new_natural_control_groups=0,protocol_frozen=False,training_started=False,independent_acceptance_accessed=False)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    result=primary_audit();a.out.write_text(json.dumps(result,indent=2)+'\n');print('Primary TRAJ58 verified; IEDB J57 conflict quarantined; clone-name hits are not TCRs')
if __name__=='__main__':main()
