"""Verify primary-source facts and exact patent/PDB sequence bridges; no identity admissions."""
from pathlib import Path
import hashlib, html, json, re
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"docs/fingerprint_np105_context_20261030"
OLD=ROOT/"docs/fingerprint_np105_validation_20261029"
def sha(data):return hashlib.sha256(data).hexdigest()
def main():
    source=BASE/"sources/imgt_gene_curl.html"
    text=" ".join(html.unescape(re.sub("<[^>]+>"," ",source.read_text())).split())
    note="Non canonical J-REGION: Phe-Ala-X-Gly instead of Phe-Gly-X-Gly."
    assert note in text and "TRBJ1-6 TRBJ1-6*01 ORF (4)" in text
    ref=json.loads((ROOT/"docs/fingerprint_np105_raw_20261029/raw_sequence_validation.json").read_text())["domains"]["TRB"]
    mapping={str(pos[0])+pos[1].strip():aa for pos,aa in ref["classic_mapping"]}
    motif="".join(mapping[str(i)] for i in range(118,122));assert motif=="FAAG"
    j=ref["J_reference_nt_exact_substrings"];assert len(j)==1 and j[0]["function"]=="ORF" and j[0]["gene"]=="TRBJ1-6*01"
    report=dict(source_URL="https://www.imgt.org/IMGTrepertoire/LocusGenes/genetable/autotable.php?group=TRBJ&species=Mus+musculus",source_SHA256=sha(source.read_bytes()),official_gene="TRBJ1-6*01",functionality="ORF",official_note4=note,observed_NP105_numbered_positions_118_to_121=motif,observed_full_J_nt=j[0]["reference_nt"],full_J_length_nt=53,observed_full_J_exact_in_original_ORF=True,reason="noncanonical J amino-acid motif; expression of a rearranged receptor does not change this germline annotation",not_claimed="ORF does not mean no expression or proven nonfunctionality of rearranged NP105",reference_annotation_unchanged=True,functional_F_gate_unchanged=True,new_identity_admissions=0,protocol_frozen=False,independent_acceptance_accessed=False)
    (BASE/"J_ORF_reason.json").write_text(json.dumps(report,indent=2)+"\n")
    official=BASE/"sources/epo_EP4223782A1.pdf";mirror=BASE/"sources/EP4223782A1.pdf"
    assert official.read_bytes()==mirror.read_bytes()
    assert sha(official.read_bytes())=="abe9aa7708099cc1b122a538d8494b4551e4a096ae35436e385976fc2645d7df"
    rows=json.loads((BASE/"patent_sequence_transcription.json").read_text())
    numbered={x["name"]:x for x in json.loads((OLD/"weekly_numbering.json").read_text())["rows"]}
    entries=json.loads((OLD/"sources/weekly_metadata/weekly_entities.json").read_text())["data"]["entries"]
    entities={x["rcsb_id"]:x for e in entries for x in e["polymer_entities"]}
    bridges=[]
    for row in rows:
        domain=numbered[row["entity"]]["boundary"]["sequence"]
        full=row["reported_full_sequence"]
        deposited=entities[row["entity"]]["entity_poly"]["pdbx_seq_one_letter_code_can"]
        assert full.count(domain)==deposited.count(domain)==1
        assert full[row["domain_offset0"]:]==domain
        assert domain==row["domain_from_prior_PDB_numbering"]
        assert sha(full.encode())==row["reported_sequence_sha256"]
        bridges.append(dict(SEQ_ID=row["SEQ_ID"],PDB_entity=row["entity"],observed_reported_length=len(full),domain_length=len(domain),reported_domain_offset0=row["domain_offset0"],exact_full_numbered_domain_match=True))
    assert len(bridges)==4
    citation=[]
    for pdb in ["9XC3","9X9I","9X87","9X85"]:
        d=json.loads((BASE/f"sources/rcsb_{pdb}.json").read_text());c=d["rcsb_primary_citation"]
        assert c["journal_abbrev"]=="To Be Published"
        assert not c.get("pdbx_database_id_DOI") and not c.get("pdbx_database_id_PubMed")
        citation.append(dict(pdb=pdb,title=c["title"],authors=c["rcsb_authors"],publication_status=c["journal_abbrev"],DOI=None,PubMed=None,polymer_entity_ids=d["rcsb_entry_container_identifiers"]["polymer_entity_ids"]))
    assert len(citation[2]["polymer_entity_ids"])==2 and len(citation[3]["polymer_entity_ids"])==3
    record=dict(patent="EP4223782A1",official_pdf_URL="https://data.epo.org/publication-server/rest/v1.0/publication-dates/20230809/patents/EP4223782NWA1/document.pdf",official_pdf_SHA256=sha(official.read_bytes()),official_and_mirror_byte_exact=True,source_page=18,source_table=1,embedded_images_manually_reviewed=True,full_page_visually_reviewed=True,sequence_bridges=bridges,structure_citations=citation,prior_exact_pair_hits_unchanged={"TCR1414_1":["9X87","9XC3"],"TCR1414_2":["9X9I"]},patent_to_PDB_full_domain_sequence_bridge_established=True,paper_to_PDB_complete_sequence_not_established_by_CDR3_alone=True,native_pair_proof_missing=True,independent_structural_study_pair_missing=True,patent_not_counted_as_independent_structural_study=True,structural_DOI_not_imputed_from_generation_paper=True,new_identity_admissions=0,new_qualified_cross_study_pairs=0,new_qualified_native_control_groups=0,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    (BASE/"patent_source_audit.json").write_text(json.dumps(record,indent=2)+"\n")
    print("J reason verified; all four patent/PDB bridges exact; qualified pairs remain zero")
if __name__=="__main__":main()
