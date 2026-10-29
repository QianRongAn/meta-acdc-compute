"""Verify context/terminal/source diagnostics without changing the original acceptance gates."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"docs/fingerprint_np105_context_20261030"
OLD=ROOT/"docs/fingerprint_np105_recovery_20261029"
def sha(data):return hashlib.sha256(data).hexdigest()
def read(name):return json.loads((BASE/name).read_text())
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--primary-originals",action="store_true");ap.add_argument("--out",type=Path);a=ap.parse_args()
    old=subprocess.run([sys.executable,"scripts/verify_np105_recovery.py"],cwd=ROOT,capture_output=True,text=True,check=True)
    assert json.loads(old.stdout)["tests_passed"]==99
    previous=json.loads((OLD/"frozen_integrity_after_storage_recovery.json").read_text())["previous_deliveries"]
    previous=previous+[dict(folder=OLD.name)]
    total=0
    for delivery in previous:
        folder=ROOT/"docs"/delivery["folder"]
        manifest=json.loads((folder/"DELIVERY_MANIFEST.json").read_text())
        for name,receipt in manifest["files"].items():
            data=(folder/name).read_bytes();assert len(data)==receipt["bytes"] and sha(data)==receipt["sha256"],str(folder/name);total+=1
    assert len(previous)==15 and total==4666
    code_paths=set()
    for folder in [ROOT/"docs/fingerprint_np105_raw_20261029",ROOT/"docs/fingerprint_np105_validation_20261029",OLD]:
        m=json.loads((folder/"DELIVERY_MANIFEST.json").read_text())
        for name,receipt in m.get("code",m.get("code_files",{})).items():
            data=(ROOT/name).read_bytes();assert sha(data)==receipt["sha256"] and len(data)==receipt["bytes"];code_paths.add(name)
    assert len(code_paths)==28
    path=BASE/"context_audit.json";before=path.read_bytes()
    subprocess.run([sys.executable,"scripts/audit_np105_context.py"],cwd=ROOT,capture_output=True,check=True)
    assert path.read_bytes()==before
    context=read("context_audit.json");assert context["n_complete_ORF_occurrences"]==394 and context["n_disagreements_retained"]==55
    assert context["groups_with_context_dependent_outcomes"]==2 and context["cropped_inputs"]==0
    terminal=read("terminal_audit.json");rows=terminal["cells"];assert len(rows)==60
    canonical={c["barcode"] for c in rows if c["original_canonical_Q30_pair"]}
    supported={c["barcode"] for c in rows if c["both_terminal_targets_supported"]}
    assert len(canonical)==54 and len(supported)==37 and len(canonical&supported)==35
    for cell in rows:
        for role,value in cell["roles"].items():
            assembled=bool(value["assembly_observations"]) and all(x["exact_observed_target"] for x in value["assembly_observations"])
            assert assembled==value["independently_assembled_all_complete_ORFs_have_exact_target"]
            assert value["terminal_supported"]==(assembled and value["also_exact_CDR3_linked_UMIs"]>0)
        assert cell["both_terminal_targets_supported"]==all(v["terminal_supported"] for v in cell["roles"].values())
    neural=json.loads((OLD/"whole_numbering_disagreement_summary.json").read_text())
    nn=set(neural["canonical_Q30_pair_barcodes_all_methods_agree"])
    joint=read("joint_diagnostic_summary.json");assert len(nn&canonical&supported)==17 and set(joint["joint_barcodes"])==nn&canonical&supported
    reference=json.loads((ROOT/"docs/fingerprint_np105_raw_20261029/raw_sequence_validation.json").read_text())["domains"]["TRB"]
    j=read("J_ORF_reason.json");assert j["functionality"]=="ORF" and j["observed_NP105_numbered_positions_118_to_121"]=="FAAG"
    assert reference["J_reference_nt_exact_substrings"][0]["function"]=="ORF"
    patent=read("patent_source_audit.json");assert patent["official_and_mirror_byte_exact"] and len(patent["sequence_bridges"])==4
    transcription=read("patent_sequence_transcription.json")
    for row in transcription:
        seq=row["reported_full_sequence"];domain=row["domain_from_prior_PDB_numbering"]
        assert seq.count(domain)==1 and seq[row["domain_offset0"]:]==domain
    for row in read("sources/rcsb_manifest.json"):
        data=(BASE/"sources"/row["file"]).read_bytes();assert row["curl_returncode"]==0 and sha(data)==row["sha256"]
    assert all(c["publication_status"]=="To Be Published" and c["DOI"] is None for c in patent["structure_citations"])
    if a.primary_originals:
        old_sources={n:(BASE/n).read_bytes() for n in ["patent_source_audit.json","J_ORF_reason.json"]}
        subprocess.run([sys.executable,"scripts/audit_np105_primary_context.py"],cwd=ROOT,capture_output=True,check=True)
        assert all((BASE/n).read_bytes()==data for n,data in old_sources.items())
    test=subprocess.run([sys.executable,"-m","unittest","tests.test_np105_context","-v"],cwd=ROOT,env=dict(os.environ,PYTHONPATH="src:."),capture_output=True,text=True)
    assert test.returncode==0 and "Ran 5 tests" in test.stderr,test.stderr
    if (BASE/"DELIVERY_MANIFEST.json").exists():
        m=read("DELIVERY_MANIFEST.json")
        for name,receipt in m["files"].items():
            data=(BASE/name).read_bytes();assert sha(data)==receipt["sha256"] and len(data)==receipt["bytes"]
        for name,receipt in m["code_files"].items():assert sha((ROOT/name).read_bytes())==receipt["sha256"]
    result=dict(status="pass",tests_passed=104,unchanged_deliveries=15,unchanged_files=total,unchanged_prior_raw_and_validation_code_files=len(code_paths),complete_ORF_context_occurrences=394,disagreements_retained=55,cropped_inputs=0,original_Q30_pairs=54,terminal_J_and_CDR3_supported_original_pairs=35,terminal_and_original_joint_numbering_pairs=17,patent_PDB_complete_domain_matches=4,primary_originals_rechecked=a.primary_originals,whole_SAM_not_replayed_by_lightweight_verifier=True,new_identity_admissions=0,new_qualified_cross_study_pairs=0,new_qualified_native_control_groups=0,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    if a.out:a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
