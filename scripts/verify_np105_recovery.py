"""Verify packaged whole-run NP105 evidence without opening independent acceptance data.
The lightweight verification rechecks archived evidence and tests; it does not rerun
all 50 million source spots or all 92 assemblies. Full replay commands are in README.
"""
from pathlib import Path
import argparse, ast, collections, gzip, hashlib, io, json, os, subprocess, sys, tarfile
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "docs/fingerprint_np105_recovery_20261029"
def sha(data):
    return hashlib.sha256(data).hexdigest()
def load(name):
    return json.loads((BASE / name).read_text())
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    reports = {}; archive_files = {}
    for row in load("archive_layout.json"):
        data = gzip.decompress((BASE / row["packaged"]).read_bytes())
        if row["kind"] == "gzip_exact":
            assert len(data) == row["bytes"] and sha(data) == row["sha256"]
            reports[row["source"]] = json.loads(data)
        else:
            with tarfile.open(fileobj=io.BytesIO(data)) as archive:
                assert {m.name for m in archive.getmembers()} == set(row["files"])
                for member in archive.getmembers():
                    assert member.isfile()
                    body = archive.extractfile(member).read()
                    receipt = row["files"][member.name]
                    assert len(body) == receipt["bytes"] and sha(body) == receipt["sha256"]
                    archive_files[member.name] = body
    selection = load("whole_selection_receipt.json")
    assert selection["complete"] and selection["source_whole_MD5_verified"]
    assert selection["n_source_pairs"] == 50653512
    assert selection["n_selected_pairs"] == 1849447
    assert selection["n_selected_barcodes"] == 92
    assert len(selection["outputs"]) == 186
    assert selection["prefix_all186_selected_outputs_byte_exact"]
    assert not selection["biological_reference_used"] and selection["bases_substituted"] == 0
    raw = load("post_storage_recovery_sra_check.json")
    assert raw == dict(bytes=6040866339, md5="bb8d681b918233fada8e0dbe875e4887", sha256="b7dc3b09865eb381d102f7be3c7ff3e129a0297897e42cfc98e8ca4935628bbb")
    assemblies = reports["whole_assemblies.json"]
    assert len(assemblies["cells"]) == 92
    assert all(x["exit_code"] == 0 for x in assemblies["cells"])
    assert len([n for n, v in archive_files.items() if n.endswith("/transcripts.fasta") and v]) == 88
    for cell in assemblies["cells"]:
        annotation = cell.get("annotation")
        if annotation and annotation.get("sha256"):
            assert sha(archive_files[cell["barcode"] + "/transcripts.fasta"]) == annotation["sha256"]
    classic = reports["whole_classic_numbering.json"]
    assert classic["n_cells"] == 92 and classic["n_CDR3_bearing_ORFs"] == 443
    assert classic["n_complete_ORFs"] == 394
    counts = collections.Counter(c["state"] for c in classic["cells"])
    assert counts == dict(both_complete=60, beta_only_complete=17, alpha_only_complete=7, neither_complete=8)
    assert len(classic["variants"]) == 4
    bases = reports["whole_base_support.json"]
    assert len(bases["cells"]) == 60
    covered = [c for c in bases["cells"] if c["all_reference_bases_have_same_cell_Q30_support"]]
    canonical = [c for c in bases["cells"] if c["observed_canonical_pair_supported"]]
    assert len(covered) == 57 and len(canonical) == 54
    for cell in bases["cells"]:
        derived = all(all(n > 0 for n in role["per_base_unambiguous_UMI_counts"]) for role in cell["roles"].values())
        assert derived == cell["all_reference_bases_have_same_cell_Q30_support"]
        assert cell["observed_canonical_pair_supported"] == (derived and cell["single_full_sequence_per_role_and_exact_prior_nt"])
    neural = reports["whole_neural_numbering.json"]
    disagreements = [x for x in neural["rows"] if x["classic_complete"] and not (x["non_gap_position_map_identical"] and x["neural_role_matches_source"])]
    assert len(neural["rows"]) == 443 and len(disagreements) == 55
    bad_cells = {x["barcode"] for x in disagreements}
    joint = {x["barcode"] for x in canonical} - bad_cells
    summary = load("whole_numbering_disagreement_summary.json")
    assert len(joint) == 26 and joint == set(summary["canonical_Q30_pair_barcodes_all_methods_agree"])
    assert len(load("cell_summary.json")) == 92
    variants = load("whole_variant_raw_support.json")
    assert len(variants["targets"]) == 4 and variants["bases_substituted"] == 0
    assert all(t["support"]["candidate"]["Q30_exact_CDR3_linked_unambiguous_UMIs"] == 0 for t in variants["targets"])
    assert load("variant_refactor_recheck.json")["report_identical"]
    replay = load("whole_predeclared3_replay.json")
    assert len(replay["replays"]) == 3 and replay["all_byte_exact"]
    for row in replay["replays"]:
        assert row["FASTA_byte_exact"] and row["returncode"] == 0
        assert sha(archive_files[row["barcode"] + "/transcripts.fasta"]) == row["sha256"]
    for prov in load("observed_variable_provenance.json"):
        cell = next(c for c in classic["cells"] if c["barcode"] == prov["barcode"])
        observed = next(v for v in cell[prov["role"]] if v["name"] == prov["ORF_name"])
        assert prov["barcode"] in joint and sha(observed["nucleotides"].encode()) == prov["sha256_nt"]
        assert observed["nucleotides"] in (BASE / "observed_variable.nt.fasta").read_text()
        assert observed["nucleotides"] in observed["ORF_nt"]
    source = load("tcr1414_source_audit.json")
    assert source["all_four_chain_queries_successful"]
    assert source["receptor_generation_DOI_not_imputed_as_structure_DOI"]
    assert source["n_new_qualified_cross_study_pairs"] == source["n_new_verified_same_source_same_V_native_groups"] == 0
    for receptor in source["receptors"]:
        assert not receptor["native_single_cell_pairing_established_by_reported_method"]
        assert not receptor["full_variable_sequence_identity_across_papers_established"]
    frozen = load("frozen_integrity_after_storage_recovery.json")
    total = 0
    for delivery in frozen["previous_deliveries"]:
        folder = ROOT / "docs" / delivery["folder"]
        previous = json.loads((folder / "DELIVERY_MANIFEST.json").read_text())
        for name, receipt in previous["files"].items():
            data = (folder / name).read_bytes()
            assert len(data) == receipt["bytes"] and sha(data) == receipt["sha256"], str(folder / name)
            total += 1
    assert total == 4625 and len(frozen["previous_deliveries"]) == 14
    manifest_path = BASE / "DELIVERY_MANIFEST.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        for name, receipt in manifest["files"].items():
            data = (BASE / name).read_bytes()
            assert sha(data) == receipt["sha256"] and len(data) == receipt["bytes"], name
        for name, receipt in manifest["code_files"].items():
            assert sha((ROOT / name).read_bytes()) == receipt["sha256"], name
    parsed = ast.parse((ROOT / "scripts/verify_np105_all_cells.py").read_text())
    modules = next(ast.literal_eval(n.value) for n in ast.walk(parsed) if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "modules" for t in n.targets))
    run = subprocess.run([sys.executable, "-m", "unittest", *modules, "tests.test_np105_recovery"], cwd=ROOT, env=dict(os.environ, PYTHONPATH="src:."), capture_output=True, text=True)
    assert run.returncode == 0 and "Ran 99 tests" in run.stderr, run.stdout + run.stderr
    result = dict(archive_payloads_verified=True, frozen_deliveries=14, unchanged_frozen_files=total, tests_passed=99, source_spots=50653512, selected_pairs=1849447, whitelist_cells=92, transcript_outputs=88, classic_paired_cells=60, canonical_Q30_pairs=54, canonical_Q30_pairs_all_complete_ORFs_numbering_agrees=26, retained_complete_ORF_numbering_disagreements=55, predeclared_cell_replays_byte_exact=3, full92_replayed=False, new_identity_admissions=0, qualified_cross_study_pairs=0, qualified_native_control_groups=0, protocol_frozen=False, model_training_started=False, independent_acceptance_accessed=False)
    if args.out:
        assert args.out.resolve() != manifest_path.resolve()
        args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
if __name__ == "__main__":
    main()
