"""Audit fresh structural leads using the existing variable-domain and atomic gates."""
import argparse
from dataclasses import replace
import json
import re

from meta_acdc.structure.graph import build_interface_graph_from_residues, parse_pdb
from meta_acdc.structure.identity import annotate, class_i_units
from meta_acdc.structure.paired_variable_identity import germline_references, numbered_domains, validate_domain
from meta_acdc.structure.receptor_identity_audit import seqres_sequences
from meta_acdc.structure.rigid_controls import atomic_screen, policy_sha256
from meta_acdc.structure.source_balanced_retrieval import author_keys
from scripts.curate_cross_source_receptors import chain_pairs, study_key
from scripts.join_native_search import pair_candidates, parse_hits
from scripts.join_terminal_rescue import BASE, ROOT, load_receipt


def verified_files(folder):
    for receipt in json.loads((folder / "manifest.json").read_text()):
        path = folder / receipt["file"]
        yield path, receipt, load_receipt(path, receipt)


def entry_map(folder):
    entries = {}
    for path, receipt, body in verified_files(folder):
        if path.suffix != ".json":
            continue
        assert receipt["status"] == "200", path
        data = json.loads(body)
        assert not data.get("errors"), path
        for entry in data["data"]["entries"]:
            assert entry is not None, path
            entries[entry["rcsb_id"].lower()] = entry
    return entries


def deposited_chains(entry):
    return {chain: entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
            for entity in entry["polymer_entities"]
            for chain in entity["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]}


def atomic_audit(path):
    residues = parse_pdb(path)
    annotate(residues, path)
    units, excluded = class_i_units(residues)
    rows = []
    for unit, meta in units:
        try:
            graph = build_interface_graph_from_residues([replace(r, plddt=None) for r in unit])
            graph_error = None
        except ValueError as error:
            graph = None
            graph_error = str(error)
        screen = atomic_screen(unit)
        rows.append(dict(chains=meta, graph_nodes=graph.n_nodes if graph else None,
                         graph_error=graph_error, screen=screen, diagnostic_pass=bool(graph and screen["passed"])))
    return dict(units=rows, scope_exclusions=excluded)


def audit():
    source_design = ROOT / "docs/fingerprint_positive_controls_20261027/source_balanced_design.json"
    old = json.loads(source_design.read_text())
    by_pdb = {}
    for row in old["rows"]:
        by_pdb.setdefault(row["pdb"], []).append(row)
    ready = ROOT / "docs/fingerprint_construct_followup_20261027/joint_readiness_final.json"
    inventory = json.loads(ready.read_text())["summary"]
    metadata = entry_map(BASE / "sources/metadata")
    metadata.update(entry_map(BASE / "sources/control"))
    metadata.update(entry_map(BASE / "sources/recent"))
    plan = json.loads((BASE / "search_plan.json").read_text())
    discovery = json.loads((BASE / "search_results.json").read_text())
    fresh = {}
    reviews = json.loads((BASE / "candidate_decisions.json").read_text())
    priors = {}
    for target, row in zip(plan["targets"], discovery["rows"], strict=True):
        assert target["target_id"] == row["target_id"]
        for pdb, matches in row["exact_entity_matches"].items():
            entry = metadata[pdb]
            citation = entry["rcsb_primary_citation"] or {}
            if pdb in by_pdb:
                priors[pdb] = sorted({r["status"] for r in by_pdb[pdb]})
                continue
            dest = fresh.setdefault(pdb, dict(title=entry["struct"]["title"], citation=citation,
                                             decision=reviews[pdb], matches=[]))
            dest["matches"].append(dict(target_id=target["target_id"], entities=matches,
                                        anchors=sorted({r["pdb"] for r in target["evidence"]})))
    assert set(reviews) == set(fresh)

    ref_dir = ROOT / "docs/fingerprint_paired_identity_20261024/sources"
    refs = germline_references(sorted(ref_dir.glob("*_TR*.fasta")))
    coordinate_checks = {}
    fasta_dir = BASE / "sources/priority_fasta"
    for path, receipt, _ in verified_files(fasta_dir):
        assert receipt["status"] == "200", path
    for folder in ("priority", "recent"):
        for path, receipt, payload in verified_files(BASE / "sources" / folder):
            assert receipt["status"] == "200", path
            if path.suffix != ".pdb":
                continue
            pdb = path.stem
            seqres = seqres_sequences(payload.decode())
            declared = deposited_chains(metadata[pdb])
            text = path.with_suffix(".html").read_text()
            pairs = chain_pairs(text)
            domains = numbered_domains(text)
            validations = []
            for alpha, beta in pairs:
                row = {"alpha": alpha, "beta": beta, "roles": {}}
                for role, chain, kind in (("alpha", alpha, "VA"), ("beta", beta, "VB")):
                    try:
                        assert domains[chain]["kind"] == kind
                        assert declared[chain] == seqres[chain]
                        info = validate_domain(seqres[chain], domains[chain],
                                               (fasta_dir / f"{pdb}_{chain}.fasta").read_text(), text, refs)
                        row["roles"][role] = {"variable": info["variable"], "V": info["v_allele"]}
                    except (OSError, ValueError, KeyError, AssertionError) as error:
                        row["roles"][role] = {"blocked": str(error) or type(error).__name__}
                validations.append(row)
            coordinate_checks[pdb] = dict(atomic_audit(path), annotated_pairs=pairs, numbering_chains=sorted(domains),
                                          checks=validations, seqres_matches_api=seqres == declared,
                                          coordinate_sha256=receipt["sha256"])

    literature = json.loads((BASE / "gpc3_primary_cells.json").read_text())
    gpc3 = {}
    for pdb, clone in literature["clones"].items():
        seq = deposited_chains(metadata[pdb])
        assert seq["C"].count(clone["alpha_cdr3"]) == seq["D"].count(clone["beta_cdr3"]) == 1
        assert seq["P"] == "FLAELAYDL" and coordinate_checks[pdb]["seqres_matches_api"]
        gpc3[pdb] = dict(clone, alpha_cdr3_offset0=seq["C"].index(clone["alpha_cdr3"]),
                        beta_cdr3_offset0=seq["D"].index(clone["beta_cdr3"]),
                        alpha_construct=seq["C"], beta_construct=seq["D"])
    a, b = gpc3["9tmm"], gpc3["9tml"]
    assert a["V"] == b["V"] and a["alpha_construct"] != b["alpha_construct"]
    assert a["beta_construct"] == b["beta_construct"] and a["beta_cdr3_nt"] != b["beta_cdr3_nt"]
    assert all(not coordinate_checks[p]["units"][0]["diagnostic_pass"] for p in gpc3)
    pmhc = [{c: deposited_chains(metadata[pdb])[c] for c in "ABP"} for pdb in gpc3]
    assert pmhc[0] == pmhc[1]

    hits = {p.stem: parse_hits(int(r["status"]), body)
            for p, r, body in verified_files(BASE / "sources/gpc3_structure")}
    motif_pairs = {name: pair_candidates(hits[f"alpha_{name}"], hits["beta_ab" if name != "c" else "beta_c"])
                   for name in "abc"}
    assert motif_pairs == {"a": ["9tmm"], "b": ["9tml"], "c": []}
    constructs = {pdb: pair_candidates(hits[f"{pdb}_3"], hits[f"{pdb}_4"]) for pdb in gpc3}
    assert all(pdbs == [pdb] for pdb, pdbs in constructs.items())

    cross = []
    for left, right in (("8i5d", "9wfs"), ("7byd", "9wnt"), ("4ftv", "9nds")):
        a, b = [metadata[p]["rcsb_primary_citation"] for p in (left, right)]
        cross.append(dict(pdbs=[left, right], study_keys=[study_key(a), study_key(b)],
                          shared_authors=sorted(author_keys({"primary_citation": a}) &
                                                author_keys({"primary_citation": b}))))
    control_dir = BASE / "sources/control"
    catalog_path = control_dir / "tcr3d_class1.html"
    catalog_receipt = next(r for r in json.loads((control_dir / "manifest.json").read_text()) if r["file"] == catalog_path.name)
    assert catalog_receipt["status"] == "200"
    catalog = json.loads(re.search(r"var data = '(\[.*?\])';", load_receipt(catalog_path, catalog_receipt).decode())[1])
    fields = ("PDB_ID", "tcr_name", "TRAV_gene", "TRBV_gene", "release_date")
    recent = [{field: row[field] for field in fields} for row in catalog if (row.get("release_date") or "") > "2026-11-01"]
    kr_comparison = {row["PDB_ID"]: {field: row[field] for field in fields} for row in catalog
                     if row["PDB_ID"] in {"8I5C", "8I5D", "9WFS", "9WFT"}}
    assert {v["TRAV_gene"] for v in kr_comparison.values()} == {"TRAV7D-2", "TRAV7D-4", "TRAV8-1"}
    for source in ("gpc3", "gpc3_supp"):
        for receipt in json.loads((BASE / "sources" / source / "manifest.json").read_text()):
            path = BASE / "sources" / source / receipt["file"]
            if path.is_file():
                load_receipt(path, receipt)
    return dict(search_coverage=discovery["counts"], statuses_in_old_design=priors,
                additional_pdbs_relative_to_old_design=fresh, coordinate_diagnostics=coordinate_checks,
                cross_source_observations=cross, natural_gpc3_controls=gpc3,
                gpc3_pmhc_constructs_equal=True, gpc3_search_entity_hits=hits,
                motif_pair_hits=motif_pairs, complete_construct_query_pair_hits=constructs,
                catalog_count=len(catalog), catalog_releases_after_previous_round=recent, kr_comparison=kr_comparison,
                SN45_companion_deposits_without_TCR=sorted(p for p in metadata if p.startswith("9wj")),
                fixed_atomic_policy_sha256=policy_sha256(), frozen_inventory=inventory, new_admissions=0,
                new_qc_cross_study_pairs=0, new_qc_same_study_natural_control_groups=0,
                protocol_frozen=False, training_started=False, independent_acceptance_accessed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    with open(args.out, "x") as stream:
        stream.write(json.dumps(audit(), indent=2, sort_keys=True) + "\n")
