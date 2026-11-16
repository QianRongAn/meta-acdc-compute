"""Replay the public PDB release audit and fixed-gate candidate checks."""
import argparse
from collections import Counter, defaultdict
from dataclasses import replace
import gzip
import json
from pathlib import Path
import re
import tempfile
from urllib.parse import parse_qs, urlsplit

from meta_acdc.structure.graph import build_interface_graph_from_residues, parse_pdb
from meta_acdc.structure.identity import annotate, class_i_units, missing_heavy_atoms
from meta_acdc.structure.paired_variable_identity import fasta_records, germline_references, numbered_domains, validate_domain
from meta_acdc.structure.primary_role_evidence import primary_roles
from meta_acdc.structure.receptor_identity_audit import is_subsequence, seqres_sequences
from meta_acdc.structure.rigid_controls import atomic_screen, policy_sha256
from meta_acdc.structure.source_balanced_retrieval import author_keys
from scripts.audit_terminal_rescue import deposited_chains, verified_files
from scripts.curate_cross_source_receptors import chain_pairs, study_key
from scripts.join_native_search import ROOT, pair_candidates, parse_hits
from scripts.validate_native_search import exact_entity_pairs

BASE = ROOT / "docs/fingerprint_weekly_controls_20261116"
TCR = re.compile(r"\btcr\b|t[ -]+(?:cell|lymphocyte)[ -]+receptor|\btra[cv]\b|\btrb[cv]\b", re.I)
MHC = re.compile(r"histocompatibility|\bhla\b|\bmhc\b|microglobulin", re.I)
FIELDS = ("PDB_ID", "tcr_name", "TRAV_gene", "TRBV_gene", "release_date")
CROSS = (("5swz", "7jwi"), ("6ulr", "7su9"), ("9pbg", "9pbh"))


def json_entries(data):
    response = json.loads(data)
    assert not response.get("errors")
    assert all(response["data"]["entries"])
    return {entry["rcsb_id"].lower(): entry for entry in response["data"]["entries"]}


def entry_ids(data):
    response = json.loads(data)
    assert response["result_type"] == "entry"
    assert len(response["result_set"]) == response["total_count"]
    return {hit["identifier"].lower() for hit in response["result_set"]}


def descriptions(entry):
    return [x["rcsb_polymer_entity"]["pdbx_description"] or "" for x in entry["polymer_entities"]]


def coordinate_units(path):
    residues = parse_pdb(path)
    annotate(residues, path)
    units, exclusions = class_i_units(residues)
    results = []
    for unit, meta in units:
        row = dict(chains=meta)
        try:
            graph = build_interface_graph_from_residues([replace(r, plddt=None) for r in unit])
            row["graph_nodes"] = graph.n_nodes
        except ValueError as error:
            row["graph_error"] = str(error)
        try:
            row["atomic_screen"] = atomic_screen(unit)
        except ValueError as error:
            row["atomic_screen_error"] = str(error)
            row["peptide_missing_atoms"] = [
                [r.chain, r.resid, r.icode, atoms] for r in unit
                if r.role == "peptide" and (atoms := missing_heavy_atoms(r))
            ]
        row["coordinate_pass"] = bool(row.get("graph_nodes") and row.get("atomic_screen", {}).get("passed"))
        results.append(row)
    return results, exclusions


def audit():
    sources = BASE / "sources"
    files, metadata, hashes, sequence_results = {}, {}, {}, {}
    for folder in sorted(sources.iterdir()):
        for path, receipt, payload in verified_files(folder):
            hashes[str(path.relative_to(BASE))] = receipt["sha256"]
            if folder.name == "sequence":
                query = json.loads(parse_qs(urlsplit(receipt["url"]).query)["json"][0])
                sequence_results[path.name] = (query, parse_hits(int(receipt["status"]), payload))
                continue
            assert receipt["status"] == "200", receipt
            if path.name.startswith("batch_"):
                metadata.update(json_entries(payload))
            else:
                assert path.name not in files
                files[path.name] = path
    releases = entry_ids(files["release.json"].read_bytes())
    revisions = entry_ids(files["revision.json"].read_bytes())
    catalog = json.loads(re.search(r"var data = '(\[.*?\])';", files["class1.html"].read_text())[1])
    recent = json.JSONDecoder().raw_decode(files["recent.html"].read_text().split("var tableData = ", 1)[1])[0]
    catalog_ids = {row["PDB_ID"].lower() for row in catalog}
    assert releases | catalog_ids | {"26jd", "7jwi", "7jwj"} == metadata.keys()
    release_audit = {}
    for pdb in sorted(releases):
        entry = metadata[pdb]
        title = (entry.get("struct") or {}).get("title") or ""
        entities = descriptions(entry)
        if TCR.search(title) or any(TCR.search(d) or MHC.search(d) for d in entities):
            release_audit[pdb] = dict(title=title, entities=entities, tcr_marker=bool(TCR.search(title) or any(TCR.search(d) for d in entities)))
    assert not any(r["tcr_marker"] for r in release_audit.values())

    groups = defaultdict(list)
    for row in catalog:
        a, b, name = (row[k] for k in ("TRAV_gene", "TRBV_gene", "tcr_name"))
        citation = metadata[row["PDB_ID"].lower()].get("rcsb_primary_citation") or {}
        source = study_key(citation) or citation.get("title")
        if a and b and name and source:
            groups[a, b, source].append({key: row[key] for key in FIELDS})
    potential_controls = [dict(V=[a, b], source=source, entries=rows)
                          for (a, b, source), rows in sorted(groups.items())
                          if len({r["tcr_name"] for r in rows}) > 1]
    refs = germline_references(sorted((ROOT / "docs/fingerprint_paired_identity_20261024/sources").glob("*_TR*.fasta")))
    coordinates, anchors = {}, {}
    for name, path in sorted(files.items()):
        if path.suffix != ".pdb":
            continue
        pdb = path.stem
        header = seqres_sequences(path.read_text())
        deposit = deposited_chains(metadata[pdb])
        text = files[f"{pdb}.html"].read_text()
        domains = numbered_domains(text)
        roles = primary_roles(path)
        units, exclusions = coordinate_units(path)
        checked_pairs = []
        for alpha, beta in chain_pairs(text):
            pair = {}
            for role, chain, kind in (("alpha", alpha, "VA"), ("beta", beta, "VB")):
                try:
                    assert roles[chain]["role"] == role and domains[chain]["kind"] == kind
                    assert header[chain] == deposit[chain]
                    info = validate_domain(header[chain], domains[chain],
                                           files[f"{pdb}_{chain}.fasta"].read_text(), text, refs)
                    pair[role] = {"sequence": info["variable"], "V_allele": info["v_allele"]}
                except (KeyError, ValueError, AssertionError) as error:
                    pair[role] = {"blocked": str(error) or type(error).__name__}
            checked_pairs.append(dict(alpha=alpha, beta=beta, roles=pair))
            if alpha in domains and beta in domains:
                seqs = [domains[c]["sequence"] for c in (alpha, beta)]
                anchors.setdefault(pdb, seqs)
                assert anchors[pdb] == seqs
        for unit in units:
            chains = unit["chains"]
            unit["observed_TCR_atoms_match_SEQRES"] = all(
                is_subsequence(chains[role], header[chains[role + "_chain"]])
                for role in ("alpha", "beta")
            )
        coordinates[pdb] = dict(seqres_matches_api=header == deposit, identity_checks=checked_pairs,
                                units=units, scope_exclusions=exclusions)
    anchors["9pbg"] = [deposited_chains(metadata["9pbg"])[ch] for ch in "DE"]
    matches = {p: sorted(q for q, row in metadata.items() if exact_entity_pairs(row, *seqs)) for p, seqs in anchors.items()}
    sequence_plan = json.loads((BASE / "sequence_plan.json").read_text())
    sequence_hits = {}
    for pdb, (alpha, beta) in sequence_plan.items():
        a, b = (sequence_results[key] for key in (alpha, beta))
        for (query, _), sequence in zip((a, b), anchors[pdb], strict=True):
            params = query["query"]["parameters"]
            assert params == dict(identity_cutoff=1.0, evalue_cutoff=.1, sequence_type="protein", value=sequence)
            assert query["return_type"] == "polymer_entity" and query["request_options"]["return_all_hits"]
        sequence_hits[pdb] = pair_candidates(a[1], b[1])

    kr_n_termini = {}
    for pdb in ("6uln", "6ulr", "7su9"):
        rows = {}
        for role, chain in (("alpha", "D"), ("beta", "E")):
            header, observed = fasta_records(files[f"{pdb}_{chain}.fasta"].read_text())[0]
            allele = header.split("|")[-1]
            full = deposited_chains(metadata[pdb])[chain]
            references = sorted({r["amino_acids"][:10] for r in refs["human", allele]})
            assert full.count(observed) == 1 and references and observed[:10] not in references
            rows[role] = dict(allele=allele, observed_start=observed[:10], reference_starts=references,
                              domain_offset0=full.index(observed), deposited_start=full[:13])
        kr_n_termini[pdb] = rows
    assemblies = {}
    for name in ("6uln", "6ulr", "7su9"):
        original = files[name + ".pdb"].read_text().splitlines()
        official = gzip.decompress(files[name + ".pdb1.gz"].read_bytes()).decode().splitlines()
        is_atom = lambda line: line.startswith(("ATOM  ", "HETATM"))
        transformed = [line for line in official if is_atom(line)]
        # In these three deposits the models partition distinct chains/atoms.
        # Keep all deposited coordinates and require the exact same atom identities.
        original_atoms = Counter(line[12:27] for line in original if is_atom(line))
        assembly_atoms = Counter(line[12:27] for line in transformed)
        assert original_atoms == assembly_atoms and max(assembly_atoms.values()) == 1
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "assembly.pdb"
            path.write_text("\n".join([line for line in original if line.startswith("COMPND")] + transformed) + "\n")
            units, exclusions = coordinate_units(path)
        assemblies[name] = dict(units=units, exclusions=exclusions, atom_identities_equal=True, atoms=len(transformed))
    assert all(assemblies[p]["units"][0]["coordinate_pass"] for p in ("6uln", "6ulr"))
    assert not assemblies["7su9"]["units"]
    assert len(coordinates["5swz"]["units"]) == 4
    assert all(not u["coordinate_pass"] for u in coordinates["5swz"]["units"])
    assert all(u["coordinate_pass"] for u in coordinates["5sws"]["units"])
    assert not any({r["V_allele"].split("*")[0] for r in pair["roles"].values() if "V_allele" in r}
                   == {"TRAV14-1", "TRBV17"} for pair in coordinates["7jwj"]["identity_checks"])
    cross = []
    for a, b in CROSS:
        ca, cb = [metadata[p]["rcsb_primary_citation"] for p in (a, b)]
        ad, bd = [deposited_chains(metadata[p]) for p in (a, b)]
        cross.append(dict(pdbs=[a, b], study_keys=[study_key(ca), study_key(cb)],
                          shared_authors=sorted(author_keys({"primary_citation": ca}) & author_keys({"primary_citation": cb})),
                          entire_deposited_chain_equal={role: ad[ch] == bd[ch] for role, ch in (("alpha", "D"), ("beta", "E"))},
                          other_entry_contains_both_query_sequences=b in matches.get(a, [])))
    prior = ROOT / "docs/fingerprint_terminal_rescue_20261101"
    checked = json.loads((prior / "search_results.json").read_text())
    known = {p for row in checked["rows"] for p in row["exact_entity_matches"]}
    revision_rows = {p: json.loads(files[p + ".json"].read_text())["rcsb_accession_info"]
                     for p in ("9tml", "9tmm", "7rtr", "7n1f", "4jfe", "6eqb", "9nu7", "9nv2", "9nv6", "9nve")}
    return dict(release_entries=len(releases), revision_entries=len(revisions),
                catalog_count=len(catalog), recent_catalog=[r for r in recent if r["release_date"] >= "2026-11-04"],
                release_metadata_complete=True, release_target_like_entities=release_audit,
                known_fragment_PDBs_revised=sorted(known & revisions),
                class1_catalog_revisions=sorted(catalog_ids & revisions),
                priority_versions=revision_rows, additional_scope_26JD=descriptions(metadata["26jd"]),
                potential_same_study_same_V_name_groups=potential_controls,
                coordinate_checks=coordinates, official_assembly_checks=assemblies, observed_sequence_matches=matches,
                same_receptor_checks=cross, RCSB_sequence_pair_hits=sequence_hits, sequence_query_count=len(sequence_results),
                KRAS_N_terminal_evidence=kr_n_termini,
                fixed_policy_sha256=policy_sha256(), source_hashes=hashes,
                new_qualified_cross_study_pairs=0, new_qualified_natural_control_groups=0,
                protocol_frozen=False, training_started=False, independent_acceptance_accessed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit()
    with args.out.open("x") as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
