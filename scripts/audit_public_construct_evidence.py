"""Replay the public clone provenance and expanded receptor discovery audit."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
from itertools import combinations
import json
from pathlib import Path
import re
import shlex
from tempfile import TemporaryDirectory
from xml.etree import ElementTree

from meta_acdc.structure.paired_variable_identity import numbered_domains
from meta_acdc.structure.rigid_controls import policy_sha256
from meta_acdc.structure.source_balanced_retrieval import author_keys
from scripts.audit_weekly_controls import coordinate_units, json_entries, verified_files
from scripts.curate_cross_source_receptors import study_key
from scripts.validate_native_search import exact_entity_pairs

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "docs/fingerprint_public_evidence_20261118"
WEEKLY_METADATA = ROOT / "docs/fingerprint_weekly_controls_20261116/sources/metadata"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_metadata():
    entries = {}
    for folder in (WEEKLY_METADATA, BASE / "broad_metadata", BASE / "metadata_retry"):
        for _, receipt, payload in verified_files(folder):
            if receipt["status"] == "200":
                entries.update(json_entries(payload))
    query = json.loads((BASE / "broad/rcsb_broad_hits.json").read_text())
    requested = {hit["identifier"].lower() for hit in query["result_set"]}
    if len(requested) != query["total_count"] or requested - entries.keys():
        raise ValueError("Broad metadata inventory is incomplete")
    return entries


def source_table(path):
    # ENA emits literal quotes in unquoted TSV titles.
    rows = list(csv.DictReader(io.StringIO(path.read_text()), delimiter="\t", quoting=csv.QUOTE_NONE))
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Malformed TSV record in {path}")
    return rows


def observed_screen(metadata):
    summary = source_table(BASE / "additional/stcrdab_summary.tsv")
    selected = [row for row in summary if row["TCRtype"] == "abTCR" and row["mhc_type"] == "MH1"
                and row["Achain"] != "NA" and row["Bchain"] != "NA"]
    groups = defaultdict(list)
    exclusions = []
    for row in selected:
        pdb = row["pdb"]
        sources = [
            ROOT / "docs/fingerprint_cross_source_extension_20261024/sources/stcrdab_all" / f"{pdb}.html",
            ROOT / "docs/fingerprint_cross_source_20261024/sources/stcrdab" / f"{pdb}.html",
        ]
        path = next((path for path in sources if path.exists()), None)
        if path is None:
            exclusions.append(dict(pdb=pdb,reason="no local numbered-domain source"))
            continue
        try:
            domains = numbered_domains(path.read_text())
            alpha, beta = domains[row["Achain"]], domains[row["Bchain"]]
            if alpha["kind"] != "VA" or beta["kind"] != "VB":
                raise ValueError("Numbered chain roles disagree with summary")
            pair = alpha["sequence"], beta["sequence"]
            if not exact_entity_pairs(metadata[pdb], *pair):
                raise ValueError("Numbered sequences absent from distinct deposited entities")
            evidence = dict(pdb=pdb, alpha_chain=row["Achain"], beta_chain=row["Bchain"],
                            annotation=str(path.relative_to(ROOT)), annotation_sha256=sha256(path))
            if evidence not in groups[pair]:
                groups[pair].append(evidence)
        except (KeyError, ValueError) as error:
            exclusions.append(dict(pdb=pdb, reason=str(error)))
    candidates = []
    for (alpha, beta), anchors in groups.items():
        hits = []
        for pdb, entry in sorted(metadata.items()):
            if not exact_entity_pairs(entry, alpha, beta):
                continue
            citation = entry.get("rcsb_primary_citation") or {}
            hits.append(dict(pdb=pdb, study=study_key(citation),
                             authors=sorted(author_keys({"primary_citation": citation}))))
        for left, right in combinations(hits, 2):
            if not left["study"] or not right["study"] or left["study"] == right["study"]:
                continue
            if not left["authors"] or not right["authors"] or set(left["authors"]) & set(right["authors"]):
                continue
            candidates.append(dict(alpha_observed=alpha, beta_observed=beta, anchors=anchors,
                                   left=left, right=right, complete_identity_unverified=True))
    unique_pairs = {tuple(sorted((row["left"]["pdb"], row["right"]["pdb"]))) for row in candidates}
    return dict(summary_rows=len(summary), summary_unique_pdbs=len({r["pdb"] for r in summary}),
                classI_pair_rows=len(selected), observed_query_pair_groups=len(groups),
                unavailable_annotation_rows=exclusions, cross_study_author_disjoint_hits=candidates,
                unique_candidate_pdb_pairs=len(unique_pairs),
                scope="Observed domain queries require complete variable-boundary adjudication.")


def complete_screen(metadata):
    plan = ROOT / "docs/fingerprint_expanded_controls_20261029/complete_search/search_plan.json"
    targets = json.loads(plan.read_text())["targets"]
    hits = []
    for target in targets:
        if not target["full_boundary_validated"]:
            raise ValueError("Incomplete anchor in complete-variable screen")
        anchor = target["anchor"]
        source = metadata[anchor].get("rcsb_primary_citation") or {}
        authors = author_keys({"primary_citation": source})
        for pdb, entry in sorted(metadata.items()):
            if pdb == anchor or not exact_entity_pairs(entry, target["sequences"]["alpha"], target["sequences"]["beta"]):
                continue
            citation = entry.get("rcsb_primary_citation") or {}
            other_authors = author_keys({"primary_citation": citation})
            if not study_key(source) or not study_key(citation) or study_key(source) == study_key(citation):
                continue
            if authors and other_authors and not authors & other_authors:
                hits.append(dict(anchor=anchor, pdb=pdb, title=entry["struct"]["title"],
                                 source_study=study_key(source), target_study=study_key(citation)))
    return dict(anchor_count=len(targets), plan_sha256=sha256(plan), hits=hits)


def assembly_checks():
    results = []
    for path, _, payload in verified_files(BASE / "candidate_assemblies"):
        pdb = path.name.split(".")[0]
        raw = (BASE / "candidate_raw" / f"{pdb}.pdb").read_text()
        lines = gzip.decompress(payload).decode().splitlines()
        atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
        identities = Counter(line[12:27] for line in atoms)
        if not identities or max(identities.values()) > 1:
            raise ValueError("Assembly atom namespace requires explicit chain-copy curation")
        with TemporaryDirectory() as temp:
            assembly = Path(temp) / "assembly.pdb"
            headers = [line for line in raw.splitlines() if line.startswith("COMPND")]
            assembly.write_text("\n".join(headers + atoms) + "\n")
            units, exclusions = coordinate_units(assembly)
        results.append(dict(pdb=pdb, n_atoms=len(atoms), units=units, exclusions=exclusions,
                            fixed_QC_pass=any(unit["coordinate_pass"] for unit in units)))
    return dict(policy_sha256=policy_sha256(), results=results)


def sample_provenance():
    runs = source_table(BASE / "retry/prjna298330.tsv")
    aliases = {}
    for path in sorted((BASE / "samples").glob("SAMN*.xml")):
        sample = ElementTree.fromstring(path.read_bytes()).find("SAMPLE")
        aliases[sample.attrib["accession"]] = sample.attrib["alias"]
    records = ElementTree.parse(BASE / "samples/3995_sra_records.xml")
    relevant = []
    rejected = []
    for package in records.iter("EXPERIMENT_PACKAGE"):
        sample = package.find("SAMPLE")
        experiment = package.find("EXPERIMENT")
        row = dict(experiment=experiment.attrib["accession"], title=experiment.findtext("TITLE"),
                   sample_alias=sample.attrib.get("alias"), sample=sample.attrib["accession"],
                   strategy=experiment.findtext(".//LIBRARY_STRATEGY"),
                   runs=[run.attrib["accession"] for run in package.iter("RUN")])
        (relevant if row["sample_alias"] == "3995" else rejected).append(row)
    project = []
    for run in runs:
        project.append(dict(run=run["run_accession"], sample=run["sample_accession"],
                            alias=aliases[run["sample_accession"]], strategy=run["library_strategy"]))
    rna_packages = ElementTree.parse(BASE / "samples/prjna298330_sra.xml")
    rna_aliases = [sample.attrib.get("alias") for sample in rna_packages.iter("SAMPLE")]
    return dict(project="PRJNA298330", project_runs=project,
                strategy_counts=dict(Counter(row["strategy"] for row in project)),
                RNA_Seq_aliases_verified_in_primary_SRA=rna_aliases,
                exact_3995_SRA_records=relevant, irrelevant_all_field_query_records=rejected,
                adjacent_project_run_counts={name: len(source_table(BASE / "history" / f"{name}.tsv"))
                                             for name in ("prjna298310", "prjna298376")},
                original_9d_paired_reads_recovered=False)


def revision_history():
    results = {}
    for path, _, payload in verified_files(BASE / "versions"):
        pdb = path.name[:4]
        text = gzip.decompress(payload).decode()
        start = text.index("_pdbx_audit_revision_history.ordinal")
        block = text[start:text.index("\n#", start)]
        data = [shlex.split(line) for line in block.splitlines() if not line.startswith("_")]
        if any(len(row) != 5 for row in data):
            raise ValueError("Unexpected revision-history columns")
        results[pdb] = dict(history=[dict(ordinal=int(row[0]),major=int(row[2]),
                                               minor=int(row[3]),date=row[4]) for row in data],
                            available_coordinate_file=path.name,
                            only_major_one=all(row[2] == "1" for row in data),
                            initial_release_bytes_available=False)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    receipts = json.loads((BASE / "retained_source_manifest.json").read_text())
    for record in receipts:
        path = ROOT / record["path"]
        if sha256(path) != record["sha256"]:
            raise ValueError(f"Source checksum mismatch for {path}")
    metadata = load_metadata()
    output = dict(metadata_entries=len(metadata), observed_discovery=observed_screen(metadata),
                  complete_variable_discovery=complete_screen(metadata), candidate_assemblies=assembly_checks(),
                  sample_provenance=sample_provenance(), revision_history=revision_history(),
                  training_started=False, protocol_frozen=False, independent_acceptance_accessed=False)
    with args.out.open("x", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
