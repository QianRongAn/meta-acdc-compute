"""Replay the 9d sequence-listing and independent-structure follow-up."""
import argparse
from collections import Counter
from itertools import product
import gzip
import hashlib
import json
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from xml.etree import ElementTree
from zipfile import ZipFile

from meta_acdc.structure.paired_variable_identity import fasta_records
from meta_acdc.structure.receptor_identity_audit import seqres_sequences
from meta_acdc.structure.rigid_controls import policy_sha256
from meta_acdc.structure.source_balanced_retrieval import author_keys
from scripts.audit_kras_clone_provenance import BASE as PREVIOUS, ROOT, WEEKLY, clone_alignment, translate
from scripts.audit_weekly_controls import coordinate_units, json_entries, verified_files
from scripts.curate_cross_source_receptors import study_key
from scripts.validate_native_search import exact_entity_pairs

BASE = ROOT / "docs/fingerprint_kras_crystal_constructs_20261117"


def sequence_results():
    page = (BASE / "au_page82.txt").read_text(encoding="utf-8")
    pdf = PREVIOUS / "followup/patent_9d_au.pdf"
    pdf_hash = hashlib.sha256(pdf.read_bytes()).hexdigest()
    if pdf_hash != "c345ad5ff130af345d0f919f4183b5a831530c79e47b6706bf2611c779810f52":
        raise ValueError("Unexpected AU patent PDF")
    prior = dict(fasta_records((PREVIOUS / "clone_sequences.fasta").read_text()))
    rows, prior_dna = clone_alignment()
    result = {}
    for role, number, protein_id, size in (("alpha", 28, 9, 387), ("beta", 29, 10, 396)):
        block = re.search(rf"^<400>\s+{number}\s*\n(.*?)^<210>", page, re.M | re.S)[1]
        seq = ""
        for line in block.splitlines():
            bases, total = re.fullmatch(r"([acgtn ]+)\s+(\d+)", line).groups()
            seq += bases.replace(" ", "").upper()
            if len(seq) != int(total):
                raise ValueError("Patent listing cumulative length mismatch")
        if len(seq) != size:
            raise ValueError("Patent listing record length mismatch")
        reference = next(seq for k, seq in prior.items() if k.startswith(f"9d_{role}_"))
        indices = [i for i, residue in enumerate(seq) if residue == "N"]
        translations = []
        for values in product("ACGT", repeat=len(indices)):
            edited = list(seq)
            for i, base in zip(indices, values, strict=True):
                edited[i] = base
            protein = translate("".join(edited))
            if len(protein) != len(reference):
                raise ValueError("Translated and protein-listed lengths do not match")
            translations.append(dict(
                choices=dict(zip((i + 1 for i in indices), values, strict=True)),
                complete_protein_exact=protein == reference,
                affected_codons=[dict(nucleotides="".join(edited[i // 3 * 3:i // 3 * 3 + 3]),
                                      amino_acid=protein[i // 3], protein_position1=i // 3 + 1,
                                      reference_amino_acid=reference[i // 3]) for i in indices],
            ))
        if not any(v["complete_protein_exact"] for v in translations):
            raise ValueError("No nucleotide interpretation matches disclosed protein")
        original = rows[f"9d_{role}"]
        leader = original["leader_length_in_disclosure"]
        start = 3 * (leader + (9 if role == "alpha" else 0))
        result[role] = dict(seq_id=number, protein_seq_id=protein_id, dna=seq, length=size,
                            unresolved_positions1=[i + 1 for i in indices], translations=translations,
                            mature_N_terminal_checked_codon=seq[start:start + 3],
                            crystal_delta=original["deposited_equal_length_segment_edits"],
                            crystal_author_position=original["difference_author_residue"],
                            crystal_modeled_residue_names=original["difference_modeled_residue_names"],
                            entire_nucleotide_prefix_matches_9a=prior_dna[role].startswith(seq))
    return dict(patent_pdf_sha256=pdf_hash, pdf_page1=82, listing_labels="Artificial Sequence; Synthetic",
                original_paired_reads_recovered=False, chains=result)


def source_check():
    paths = {}
    receipts = json.loads((BASE / "retrieval_receipts.json").read_text())
    for r in receipts:
        if r.get("retained_path"):
            path = BASE / r["retained_path"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != r["sha256"]:
                raise ValueError(f"Source hash mismatch: {path}")
            paths[path.name] = path
    xml = ElementTree.parse(paths["elife.xml"])
    methods = [" ".join("".join(e.itertext()).split()) for e in xml.iter("p")]
    constructs = [text for text in methods if "pET30a via NdeI and XhoI (GenScript)" in text]
    models = [text for text in methods if "as a search model (PDB: 6ULR)" in text]
    if len(constructs) != 1 or len(models) != 1:
        raise ValueError("Missing construct or search-model methods")
    reagents = {}
    for row in xml.iter("tr"):
        cells = [" ".join("".join(e.itertext()).split()) for e in row]
        if len(cells) > 3 and cells[1] in ("TCR9a alpha chain", "TCR9a beta chain"):
            reagents[cells[1]] = dict(source=cells[2], identifier=cells[3] or None)
    if set(reagents) != {"TCR9a alpha chain", "TCR9a beta chain"}:
        raise ValueError("Missing α or β key reagent")
    with ZipFile(paths["elife_supp.docx"]) as archive:
        content = ElementTree.fromstring(archive.read("word/document.xml"))
    supplement = " ".join(e.text.strip() for e in content.iter() if e.tag.endswith("}t") and e.text)
    if not all(x in supplement for x in ("7SU9", "Data collection", "Refinement")):
        raise ValueError("Unexpected eLife supplement")
    latest = json.loads((PREVIOUS / "audit_final.json").read_text())
    old = seqres_sequences((WEEKLY / "sources/coordinates/6ulr.pdb").read_text())
    new = seqres_sequences((WEEKLY / "sources/cross/7su9.pdb").read_text())
    return dict(reagents=reagents, followup_PDB="7su9", molecular_replacement_model="6ulr",
                supplementary_file_character_count=len(supplement),
                supplementary_listed_construct_keywords=[k for k in ("pET30a", "GenScript", "NdeI", "plasmid")
                                                          if k.casefold() in supplement.casefold()],
                identical_6ulr_7su9_TCR_SEQRES={chain: bool(old[chain]) and old[chain] == new[chain] for chain in "DE"},
                upstream_release_search_reused_from="fingerprint_kras_clones_20261117/audit_final.json",
                upstream_release_keys=[k for k in latest if "search" in k])


def validate_assemblies(pdbs):
    receipts = {r["file"]: r for r in json.loads((BASE / "retrieval_receipts.json").read_text())}
    former = ROOT / "docs/fingerprint_np_adjudication_20261029/sources/followup"
    result = {}
    for pdb in pdbs:
        r = receipts[f"{pdb}.pdb1.gz"]
        data = (BASE / r["retained_path"]).read_bytes()
        if hashlib.sha256(data).hexdigest() != r["sha256"]:
            raise ValueError("Assembly checksum mismatch")
        atoms = lambda lines: [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
        assembly = atoms(gzip.decompress(data).decode().splitlines())
        previous = gzip.decompress((former / f"{pdb}.pdb.gz").read_bytes()).decode().splitlines()
        original_ids = Counter(line[12:27] for line in atoms(previous))
        assembly_ids = Counter(line[12:27] for line in assembly)
        if original_ids != assembly_ids or max(assembly_ids.values()) != 1:
            raise ValueError("Assembly does not preserve each deposited atom exactly once")
        with TemporaryDirectory() as temp:
            path = Path(temp) / "assembly.pdb"
            headers = [line for line in previous if line.startswith("COMPND")]
            path.write_text("\n".join(headers + assembly) + "\n")
            units, exclusions = coordinate_units(path)
        result[pdb] = dict(atom_count=len(assembly), atom_identities_match_prior=True,
                           fixed_QC_pass=any(row["coordinate_pass"] for row in units),
                           graph_failures=[row["graph_error"] for row in units if "graph_error" in row],
                           unit_exclusions=exclusions)
    return result


def structure_screen():
    plan = ROOT / "docs/fingerprint_expanded_controls_20261029/complete_search/search_plan.json"
    targets = json.loads(plan.read_text())["targets"]
    metadata = {}
    for _, receipt, payload in verified_files(WEEKLY / "sources/metadata"):
        if receipt["status"] != "200":
            raise ValueError("Incomplete weekly metadata")
        metadata.update(json_entries(payload))
    classifications, leads = Counter(), []
    for t in targets:
        if not t["full_boundary_validated"]:
            raise ValueError("Incomplete variable anchor")
        anchor = t["anchor"]
        # Some older complete anchors lie outside the current metadata snapshot.
        prior_citation = (metadata.get(anchor, {}).get("rcsb_primary_citation") or {})
        for pdb, entry in sorted(metadata.items()):
            pairs = exact_entity_pairs(entry, t["sequences"]["alpha"], t["sequences"]["beta"])
            if not pairs or pdb == anchor:
                continue
            citation = entry.get("rcsb_primary_citation") or {}
            study, source = study_key(citation), t["anchor_study"]
            left, right = (author_keys({"primary_citation": item}) for item in (prior_citation, citation))
            if not source or not study:
                status = "unknown_study"
            elif source == study:
                status = "same_study"
            elif not left or not right or study_key(prior_citation) != source:
                status = "different_study_author_unresolved"
            elif right & left:
                status = "different_study_overlapping_authors"
            else:
                status = "different_study_disjoint_authors_requires_coordinate_and_identity_gates"
            classifications[status] += 1
            if status != "same_study":
                leads.append(dict(anchor=anchor, pdb=pdb, anchor_study=source, study=study,
                                  status=status, shared_author_keys=sorted(left & right),
                                  entity_pairs=[[a["rcsb_id"], b["rcsb_id"]] for a, b in pairs]))
    for lead in leads:
        if lead["status"] != "different_study_disjoint_authors_requires_coordinate_and_identity_gates":
            continue
        pdbs = [lead["anchor"], lead["pdb"]]
        lead["official_assembly_QC"] = validate_assemblies(pdbs)
        if all(row["fixed_QC_pass"] for row in lead["official_assembly_QC"].values()):
            raise ValueError("New twice-passing candidate requires full independent identity adjudication")
    return dict(anchor_count=len(targets), metadata_entries=len(metadata), classifications=dict(classifications),
                non_same_study_leads=leads, fixed_atomic_policy_sha256=policy_sha256(),
                screen_is_not_identity_or_atom_QC=True, n_new_qualified_cross_study_pairs=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    data = dict(sequence_listing=sequence_results(), crystallography=source_check(), structure_scan=structure_screen(),
                independent_acceptance_accessed=False, protocol_frozen=False, model_training_started=False)
    with args.out.open("x", encoding="utf-8") as output:
        json.dump(data, output, indent=2, sort_keys=True)
        output.write("\n")


if __name__ == "__main__":
    main()
