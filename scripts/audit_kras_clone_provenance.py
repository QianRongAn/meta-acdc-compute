"""Compare the disclosed clone with the crystal construct and bounded public reads."""
import argparse
from collections import Counter, defaultdict
import csv
from difflib import SequenceMatcher
from itertools import combinations, product
import gzip
import hashlib
import json
from pathlib import Path
import re

from meta_acdc.structure.graph import AA3TO1
from meta_acdc.structure.paired_variable_identity import _CODONS, fasta_records
from meta_acdc.structure.receptor_identity_audit import seqres_sequences

from scripts.audit_weekly_controls import json_entries, verified_files
from scripts.join_terminal_rescue import load_receipt
from scripts.curate_cross_source_receptors import study_key
from scripts.join_native_search import parse_hits
from meta_acdc.structure.source_balanced_retrieval import author_keys

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "docs/fingerprint_kras_clones_20261117"
WEEKLY = ROOT / "docs/fingerprint_weekly_controls_20261116"
DNA = set("ACGT")
RC = str.maketrans("ACGTN", "TGCAN")


def translate(sequence):
    return "".join(_CODONS[sequence[i:i + 3]] for i in range(0, len(sequence) - 2, 3))


def clone_alignment():
    rows, dna = {}, {}
    receipts = {(r["clone"], r["id"]): r for r in json.loads((BASE / "patent_listing_receipts.json").read_text())}
    for header, seq in fasta_records((BASE / "clone_sequences.fasta").read_text()):
        clone, role = header.split("_")[:2]
        pdb = {"9a": "6ulr", "9d": "6uln"}[clone]
        chain = {"alpha": "D", "beta": "E"}[role]
        if header.endswith("DNA"):
            seq = seq.upper()
            assert len(seq) == {"alpha": 388, "beta": 400}[role]
            origin = receipts[(clone, {"alpha": 63, "beta": 64}[role])]
            assert seq == origin["sequence"]
            dna[role] = seq
            sequence = translate(seq)
        else:
            sequence = seq
            origin = receipts[(clone, {"alpha": 9, "beta": 10}[role])]
            assert "".join(AA3TO1[token.upper()] for token in origin["residue_tokens"]) == sequence
        expected_start = {"alpha": "LAKTT", "beta": "DAGVT"}[role]
        start = sequence.index(expected_start)
        segment = sequence[start:]
        cdr3 = {"alpha": "CLVGDMDQAGTALIF", "beta": {
            "9a": "CASSLGEGRVDGYTF", "9d": "CASSLGQTNYGYTF"}[clone]}[role]
        assert len(sequence) == {"9a": (129, 133), "9d": (129, 132)}[clone][role == "beta"]
        numbered = fasta_records((WEEKLY / "sources/followup" / f"{pdb}_{chain}.fasta").read_text())[0][1]
        pdb_text = (WEEKLY / "sources/coordinates" / f"{pdb}.pdb").read_text()
        full = seqres_sequences(pdb_text)[chain]
        assert sequence.count(cdr3) == numbered.count(cdr3) == 1
        assert full.count(numbered) == 1
        offset = full.index(numbered)
        deposited_segment = full[offset:offset + len(segment)]
        edits = [dict(operation=op, clone_bounds0=[i, j], pdb_bounds0=[k, m],
                      clone=segment[i:j], pdb=numbered[k:m])
                 for op, i, j, k, m in SequenceMatcher(a=segment, b=numbered, autojunk=False).get_opcodes()
                 if op != "equal"]
        padded_edits = [dict(position1=i + 1, clone=a, pdb=b)
                        for i, (a, b) in enumerate(zip(segment, deposited_segment, strict=True)) if a != b]
        assert padded_edits == ([dict(position1=10, clone="M", pdb="V")] if role == "alpha"
                                else [dict(position1=1, clone="D", pdb="M")])
        author_pos = offset + padded_edits[0]["position1"]
        modeled = sorted({line[17:20] for line in pdb_text.splitlines()
                          if line.startswith("ATOM  ") and line[21] == chain and line[22:26].strip() == str(author_pos)})
        assert modeled == (["VAL"] if role == "alpha" else [])
        rows[clone + "_" + role] = dict(
            difference_author_residue=author_pos, difference_modeled_residue_names=modeled,
            pdb=pdb, chain=chain, seq_id=header, clone_protein_length=len(sequence),
            leader_length_in_disclosure=start, clone_mature_segment=segment,
            numbered_segment=numbered, deposited_segment=deposited_segment,
            deposited_prefix_outside_numbered=full[:offset], numbered_segment_edits=edits,
            deposited_equal_length_segment_edits=padded_edits, cdr3=cdr3,
        )
    return rows, dna


def screen_reads(path, rows, dna):
    """Preserve observed codons and tie separate reads by exact barcode and UMI."""
    targets = {}
    for role, sequence in dna.items():
        r = rows["9a_" + role]
        protein = translate(sequence)
        begin = protein.index(r["cdr3"]) * 3
        end = begin + 3 * len(r["cdr3"])
        targets[role + "_cdr3"] = (sequence[begin:end], None)
        start = (r["leader_length_in_disclosure"] + (9 if role == "alpha" else 0)) * 3
        targets[role + "_start"] = (sequence[start - 12:start + 21], 12)
    molecules = defaultdict(lambda: defaultdict(set))
    witnesses = defaultdict(dict)
    spot_counts = Counter()
    total = 0
    with gzip.open(path, "rt") as stream:
        for total, record in enumerate(stream, 1):
            bases, lengths, quals = record.rstrip("\n").split("\t")
            lengths = [int(value) for value in lengths.split(",")]
            assert lengths == [150, 150, 8] and len(bases) == 308
            if set(bases[:26]) - DNA:
                continue
            barcode, umi = bases[:16], bases[16:26]
            quality = [int(value) for value in quals.split(",")]
            assert len(quality) == len(bases)
            if min(quality[:26]) < 20:
                continue
            spot_counts["Q20_barcode_UMI"] += 1
            observations = set()
            for left, right in ((26, 150), (150, 300)):
                for observed, scores in ((bases[left:right], quality[left:right]),
                                         (bases[left:right].translate(RC)[::-1], quality[left:right][::-1])):
                    for name, (target, codon) in targets.items():
                        if codon is None:
                            pos = observed.find(target)
                            end = pos + len(target)
                            if pos >= 0 and min(scores[pos:end]) >= 20:
                                observations.add(name)
                            continue
                        pos = observed.find(target[:codon])
                        if pos < 0 or not observed.startswith(target[codon + 3:], pos + codon + 3):
                            continue
                        start = pos + codon
                        variant = observed[start:start + 3]
                        if variant in _CODONS and min(scores[pos:pos + len(target)]) >= 20 and min(scores[start:start + 3]) >= 30:
                            observations.add(name + ":" + variant)
            for observation in observations:
                spot_counts[observation] += 1
                molecules[barcode][observation].add(umi)
                witnesses[barcode].setdefault(observation, total)
    paired = []
    for barcode, features in sorted(molecules.items()):
        if not {"alpha_cdr3", "beta_cdr3"} <= features.keys():
            continue
        signature = hashlib.sha256(barcode.encode()).hexdigest()
        paired.append(dict(barcode_sha256=signature,
                           distinct_umis={k: len(v) for k, v in sorted(features.items())},
                           first_spot=witnesses[barcode]))
    return dict(total_spots=total, counts=dict(spot_counts), targets=targets, paired_barcode_observations=paired,
                scope="Exact raw barcode co-occurrence and targeted codons, without whole transcript assembly.")


def structure_search(rows):
    catalog_text = (WEEKLY / "sources/discovery/class1.html").read_text()
    catalog = json.loads(re.search(r"var data = '(\[.*?\])';", catalog_text)[1])
    metadata = {}
    for _, receipt, body in verified_files(WEEKLY / "sources/metadata"):
        assert receipt["status"] == "200"
        metadata.update(json_entries(body))
    groups = defaultdict(list)
    for row in catalog:
        pdb = row["PDB_ID"].lower()
        alpha, beta = [], []
        for entity in metadata[pdb]["polymer_entities"]:
            label = entity["rcsb_polymer_entity"]["pdbx_description"] or ""
            if not re.search(r"TCR|T.cell", label, re.I):
                continue
            sequence = entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
            if re.search(r"alpha|\btra[vjc]|v.alpha", label, re.I):
                alpha.append(sequence)
            if re.search(r"beta|\btrb[vjc]|v.beta", label, re.I):
                beta.append(sequence)
        for pair in set(product(alpha, beta)):
            groups[pair].append(pdb)
    independent = set()
    for members in groups.values():
        for left, right in combinations(members, 2):
            a, b = [metadata[p]["rcsb_primary_citation"] or {} for p in (left, right)]
            studies = study_key(a), study_key(b)
            aa, bb = [author_keys({"primary_citation": source}) for source in (a, b)]
            if all(studies) and studies[0] != studies[1] and aa and bb and not aa & bb:
                independent.add(tuple(sorted((left, right))))
    current, clone_hits = {}, {}
    survey = BASE / "survey"
    for receipt in json.loads((survey / "manifest.json").read_text()):
        path = survey / receipt["file"]
        if path.stem not in rows and path.stem not in ("new_release", "new_revision"):
            continue
        body = load_receipt(path, receipt)
        if path.stem in rows:
            clone_hits[path.stem] = parse_hits(int(receipt["status"]), body)
        elif path.stem in ("new_release", "new_revision"):
            assert receipt["status"] == "204" and not body
            current[path.stem] = []
    assert clone_hits.keys() == rows.keys() and len(current) == 2
    return dict(catalog_entries=len(catalog), exact_full_deposited_pair_groups=len(groups),
                author_disjoint_distinct_study_pairs=sorted(independent),
                natural_variable_sequence_local_search_hits=clone_hits, date_hits=current,
                scope="Description-based full deposited chain screen. Changed constant domains can hide variable-sequence matches.")


def project_availability():
    output = {}
    for label, path in (("GSE136394", BASE / "followup/ena_runs.txt"),
                        ("PRJNA342632", BASE / "final_sources/tran_project_runs.txt")):
        with path.open() as stream:
            records = list(csv.DictReader(stream, delimiter="\t"))
        output[label] = dict(runs=len(records),
                             titled_4095=sum("4095" in row["experiment_title"] for row in records),
                             titled_3995=sum("3995" in row["experiment_title"] for row in records))
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows, dna = clone_alignment()
    result = dict(clone_alignments=rows, structure_search=structure_search(rows), public_projects=project_availability(),
                  protocol_frozen=False, independent_acceptance_accessed=False)
    if args.reads:
        result["single_cell_evidence"] = screen_reads(args.reads, rows, dna)
    with args.out.open("x") as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
