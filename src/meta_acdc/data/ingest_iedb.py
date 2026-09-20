"""Normalize IEDB mhc_ligand_full CSV into the unified schema.

The current IEDB export uses a TWO-ROW header (group row + column-name row),
e.g. groups: Assay, Epitope, Host, MHC Restriction, Reference, ...

Rows kept: MHC ligand / binding assays (Method contains mhc/binding/mass
spectrometry/elution, NOT T cell assays) on human data. T cell assay rows are
excluded — activation data comes from VDJdb/ACDC (see schema.py rationale).

Labels: Qualitative Measurement Positive -> 1, Negative -> 0.
Score: -log10(IC50 nM) from Quantitative measurement when parseable.
Streaming: the CSV is ~9 GB; rows are deduplicated on (peptide, allele) with a
seen-set and written incrementally (memory-safe).

Usage:
    python ingest_iedb.py --in mhc_ligand_full.csv --out processed/iedb.clean.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

from schema import DataSource, DatasetStats, InteractionRecord, SCHEMA_TSV

METHOD_KEEP = ("mhc", "mass spectr", "binding", "elution")
METHOD_DROP = ("t cell",)


def parse_ic50(raw: str) -> float | None:
    try:
        v = float(raw)
        return -math.log10(max(v, 1e-12))
    except (TypeError, ValueError):
        return None


def ingest(path: Path, out_path: Path, max_rows: int | None = None) -> DatasetStats:
    stats = DatasetStats()
    seen: set[tuple[str, str]] = set()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, newline="", encoding="utf-8", errors="replace") as fh, \
         open(out_path, "w", newline="") as out:
        raw = csv.reader(fh)
        h1 = next(raw)
        h2 = next(raw)
        fieldnames = [f"{g}|{n}" for g, n in zip(h1, h2)]

        # positionally index the fields we need (faster than DictReader on 9GB)
        def idx(group: str, name: str) -> int:
            return fieldnames.index(f"{group}|{name}")

        i_method = idx("Assay", "Method")
        i_epitope = idx("Epitope", "Name")
        i_allele = idx("MHC Restriction", "Name")
        i_qual = idx("Assay", "Qualitative Measurement")
        i_quant = idx("Assay", "Quantitative measurement")
        i_host = idx("Host", "Name")
        i_pmid = idx("Reference", "PMID")

        out.write(SCHEMA_TSV + "\n")
        n_written = 0
        n_labeled = 0
        for line in raw:
            if len(line) <= i_method:
                continue
            method = line[i_method].strip().lower()
            if not method:
                continue
            if any(k in method for k in METHOD_DROP):
                continue
            if not any(k in method for k in METHOD_KEEP):
                continue
            host = line[i_host].strip().lower()
            if host and "homo sapiens" not in host and "human" not in host:
                continue
            peptide = line[i_epitope].strip()
            allele_raw = line[i_allele].strip()
            if not peptide or not allele_raw:
                continue
            allele = allele_raw.split(",")[0].strip()

            qualitative = line[i_qual].strip().lower()
            label: int | None = None
            if qualitative.startswith("positive"):
                label = 1
            elif qualitative.startswith("negative"):
                label = 0
            score = parse_ic50(line[i_quant].strip())
            if label is None and score is None:
                continue
            key = (peptide, allele)
            if key in seen:
                continue
            seen.add(key)

            rec = InteractionRecord(
                source=DataSource.IEDB,
                peptide=peptide,
                mhc_allele=allele,
                label=label,
                score=score,
                reference_id=(line[i_pmid].strip() or None),
                organism="HomoSapiens",
                method=f"iedb:{method}" if method else "iedb",
            )
            out.write("\t".join(rec.to_row()) + "\n")
            n_written += 1
            n_labeled += 1 if label is not None else 0
            if max_rows is not None and n_written >= max_rows:
                break

    stats.n_records = n_written
    stats.n_unique_peptides = len({k[0] for k in seen})
    stats.n_alleles = len({k[1] for k in seen if k[1]})
    stats.n_labeled = n_labeled
    stats.n_unlabeled = n_written - n_labeled
    if n_written == 0:
        stats.warnings.append("no records ingested — check column layout")
    return stats


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None)
    args = ap.parse_args()
    if not args.in_path.exists():
        print(f"input not found: {args.in_path}", file=sys.stderr)
        return 1
    stats = ingest(args.in_path, args.out, args.max_rows)
    print(stats.summary())
    return 0


if __name__ == "__main__":
    sys.exit(main())
