"""Normalize IEDB mhc_ligand_full CSV into the unified schema.

The IEDB export is a single CSV with columns including:
    object_type, epitope name, mhc allele names, mhc class, assay group,
    qualitative_measure, quantitative_measurement, host organism name,
    pubMed id, references, ...

We keep MHC ligand/binding assays on human data as presentation-level
records (used for NetMHCpan-style pre-filtering benchmarks). T cell assay
rows are excluded — activation data comes from VDJdb/ACDC.

Labels: qualitative_measure Positive -> 1, Negative -> 0.
Score: -log10(IC50 nM) from quantitative_measurement when parseable.

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

ASSAY_GROUPS = {"binding", "mhc ligand", "ligand", "mhc ligand assay"}


def parse_ic50(raw: str) -> float | None:
    try:
        v = float(raw)
        return -math.log10(max(v, 1e-12))
    except (TypeError, ValueError):
        return None


def ingest(path: Path, out_path: Path, max_rows: int | None = None) -> DatasetStats:
    stats = DatasetStats()
    seen: set[tuple[str, str]] = set()
    rows: list[list[str]] = []

    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        fieldnames = {f.lower(): f for f in (reader.fieldnames or [])}

        def get(line: dict, key: str) -> str:
            return (line.get(fieldnames[key]) or "").strip()

        for line in reader:
            if max_rows is not None and len(rows) >= max_rows:
                break
            assay_group = get(line, "assay group").lower()
            if assay_group and assay_group not in ASSAY_GROUPS:
                continue
            host = get(line, "host organism name").lower()
            if host and "homo sapiens" not in host and "human" not in host:
                continue

            peptide = get(line, "epitope name")
            allele_raw = get(line, "mhc allele names")
            if not peptide or not allele_raw:
                continue
            allele = allele_raw.split(",")[0].strip()

            qualitative = get(line, "qualitative_measure").lower()
            label: int | None = None
            if qualitative.startswith("positive"):
                label = 1
            elif qualitative.startswith("negative"):
                label = 0

            score = parse_ic50(get(line, "quantitative_measurement"))
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
                reference_id=get(line, "pubMed id") or None,
                organism="HomoSapiens",
                method=f"iedb:{assay_group}" if assay_group else "iedb",
            )
            rows.append(rec.to_row())

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        fh.write(SCHEMA_TSV + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")

    stats.n_records = len(rows)
    stats.n_unique_tcr_beta = 0  # IEDB ligand data has no TCR
    stats.n_unique_peptides = len({r[7] for r in rows})
    stats.n_alleles = len({r[8] for r in rows if r[8]})
    stats.n_labeled = sum(1 for r in rows if r[9])
    stats.n_unlabeled = stats.n_records - stats.n_labeled
    if stats.n_records == 0:
        stats.warnings.append("no records ingested — check input columns")
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
