"""Normalize VDJdb TSV into the unified InteractionRecord schema.

VDJdb 2026 file format (vdjdb.txt): TSV with header
    complex.id  gene  cdr3  v.segm  j.segm  species  mhc.a  mhc.b  mhc.class
    antigen.epitope  antigen.gene  antigen.species  reference.id  vdjdb.score
    TCR_hash  method  meta  cdr3fix  ...

We import rows with CDR3 + antigen.epitope present, keep vdjdb.score
(0-3: 0 no response .. 3 strong response) as the activation signal.
Deduplication key: (cdr3, antigen.epitope, mhc.a), keeping the max-score row.

Usage:
    python ingest_vdjdb.py --in vdjdb.txt --out processed/vdjdb.clean.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from schema import DataSource, DatasetStats, InteractionRecord, SCHEMA_TSV

# friendly VDJdb column name -> record field; None = skipped
COLUMN_MAP = {
    "complex.id": "reference_id",
    "gene": None,  # "TRB" — V/J come from v.segm / j.segm
    "cdr3": "cdr3_beta",
    "v.segm": "trb_v",
    "j.segm": "trb_j",
    "species": "organism",
    "mhc.a": "mhc_allele",
    "mhc.b": None,
    "mhc.class": None,
    "antigen.epitope": "peptide",
    "antigen.gene": None,
    "antigen.species": None,
    "reference.id": "reference_id",
    "vdjdb.score": "score",
    "tcr_hash": None,
    "method": "method",
    "meta": None,
    "cdr3fix": None,
    "web.method": None,
    "web.method.seq": None,
    "web.cdr3fix.nc": None,
    "web.cdr3fix.unmp": None,
}


def parse_score(raw: str) -> float | None:
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


def ingest(path: Path, out_path: Path) -> DatasetStats:
    stats = DatasetStats()
    best: dict[tuple[str, str, str], InteractionRecord] = {}

    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = [h.lower() for h in (reader.fieldnames or [])]
        unknown = [h for h in header if h not in COLUMN_MAP]
        if unknown:
            stats.warnings.append(f"unknown columns: {unknown}")

        for line in reader:
            rec = InteractionRecord(source=DataSource.VDJDB)
            for col in header:
                field = COLUMN_MAP.get(col)
                if field is None:
                    continue
                setattr(rec, field, (line.get(col) or "").strip() or None)

            if rec.cdr3_beta is None or rec.peptide is None:
                continue
            # normalize HLA allele style: HLA-A*02:01 -> A*02:01 handled downstream;
            # keep raw for now, provenance preserved
            rec.score = parse_score(rec.score)
            if rec.score is not None:
                rec.label = 1 if rec.score > 0 else 0
            key = (rec.cdr3_beta, rec.peptide, rec.mhc_allele or "")
            prev = best.get(key)
            if prev is None or (rec.score or -1) > (prev.score or -1):
                best[key] = rec

    rows = [r.to_row() for r in best.values()]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        fh.write(SCHEMA_TSV + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")

    stats.n_records = len(rows)
    stats.n_unique_tcr_beta = len({r[2] for r in rows})
    stats.n_unique_peptides = len({r[7] for r in rows})
    stats.n_alleles = len({r[8] for r in rows if r[8]})
    stats.n_labeled = sum(1 for r in rows if r[9])
    stats.n_unlabeled = stats.n_records - stats.n_labeled
    if stats.n_records == 0:
        stats.warnings.append("no records ingested — check the input file format")
    return stats


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if not args.in_path.exists():
        print(f"input not found: {args.in_path}", file=sys.stderr)
        return 1
    stats = ingest(args.in_path, args.out)
    print(stats.summary())
    return 0


if __name__ == "__main__":
    sys.exit(main())
