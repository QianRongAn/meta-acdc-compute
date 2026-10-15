"""Normalize VDJdb TSV into the unified InteractionRecord schema.

VDJdb 2026 file format (vdjdb.txt): TSV with header
    complex.id  gene  cdr3  v.segm  j.segm  species  mhc.a  mhc.b  mhc.class
    antigen.epitope  antigen.gene  antigen.species  reference.id  vdjdb.score
    TCR_hash  method  meta  cdr3fix  ...

We import rows with CDR3 + antigen.epitope present and keep vdjdb.score as
**evidence confidence** (0 = no evidence .. 3 = strong evidence). It is NOT a
binding label. The `score` column holds that confidence value unchanged, and
a duplicate `confidence` column is written so that no downstream reader has to
assume `score` means affinity. The `label` column is deliberately always
empty (AUDIT.md #1). Deduplication key: (cdr3, gene, epitope, mhc.a), keeping
the max-confidence row.

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
    # `gene` says which chain the CDR3 belongs to. It used to be discarded and
    # every CDR3 written to cdr3_beta, which merged alpha and beta chains and
    # let an alpha-only match count as the same receptor (AUDIT.md #7).
    "gene": "gene",
    "cdr3": "cdr3",
    "cdr3.alpha": "cdr3",
    "cdr3.beta": "cdr3",
    # v.segm / j.segm land in a staging field and are routed to tra_* or
    # trb_* by `gene` below. Writing them straight into trb_v/trb_j labelled
    # an alpha-chain row as beta (AUDIT.md #7).
    "v.segm": "v_segm",
    "j.segm": "j_segm",
    "species": "organism",
    "mhc.a": "mhc_allele",
    "mhc.b": None,
    "mhc.class": None,
    "antigen.epitope": "peptide",
    "antigen.gene": None,
    "antigen.species": None,
    "reference.id": "reference_id",
    "vdjdb.score": "score",  # evidence confidence, NOT a binding label
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

            # Route the CDR3 to its own chain column. Rows whose `gene` is
            # neither TRA nor TRB are dropped rather than silently filed as
            # beta.
            gene = (rec.gene or "").upper()
            if gene in ("TRA", "A"):
                rec.cdr3_alpha, rec.cdr3_beta = rec.cdr3, None
                rec.tra_v, rec.tra_j = rec.v_segm, rec.j_segm
                rec.trb_v = rec.trb_j = None
            elif gene in ("TRB", "B"):
                rec.cdr3_beta, rec.cdr3_alpha = rec.cdr3, None
                rec.trb_v, rec.trb_j = rec.v_segm, rec.j_segm
                rec.tra_v = rec.tra_j = None
            else:
                stats.n_skipped_no_gene += 1
                continue

            if (rec.cdr3_beta is None and rec.cdr3_alpha is None) \
                    or rec.peptide is None:
                continue
            # vdjdb.score is the database's EVIDENCE CONFIDENCE, not a binding
            # measurement: 0 means "low confidence / no information", not
            # "does not bind". The old importer turned score>0 into a positive
            # label and score==0 into a negative, which manufactured 185,117
            # "non-binders" (92.9% of the set) and every number derived from
            # them (AUDIT.md #1). The label column is now left empty on
            # purpose; the value is written twice, as `score` (for the
            # existing consumers) and as `confidence` (so the meaning is
            # explicit at the point of use).
            rec.score = parse_score(rec.score)
            rec.confidence = rec.score
            rec.v_segm = rec.j_segm = None
            rec.label = None
            chain = rec.cdr3_beta or rec.cdr3_alpha or ""
            key = (chain, gene, rec.peptide, rec.mhc_allele or "")
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
    stats.n_unique_tcr_beta = len({r[2] for r in rows if r[2]})
    stats.n_unique_tcr_alpha = len({r[1] for r in rows if r[1]})
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
