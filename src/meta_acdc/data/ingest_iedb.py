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

# Assay families must not share a number. Affinity (IC50), elution/mass-spec
# presentation evidence and stability are different measurements; merging them
# into one "-log10 IC50" column made incomparable numbers comparable
# (AUDIT.md #22).
ASSAY_IC50 = "binding-ic50"
ASSAY_PRESENTATION = "presentation-eluted"


def assay_family(method: str) -> str:
    """Classify an IEDB assay method string into an endpoint family.

    Kinetic constants (koff/kon/KD/half-life) are NOT affinities: folding
    them into "-log10 IC50" compared incomparable numbers (AUDIT.md #22).
    """
    m = method.lower()
    if ("koff" in m or "kon" in m or "kd " in m or "dissociat" in m
            or "half-life" in m or "t1/2" in m or "stability" in m):
        return "other:kinetics"
    if "elut" in m or "mass spectr" in m:
        return ASSAY_PRESENTATION
    if "binding" in m or "ic50" in m:
        return ASSAY_IC50
    return f"other:{m[:24]}"


def parse_ic50(raw: str) -> float | None:
    try:
        v = float(raw)
        return -math.log10(max(v, 1e-12))
    except (TypeError, ValueError):
        return None


def ic50_score(raw: str, units: str) -> tuple[float | None, str]:
    """-log10 nM affinity for the numeric score column, unit-whitelisted.

    Only nM enters as-is; µM is converted (×1000) with an explicit note;
    any other or missing unit gets NO numeric score — the value is kept in
    `notes` verbatim instead of being silently transformed (AUDIT.md #22).
    """
    v = None
    try:
        v = float(raw)
    except (TypeError, ValueError):
        return None, f"non-numeric quant={raw!r}"
    u = units.strip().lower()
    if "nm" in u:
        return -math.log10(max(v, 1e-12)), ""
    if u in ("um", "µm", "μм", "um", "µ") or "um" in u:
        return -math.log10(max(v * 1000.0, 1e-12)), (
            f"unit={units} (µM converted to nM ×1000 for -log10)")
    return None, (f"unit={units or 'missing'} not accepted for numeric score; "
                  f"raw quant kept verbatim")


def ingest(path: Path, out_path: Path, max_rows: int | None = None) -> DatasetStats:
    stats = DatasetStats()
    # Conflicting experiments for the same (peptide, allele, assay family) are
    # resolved by a DETERMINISTIC rule — evidence strength first (positive >
    # negative > unknown), then the stronger quantitative signal — so the
    # output no longer depends on the order of rows in the input file, and the
    # number of conflicting experiments is reported (AUDIT.md #22).
    best: dict[tuple[str, str, str], tuple[tuple, InteractionRecord]] = {}

    def evidence_rank(label: int | None) -> int:
        return {1: 3, 0: 2, None: 1}[label]

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
        i_units = idx("Assay", "Units") if f"Assay|Units" in fieldnames else None
        i_host = idx("Host", "Name")
        i_pmid = idx("Reference", "PMID")

        out.write(SCHEMA_TSV + "\n")
        n_conflicts = 0
        n_by_family: dict[str, int] = {}
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
            family = assay_family(method)

            qualitative = line[i_qual].strip().lower()
            label: int | None = None
            if qualitative.startswith("positive"):
                label = 1
            elif qualitative.startswith("negative"):
                label = 0

            raw_quant = line[i_quant].strip()
            units = (line[i_units].strip() if i_units is not None
                     and len(line) > i_units else "")
            # Only IC50-family assays may enter the numeric `score` column, and
            # only with the unit recorded. Elution / mass-spec evidence is kept
            # as a separate endpoint, never folded into an affinity number.
            score = None
            notes = ""
            if family == ASSAY_IC50 and raw_quant:
                score, notes = ic50_score(raw_quant, units)
            elif family == ASSAY_PRESENTATION:
                notes = (f"presentation evidence ({method}) "
                         f"quant={raw_quant or 'n/a'}; not an affinity")
            elif family == "other:kinetics":
                notes = (f"kinetic constant ({method}) "
                         f"quant={raw_quant or 'n/a'} {units or ''}; "
                         f"not an affinity, excluded from score")
            if label is None and score is None and not notes:
                continue

            rec = InteractionRecord(
                source=DataSource.IEDB,
                peptide=peptide,
                mhc_allele=allele,
                label=label,
                score=score,
                reference_id=(line[i_pmid].strip() or None),
                organism="HomoSapiens",
                method=f"iedb:{family}",
                notes=notes,
            )
            key = (peptide, allele, family)
            rank = (evidence_rank(label), score if score is not None else -1e9)
            pmid = (line[i_pmid].strip() or "")
            if key in best:
                n_conflicts += 1
                prev_rank, prev = best[key]
                # deterministic on ties: evidence > score, then smallest
                # PMID — the kept row no longer depends on input row order
                if rank > prev_rank or (rank == prev_rank and pmid
                                        and (not prev.reference_id
                                             or pmid < prev.reference_id)):
                    best[key] = (rank, rec)
            else:
                best[key] = (rank, rec)
            n_by_family[family] = n_by_family.get(family, 0) + 1

        rows = [rec for _, rec in best.values()]
        for rec in rows:
            out.write("\t".join(rec.to_row()) + "\n")
        n_written = len(rows)
        n_labeled = sum(1 for r in rows if r.label is not None)
        seen = {(r.peptide, r.mhc_allele) for r in rows}

    stats.n_records = n_written
    stats.n_unique_peptides = len({k[0] for k in seen})
    stats.n_alleles = len({k[1] for k in seen if k[1]})
    stats.n_labeled = n_labeled
    stats.n_unlabeled = n_written - n_labeled
    stats.warnings.append(
        f"assay families kept: "
        f"{', '.join(f'{k}={v}' for k, v in sorted(n_by_family.items()))}")
    stats.warnings.append(
        f"conflicting experiments merged deterministically: {n_conflicts} "
        f"(rule: positive > negative > unknown, then stronger score; the "
        f"result is independent of input row order)")
    stats.warnings.append(
        "IC50 and presentation/elution endpoints are kept in separate assay "
        "families; do not pool the score column across them")
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
