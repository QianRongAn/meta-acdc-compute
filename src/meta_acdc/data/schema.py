"""Unified schema for TCR-pMHC interaction data (Year 1 Q1, data pipeline v0).

All ingested datasets (VDJdb, IEDB, 10x, in-house ACDC) are normalized to
`InteractionRecord` and stored as TSV + SQLite. Schema decisions:
- TCR identity: CDR3 beta (and alpha when available) + V/J gene. Sequence-based
  identity only at v0; structure comes later via AlphaFold in `structure/`.
- Peptide identity: sequence (and MHC allele context for presentation).
- Label: binding signal. We distinguish *presentation* (IEDB) from *activation*
  (VDJdb / ACDC) — the latter is the target variable of the model.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class DataSource(str, Enum):
    VDJDB = "vdjdb"
    IEDB = "iedb"
    TENX = "10x"
    ACDC = "acdc"  # in-house ACDC platform
    MCPAS = "mcpas"


@dataclass
class InteractionRecord:
    """One TCR-pMHC interaction datapoint in the unified schema."""

    source: DataSource
    # TCR
    cdr3_alpha: Optional[str] = None
    cdr3_beta: Optional[str] = None
    tra_v: Optional[str] = None
    tra_j: Optional[str] = None
    trb_v: Optional[str] = None
    trb_j: Optional[str] = None
    # peptide / MHC
    peptide: Optional[str] = None
    mhc_allele: Optional[str] = None  # e.g. "HLA-A*02:01"
    # labels
    label: Optional[float] = None  # 1/0 binding; None = unlabeled
    score: Optional[float] = None  # continuous signal (enrichment, intensity)
    # provenance
    reference_id: Optional[str] = None  # original DB id / paper PMID
    organism: Optional[str] = None
    method: Optional[str] = None  # experimental assay
    notes: str = field(default="", repr=False)

    def to_row(self) -> list[str]:
        return [
            self.source.value,
            self.cdr3_alpha or "",
            self.cdr3_beta or "",
            self.tra_v or "",
            self.tra_j or "",
            self.trb_v or "",
            self.trb_j or "",
            self.peptide or "",
            self.mhc_allele or "",
            "" if self.label is None else str(int(self.label)),
            "" if self.score is None else f"{self.score:.6g}",
            self.reference_id or "",
            self.organism or "",
            self.method or "",
            self.notes,
        ]


SCHEMA_HEADER = [
    "source",
    "cdr3_alpha",
    "cdr3_beta",
    "tra_v",
    "tra_j",
    "trb_v",
    "trb_j",
    "peptide",
    "mhc_allele",
    "label",
    "score",
    "reference_id",
    "organism",
    "method",
    "notes",
]

SCHEMA_TSV = "\t".join(SCHEMA_HEADER)


@dataclass
class DatasetStats:
    """Sanity-check statistics after ingestion."""

    n_records: int = 0
    n_unique_tcr_beta: int = 0
    n_unique_peptides: int = 0
    n_alleles: int = 0
    n_labeled: int = 0
    n_unlabeled: int = 0
    warnings: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"records: {self.n_records}",
            f"unique CDR3-beta: {self.n_unique_tcr_beta}",
            f"unique peptides: {self.n_unique_peptides}",
            f"MHC alleles: {self.n_alleles}",
            f"labeled: {self.n_labeled}, unlabeled: {self.n_unlabeled}",
        ]
        if self.warnings:
            lines.append("warnings:")
            lines.extend(f"  - {w}" for w in self.warnings)
        return "\n".join(lines)
