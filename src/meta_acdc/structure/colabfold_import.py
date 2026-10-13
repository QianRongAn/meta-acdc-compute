"""Import ColabFold (AF2-multimer) TCR-pMHC predictions into the KN-5 pipeline.

ColabFold writes PDBs; the scoring pipeline (score_predictions.py) reads
mmCIFs via structure/cif.parse_cif. This script converts each ColabFold PDB
into the same minimal _atom_site mmCIF layout that AF3 CIFs use, so the
downstream stack (interface graphs, EGNN scoring, ensemble) runs unchanged.

- chain IDs / residue numbering are preserved as written by ColabFold
  (5-chain FASTA order: alpha / beta / peptide / mhc / b2m -> chains A-E).
  Chain ROLES are re-derived by content in classify_chains(), so IDs are
  irrelevant for the TCR/pMHC split.
- B-factors (ColabFold pLDDT) land in _atom_site.B_iso_or_equiv, the same
  slot parse_cif maps to pLDDT.
- one ColabFold model per job -> output file {job}_model_1.cif, so the
  per-job aggregation in score_predictions.py (rsplit("_model_", 1)) works
  and n_models=1.

Usage:
    .venv/bin/python src/meta_acdc/structure/colabfold_import.py \
        --src data/raw/colabfold_run1 \
        --out data/raw/colabfold_predictions \
        --report data/processed/colabfold_run1_import_report.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from meta_acdc.structure.graph import parse_pdb, classify_chains, Residue

# same column names parse_cif() resolves (see structure/cif.py)
CIF_HEADER = """data_colabfold
#
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.auth_asym_id
_atom_site.auth_comp_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
"""


def job_id_of(path: Path) -> str:
    """5c07_WGPDPAAA_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb
    -> 5c07_WGPDPAAA"""
    stem = path.stem
    for marker in ("_unrelaxed_", "_relaxed_"):
        if marker in stem:
            return stem.split(marker)[0]
    raise ValueError(f"{path.name}: not a ColabFold rank file name")


def write_cif(residues: list[Residue], out_path: Path) -> None:
    """Minimal mmCIF _atom_site loop compatible with structure/cif.parse_cif."""
    lines = [CIF_HEADER]
    for r in residues:
        b = r.plddt if r.plddt is not None else 0.0
        b = f"{b:.2f}"
        for name, x, y, z in r.atoms:
            lines.append(
                f"ATOM {name:<4} {r.chain} {r.resname:<3} "
                f"{r.resid} {x:.3f} {y:.3f} {z:.3f} {b} 1"
            )
    out_path.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, required=True,
                    help="directory of ColabFold PDB outputs (run1/run2 dumps)")
    ap.add_argument("--out", type=Path,
                    default=Path("data/raw/colabfold_predictions"))
    ap.add_argument("--report", type=Path,
                    default=Path("data/processed/colabfold_import_report.tsv"))
    args = ap.parse_args()

    pdb_files = sorted(args.src.glob("*.pdb"))
    if not pdb_files:
        print(f"{args.src}: no .pdb files", file=sys.stderr)
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    report_rows = []
    n_ok = n_skip = n_fail = 0
    for pdb in pdb_files:
        try:
            job = job_id_of(pdb)
        except ValueError as e:
            n_skip += 1
            report_rows.append({"job_id": pdb.stem, "n_residues": "",
                                "chains": "", "status": f"SKIP: {e}"})
            continue
        out_cif = args.out / f"{job}_model_1.cif"
        try:
            residues = parse_pdb(pdb)
            if len(residues) < 10:
                raise ValueError(f"too few residues ({len(residues)})")
            roles = classify_chains(residues)
            per_chain: dict[str, int] = {}
            for r in residues:
                per_chain[r.chain] = per_chain.get(r.chain, 0) + 1
            chain_summary = " ".join(
                f"{c}:{per_chain[c]}:{roles.get(c, '?')}"
                for c in sorted(per_chain))
            write_cif(residues, out_cif)
            n_ok += 1
            report_rows.append({"job_id": job, "n_residues": len(residues),
                                "chains": chain_summary, "status": "OK"})
            print(f"{job}: {len(residues)} residues [{chain_summary}]",
                  flush=True)
        except Exception as e:  # noqa: BLE001 — per-file QC loop
            n_fail += 1
            report_rows.append({"job_id": job, "n_residues": "",
                                "chains": "", "status": f"FAIL: {e}"})
            print(f"{job}: FAIL {e}", file=sys.stderr, flush=True)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["job_id", "n_residues", "chains",
                                           "status"], delimiter="\t")
        w.writeheader()
        w.writerows(report_rows)
    print(f"imported {n_ok} OK / {n_skip} skipped / {n_fail} failed "
          f"-> {args.report}")
    return 0 if n_fail == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
