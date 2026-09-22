"""Re-label AF3 instance CIFs by matching their TCR chains to crystal structures.

Why: kn5_submission_list.tsv has its tcr_a/tcr_b columns swapped relative to
the AF3 job FASTA (tcr_a holds the beta chain), so the original import could
never match the alpha chain and silently fell back to "first candidate row"
for the peptide. For peptides shared by several TCR families (the A6/B7, 2C,
JM22, MART-1 sets) this mislabeled the structures (e.g. a 1qrn/B7 structure
was written as 1ao7/B7).

This script recovers the TRUE identity from the coordinates: it parses each
CIF's TCR alpha/beta sequences, matches them against the 291 crystal
structures' TCR signatures, and renames the CIF/confidence files to the
correct {pdb}_{peptide}. When a signature is shared by several crystals, it
prefers one that is a KN-5 candidate for that peptide.

Usage:
    .venv/bin/python src/meta_acdc/structure/relabel_instances_by_chains.py \
        --instances data/raw/af3_instances \
        --legacy data/raw/af3_predictions \
        --manifest data/processed/af3_relabel_map.tsv
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
from collections import defaultdict
from pathlib import Path

from meta_acdc.structure.cif import parse_cif
from meta_acdc.structure.graph import TCR, classify_chains, parse_pdb
from meta_acdc.structure.vdjdb_structure_map import aa3to1

STRUCT = Path("data/raw/structures")


def tcr_signature(residues) -> tuple[str, ...]:
    roles = classify_chains(residues)
    tcr = [c for c in {r.chain for r in residues} if roles.get(c) == TCR]
    out = []
    for c in tcr:
        rs = sorted([r for r in residues if r.chain == c], key=lambda x: x.resid)
        out.append("".join(aa3to1(r.resname) for r in rs))
    return tuple(sorted(out))


def peptide_of(residues) -> str:
    from meta_acdc.structure.graph import PEPTIDE
    roles = classify_chains(residues)
    peps = [c for c in {r.chain for r in residues} if roles.get(c) == PEPTIDE]
    if not peps:
        return ""
    c = min(peps, key=lambda ch: sum(1 for r in residues if r.chain == ch))
    rs = sorted([r for r in residues if r.chain == c], key=lambda x: x.resid)
    return "".join(aa3to1(r.resname) for r in rs)


def build_crystal_index() -> dict[tuple, list[str]]:
    idx: dict[tuple, list[str]] = defaultdict(list)
    for p in sorted(STRUCT.glob("*.pdb")):
        try:
            s = tcr_signature(parse_pdb(p))
        except Exception:  # noqa: BLE001
            continue
        if len(s) == 2:
            idx[s].append(p.stem)
    return idx


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--instances", type=Path, default=Path("data/raw/af3_instances"))
    ap.add_argument("--legacy", type=Path, default=Path("data/raw/af3_predictions"))
    ap.add_argument("--list", type=Path,
                    default=Path("data/processed/kn5_submission_list.tsv"))
    ap.add_argument("--manifest", type=Path,
                    default=Path("data/processed/af3_relabel_map.tsv"))
    ap.add_argument("--apply", action="store_true",
                    help="actually rename files (default: dry run)")
    args = ap.parse_args()

    crystal = build_crystal_index()
    cand_pdbs: dict[str, set[str]] = defaultdict(set)
    with open(args.list, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cand_pdbs[row["peptide"]].add(row["pdb"])

    mapping = []

    def relabel_dir(d: Path, legacy: bool) -> None:
        # group files by current job id
        jobs = sorted({os.path.basename(f).split("_model_")[0].split("_conf")[0]
                       for f in glob.glob(str(d / "*_model_0.cif"))})
        for old in jobs:
            f0 = d / f"{old}_model_0.cif"
            try:
                res = parse_cif(f0)
            except Exception:  # noqa: BLE001
                continue
            sig = tcr_signature(res)
            pep = peptide_of(res)
            options = crystal.get(sig, [])
            prefer = [p for p in options if p in cand_pdbs.get(pep, set())]
            true = sorted(prefer or options)
            new = f"{true[0]}_{pep}" if true else old
            status = "ok" if new == old else ("unknown" if not true else "changed")
            mapping.append({"instance": d.name, "old_job": old, "new_job": new,
                            "peptide": pep, "true_pdbs": ";".join(true),
                            "status": status})
            if not args.apply or new == old:
                continue
            for src in sorted(d.glob(f"{old}_model_*.cif")):
                tag = src.name.split("_model_")[1]
                dst = d / f"{new}_model_{tag}"
                if dst.exists():
                    print(f"  {d.name}: {src.name} -> {dst.name} EXISTS, skip",
                          flush=True)
                    continue
                src.rename(dst)
                if legacy:
                    cs = d / f"{old}_conf.json"
                    if cs.exists():
                        cs.rename(d / f"{new}_conf.json")
                else:
                    for conf in d.glob(f"{old}_conf_*.json"):
                        conf.rename(d / f"{new}_conf_{conf.name.split('_conf_')[1]}")

    for d in sorted(args.instances.iterdir()):
        if d.is_dir():
            relabel_dir(d, legacy=False)
    if args.legacy.is_dir():
        relabel_dir(args.legacy, legacy=True)

    with open(args.manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "instance", "old_job", "new_job", "peptide", "true_pdbs", "status"])
        w.writeheader()
        for m in mapping:
            w.writerow(m)

    n = len(mapping)
    changed = sum(1 for m in mapping if m["status"] == "changed")
    unknown = sum(1 for m in mapping if m["status"] == "unknown")
    print(f"records: {n}  changed: {changed}  unknown: {unknown}  "
          f"({'APPLIED' if args.apply else 'dry run'})")
    print(f"map written to {args.manifest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
