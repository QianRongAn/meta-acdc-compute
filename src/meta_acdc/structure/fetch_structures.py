"""Fetch and auto-classify TCR-pMHC complexes from STCRDab + RCSB (KN-5 prep).

Pipeline:
1. Scrape the STCRDab browser page for curated TCR structure PDB IDs.
2. Download each PDB from RCSB (small files).
3. Auto-classify TCR-pMHC complexes by chain composition:
   - TCR side: chains A/B (standard TCR convention)
   - pMHC side: one chain > 100 residues (MHC heavy chain) AND one chain with
     5-20 residues (peptide)
4. Write the curated complex list to structures/complexes.tsv

Usage:
    .venv/bin/python src/meta_acdc/structure/fetch_structures.py \
        --out data/raw/structures --workers 8
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

STCRDAB_URL = "https://opig.stats.ox.ac.uk/webapps/stcrdab/Browser?all=true"
RCSB_URL = "https://files.rcsb.org/download/{pdb}.pdb"


def scrape_ids(url: str) -> list[str]:
    req = urllib.request.Request(url, headers={"User-Agent": "meta-acdc/0.1"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        html = resp.read().decode("utf-8", "replace")
    ids = re.findall(r">([0-9][a-z0-9]{3})<", html)
    return sorted(set(ids))


def download_pdb(pdb: str, out_dir: Path) -> Path | None:
    dest = out_dir / f"{pdb.lower()}.pdb"
    if dest.exists() and dest.stat().st_size > 10_000:
        return dest
    try:
        req = urllib.request.Request(
            RCSB_URL.format(pdb=pdb.upper()), headers={"User-Agent": "meta-acdc/0.1"}
        )
        with urllib.request.urlopen(req, timeout=60) as resp, open(dest, "wb") as fh:
            fh.write(resp.read())
        return dest if dest.stat().st_size > 10_000 else None
    except Exception:
        return None


def chain_sizes(pdb_path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with open(pdb_path, errors="replace") as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                chain = line[21].strip()
                sizes[chain] = sizes.get(chain, 0) + 1
    return sizes


def is_tcr_pmhc(pdb_path: Path) -> tuple[bool, dict[str, int]]:
    """Classify: TCR chains A/B + pMHC side with MHC-like (>100) + peptide (5-20)."""
    sizes = chain_sizes(pdb_path)
    tcr_chains = [c for c in ("A", "B") if c in sizes]
    other = {c: n for c, n in sizes.items() if c not in ("A", "B")}
    has_mhc = any(n > 100 for n in other.values())
    has_peptide = any(5 <= n <= 20 for n in other.values())
    return len(tcr_chains) >= 1 and has_mhc and has_peptide, sizes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("data/raw/structures"))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max", type=int, default=0, help="limit downloads (0=all)")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    ids = scrape_ids(STCRDAB_URL)
    print(f"STCRDab entries: {len(ids)}")
    if args.max:
        ids = ids[: args.max]

    print(f"downloading {len(ids)} PDBs with {args.workers} workers ...")
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(lambda p: download_pdb(p, args.out), ids))

    complexes = []
    n_ok = 0
    for pdb, path in zip(ids, results):
        if path is None:
            continue
        n_ok += 1
        ok, sizes = is_tcr_pmhc(path)
        if ok:
            complexes.append((pdb, sizes))

    print(f"downloaded: {n_ok}/{len(ids)}; TCR-pMHC complexes: {len(complexes)}")
    with open(args.out / "complexes.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pdb", "chains"])
        for pdb, sizes in complexes:
            w.writerow([pdb, ";".join(f"{c}:{n}" for c, n in sorted(sizes.items()))])
    return 0


if __name__ == "__main__":
    sys.exit(main())
