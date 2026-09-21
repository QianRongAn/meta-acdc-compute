"""Batch TCRmodel2 submissions for domain-adaptation data (423 complexes).

Why: the EGNN ranking on AF3-predicted structures is not reproducible across
model instances (Rashomon effect). Fix: include TCRmodel2-style predicted
structures in TRAINING. TCRmodel2 (tcrmodel.ibbr.umd.edu) has no daily
quota and a scriptable form (POST /tcrpmhc_submitjob1/1); ~15 min per
complex. We submit the 423 training complexes' native sequences.

Politeness: sequential submissions with delay; jobs tracked in a TSV.

Usage:
    .venv/bin/python src/meta_acdc/structure/tcrmodel_batch.py \
        --struct-dir data/raw/structures \
        --state data/processed/tcrmodel_jobs.tsv \
        --limit 50 --delay 10
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

from meta_acdc.structure.graph import AA3TO1, classify_chains, parse_pdb

SUBMIT_URL = "https://tcrmodel.ibbr.umd.edu/tcrpmhc_submitjob1/1"
ALLELE_BY_MHC_PREFIX = {
    "GSHSM": "HLA-A*02:01", "GSHSL": "HLA-C*07:02",
    "GPHSLRYFVTAVSRP": "HLA-E*01:03",
}


def extract_components(pdb_path: Path):
    res = parse_pdb(pdb_path)
    roles = classify_chains(res)
    seqs = {}
    for c, r in roles.items():
        seq = "".join(AA3TO1.get(x.resname[:3].upper(), "X")
                      for x in res if x.chain == c)
        seqs.setdefault(r, []).append(seq)
    tcr = seqs.get("tcr", [])
    alpha = next((s for s in tcr if any(m in s for m in
                 ("YFCAV", "YFCAI", "YFCAL", "YLCAV", "YLCAI", "YLCAL"))), None)
    beta = next((s for s in tcr if any(m in s for m in ("YFCAS", "YLCAS"))), None)
    pep = (seqs.get("peptide") or [None])[0]
    mhc = (seqs.get("mhc") or [None])[0]
    if not all([alpha, beta, pep, mhc]):
        return None
    allele = next((a for k, a in ALLELE_BY_MHC_PREFIX.items()
                   if mhc.startswith(k)), "HLA-A*02:01")
    return {"alpha": alpha, "beta": beta, "pep": pep, "mhc": mhc,
            "allele": allele}


def submit(comp: dict) -> str:
    fields = {
        "alphachain": comp["alpha"], "betachain": comp["beta"],
        "pepchain": comp["pep"], "mhc1speciestype": "human",
        "mhc1a": comp["allele"], "mhc1aseq": comp["mhc"], "ar": "on",
        "pdbblacklist": "None", "pdbincludelist": "None",
        "datecutoff": "2019-07-02",
    }
    data = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(SUBMIT_URL, data=data,
                                 headers={"User-Agent": "Mozilla/5.0"})
    r = urllib.request.urlopen(req, timeout=60)
    final = r.geturl()
    m = re.search(r"/rtcr/([A-Za-z0-9_]+)", final)
    if not m:
        raise RuntimeError(f"no job id in {final}")
    return m.group(1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct-dir", type=Path, default=Path("data/raw/structures"))
    ap.add_argument("--state", type=Path,
                    default=Path("data/processed/tcrmodel_jobs.tsv"))
    ap.add_argument("--limit", type=int, default=50)
    ap.add_argument("--delay", type=int, default=10)
    args = ap.parse_args()

    done: set[str] = set()
    if args.state.exists():
        with open(args.state) as fh:
            done = {r["pdb"] for r in csv.DictReader(fh, delimiter="\t")}

    with open(args.struct_dir / "complexes.tsv") as fh:
        complexes = [r["pdb"] for r in csv.DictReader(fh, delimiter="\t")]

    n_submitted = 0
    with open(args.state, "a", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        if not args.state.exists() or args.state.stat().st_size == 0:
            w.writerow(["pdb", "job_id"])
        for pdb in complexes:
            if pdb in done or n_submitted >= args.limit:
                continue
            comp = extract_components(args.struct_dir / f"{pdb.lower()}.pdb")
            if comp is None:
                continue
            try:
                job_id = submit(comp)
            except Exception as e:
                print(f"{pdb}: submit failed ({e})", flush=True)
                time.sleep(args.delay)
                continue
            w.writerow([pdb, job_id])
            out.flush()
            done.add(pdb)
            n_submitted += 1
            print(f"{n_submitted:3d}. {pdb} -> {job_id}", flush=True)
            time.sleep(args.delay)
    print(f"submitted {n_submitted} jobs; state: {args.state}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
