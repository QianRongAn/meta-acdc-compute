"""NetMHCpan-4.2 pre-filter benchmark (paper 2, experiment 1).

Validates the proposal's claim: pre-filtering with NetMHCpan compresses the
proteome-scale peptide search space to ~5M candidates while retaining the
true binders.

Protocol: 50 known HLA-A*02:01 9-mer binders (IEDB positives) + 50 shuffled
decoys -> one NetMHCpan-4.2 submission -> %Rank_EL analysis:
- compression: fraction of candidates passing %rank<2 (the standard weak
  binder threshold)
- recall: fraction of known binders passing %rank<2

Usage:
    .venv/bin/python src/meta_acdc/data/netmhcpan_bench.py \
        --iedb data/processed/iedb.clean.tsv --n 50
"""

from __future__ import annotations

import argparse
import random
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

SERVICE = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
CONFIG = "/var/www/html/services/NetMHCpan-4.2/webface.cf"
ALLELE = "HLA-A*02:01"


def load_binders(path: Path, allele: str, n: int, seed: int = 0) -> list[str]:
    import csv
    rng = random.Random(seed)
    pep9 = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pep = row["peptide"]
            if (row["label"] == "1" and row["mhc_allele"] == allele
                    and len(pep) == 9 and pep.isalpha() and pep.isupper()):
                pep9.append(pep)
    rng.shuffle(pep9)
    return pep9[:n]


def shuffle_peptides(peps: list[str], seed: int = 1) -> list[str]:
    rng = random.Random(seed)
    out = []
    for p in peps:
        s = list(p)
        rng.shuffle(s)
        out.append("".join(s))
    return out


def submit(peptides: list[str], length: int = 9, wait: int = 5) -> str:
    """Submit peptide list; return the jobid (from the redirect URL)."""
    paste = "\n".join(peptides)
    data = urllib.parse.urlencode({
        "configfile": CONFIG,
        "PEPPASTE": paste,
        "allele": ALLELE,
    }).encode()
    req = urllib.request.Request(f"{SERVICE}?wait={wait}", data=data,
                                 headers={"User-Agent": "meta-acdc/0.1"})
    with urllib.request.urlopen(req, timeout=wait + 60) as resp:
        final_url = resp.geturl()
    jobid = urllib.parse.parse_qs(urllib.parse.urlparse(final_url).query).get("jobid")
    if not jobid:
        raise RuntimeError(f"no jobid in redirect URL: {final_url}")
    return jobid[0]


def poll_result(jobid: str, timeout_s: int = 900) -> str:
    """Poll the status page until the job finishes; return the results page."""
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        url = f"{SERVICE}?jobid={jobid}&wait=20"
        with urllib.request.urlopen(urllib.request.Request(
                url, headers={"User-Agent": "meta-acdc/0.1"}), timeout=60) as resp:
            text = resp.read().decode("utf-8", "replace")
        if "Status: finished" in text or ("Pos" in text and "Allele" in text):
            return text
        if "failed" in text.lower():
            raise RuntimeError(f"job failed: {text[:300]}")
        time.sleep(15)
    raise TimeoutError(f"job {jobid} not finished within {timeout_s}s")


def parse_ranks(text: str) -> list[float]:
    """Extract %Rank_EL from the results page/table.

    The results table has a header 'Pos MHC Allele Peptide ... %Rank_EL'
    followed by rows; parse any token that parses as a float after the
    peptide token per line.
    """
    import html as H
    clean = re.sub(r"<[^>]+>", " ", text)
    lines = [ln.split() for ln in clean.splitlines() if ln.strip()]
    ranks = []
    header_seen = False
    for parts in lines:
        if "Pos" in parts and "Allele" in parts:
            header_seen = True
            continue
        if not header_seen or len(parts) < 4:
            continue
        # peptide token then numbers; last float is typically %Rank
        floats = []
        for p in parts:
            try:
                floats.append(float(p))
            except ValueError:
                continue
        if floats:
            ranks.append(floats[-1])
    return ranks


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--iedb", type=Path, default=Path("data/processed/iedb.clean.tsv"))
    ap.add_argument("--n", type=int, default=50)
    args = ap.parse_args()

    binders = load_binders(args.iedb, ALLELE, args.n)
    decoys = shuffle_peptides(binders)
    print(f"submitting {len(binders) * 2} peptides (HLA-A*02:01, 9-mers) ...",
          flush=True)
    jobid = submit(binders + decoys)
    print(f"job {jobid} submitted, polling ...", flush=True)
    result_page = poll_result(jobid)
    ranks = parse_ranks(result_page)
    if len(ranks) != len(binders) * 2:
        print(f"warning: parsed {len(ranks)} ranks, expected {len(binders) * 2}")
    r_b = ranks[: len(binders)]
    r_d = ranks[len(binders):]
    for thresh in (0.5, 2.0):
        recall = sum(r < thresh for r in r_b) / len(r_b)
        decoy_pass = sum(r < thresh for r in r_d) / len(r_d)
        print(f"%rank < {thresh}: binder recall {recall:.2f}, "
              f"decoy pass-rate {decoy_pass:.2f}")
    overall = sum(r < 2.0 for r in ranks) / len(ranks)
    print(f"compression at %rank<2: {overall:.2f} of candidates retained")
    return 0


if __name__ == "__main__":
    sys.exit(main())
