"""Poll TCRmodel2 jobs and download finished results.

Reads the submission state TSV (pdb, job_id, downloaded), polls each job
page until result links appear, downloads the modeled PDB to
data/raw/tcrmodel/{pdb}.pdb (TCRmodel2 output = PDB format, which feeds the
existing PDB pipeline directly).

Usage:
    .venv/bin/python src/meta_acdc/structure/tcrmodel_poll.py \
        --state data/processed/tcrmodel_jobs.tsv \
        --out data/raw/tcrmodel \
        --rounds 60 --interval 60
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import time
import urllib.request
from pathlib import Path

JOB_URL = "https://tcrmodel.ibbr.umd.edu/rtcr/{job_id}"


def get_page(job_id: str, timeout: int = 45) -> str:
    req = urllib.request.Request(JOB_URL.format(job_id=job_id),
                                 headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def find_download_links(html: str) -> list[str]:
    return re.findall(r'href="([^"]*(?:download|/static/result)[^"]*)"',
                      html, re.I)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", type=Path,
                    default=Path("data/processed/tcrmodel_jobs.tsv"))
    ap.add_argument("--out", type=Path, default=Path("data/raw/tcrmodel"))
    ap.add_argument("--rounds", type=int, default=60)
    ap.add_argument("--interval", type=int, default=60)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    for rnd in range(args.rounds):
        rows = []
        with open(args.state, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        pending = [r for r in rows
                   if not r.get("downloaded") and not r.get("failed")]
        done_this_round = 0
        for r in pending:
            pdb, job_id = r["pdb"], r["job_id"]
            try:
                html = get_page(job_id)
            except Exception as e:
                print(f"{pdb}: poll error {e}", flush=True)
                continue
            links = find_download_links(html)
            if links:
                # download first PDB-ish link
                url = links[0] if links[0].startswith("http") \
                    else f"https://tcrmodel.ibbr.umd.edu{links[0]}"
                try:
                    req = urllib.request.Request(
                        url, headers={"User-Agent": "Mozilla/5.0"})
                    with urllib.request.urlopen(req, timeout=120) as resp:
                        content = resp.read()
                    dest = args.out / f"{pdb}.pdb"
                    dest.write_bytes(content)
                    r["downloaded"] = "1"
                    print(f"{pdb}: DOWNLOADED ({len(content)} bytes)", flush=True)
                    done_this_round += 1
                except Exception as e:
                    print(f"{pdb}: download error {e}", flush=True)
            elif "failed" in html.lower() or "error" in html.lower():
                r["failed"] = "1"
                print(f"{pdb}: job failed", flush=True)
        if done_this_round:
            # rewrite state with updated flags
            with open(args.state, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()),
                                   delimiter="\t")
                w.writeheader()
                for r in rows:
                    w.writerow(r)
        n_pending = len([r for r in rows
                         if not r.get("downloaded") and not r.get("failed")])
        if n_pending == 0:
            print("all jobs resolved")
            break
        print(f"round {rnd + 1}: {done_this_round} downloaded, "
              f"{n_pending} pending", flush=True)
        if rnd < args.rounds - 1:
            time.sleep(args.interval)
    return 0


if __name__ == "__main__":
    sys.exit(main())
