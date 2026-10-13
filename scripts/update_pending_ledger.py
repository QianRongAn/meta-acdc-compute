#!/usr/bin/env python3
"""Update af3_pending_upload/already_uploaded.tsv with a batch of new jobs.

- reads --jobs (one job_id per line)
- for job_ids already present, appends --append-note to their status
- for new job_ids, adds rows "job_id\\t<--status>"
- rewrites the ledger with LF endings (fixes mixed CRLF from the
  Windows-side tool while keeping content identical otherwise)

Usage:
    .venv/bin/python scripts/update_pending_ledger.py \
        --jobs /tmp/run1_jobs.txt \
        --status "results on disk (colabfold run1, 2026-10-13)" \
        --append-note " [+colabfold run1 2026-10-13]"
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=Path, required=True)
    ap.add_argument("--ledger", type=Path,
                    default=Path("af3_pending_upload/already_uploaded.tsv"))
    ap.add_argument("--status", required=True,
                    help="status for newly added rows")
    ap.add_argument("--append-note", default="",
                    help="note appended to the status of already-present rows")
    args = ap.parse_args()

    jobs = [ln.strip().rstrip("\r") for ln in args.jobs.read_text().splitlines()
            if ln.strip()]
    ledger = args.ledger.read_text(encoding="utf-8", errors="replace")
    lines = ledger.splitlines()  # drops \r and \n uniformly -> LF on rewrite
    if not lines:
        print(f"{args.ledger}: empty", file=sys.stderr)
        return 1

    existing: dict[str, str] = {}
    for ln in lines[1:]:
        if "\t" in ln:
            jid, status = ln.split("\t", 1)
            existing[jid] = status

    n_new = n_note = 0
    for jid in jobs:
        if jid in existing:
            if args.append_note and args.append_note not in existing[jid]:
                existing[jid] += args.append_note
                n_note += 1
        else:
            existing[jid] = args.status
            n_new += 1

    out = [lines[0]] + [f"{jid}\t{st}" for jid, st in sorted(existing.items())]
    args.ledger.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"ledger {args.ledger}: {n_new} new rows, {n_note} notes appended; "
          f"total rows {len(existing)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
