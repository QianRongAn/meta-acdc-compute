"""Download VDJdb (TCR-pMHC specificity database).

Sources tried in order (stdlib-only, no external deps):
1. GitHub releases of antigenomics/vdjdb-db (community mirror of the database)
2. Official site download (vdjdb.cdr3.net)

Usage:
    python download_vdjdb.py [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

GITHUB_API = "https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest"
DEFAULT_OUT = Path(__file__).resolve().parents[3] / "data" / "raw" / "vdjdb"


def fetch_json(url: str, timeout: int = 30) -> dict:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "meta-acdc/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def download(url: str, dest: Path, timeout: int = 120) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {url} -> {dest}")
    req = urllib.request.Request(url, headers={"User-Agent": "meta-acdc/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    size = dest.stat().st_size
    print(f"done: {dest.name} ({size / 1e6:.1f} MB)")
    if size < 1024:
        raise RuntimeError(f"suspiciously small file ({size} B) — download likely failed")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    ok = False
    # 1) GitHub release assets of the community DB mirror
    try:
        rel = fetch_json(GITHUB_API)
        for asset in rel.get("assets", []):
            name = asset.get("name", "")
            if name.endswith((".zip", ".txt", ".tsv", ".csv", ".gz")):
                download(asset["browser_download_url"], args.out / name)
                ok = True
    except Exception as e:  # noqa: BLE001
        print(f"GitHub mirror failed: {e}", file=sys.stderr)

    # 2) Official site (fallback, may change)
    if not ok:
        candidates = [
            "https://vdjdb.cdr3.net/search/downloads",
        ]
        for url in candidates:
            try:
                download(url, args.out / "vdjdb_download.html")
                ok = True
                break
            except Exception as e:  # noqa: BLE001
                print(f"{url} failed: {e}", file=sys.stderr)
    if not ok:
        print("All VDJdb download sources failed.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
