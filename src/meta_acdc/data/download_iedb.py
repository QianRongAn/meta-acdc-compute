"""Download IEDB MHC ligand data (presentation-level training/benchmark data).

IEDB offers a well-known static export (mhc_ligand_full) used by NetMHCpan and
many TCR-pMHC papers. We attempt the legacy direct downloader URL first, then
fall back to the database-export page for manual guidance.

Usage:
    python download_iedb.py [--out DIR]
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

DEFAULT_OUT = Path(__file__).resolve().parents[3] / "data" / "raw" / "iedb"

# Legacy direct download (stable for many years, used by NetMHCpan benchmarks).
# CSV variant of mhc_ligand_full; the zip contains the full single file.
CANDIDATES = [
    "https://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full_single_file.zip",
    "https://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full.zip",
]


def download(url: str, dest: Path, timeout: int = 600) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {url} -> {dest}")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (meta-acdc/0.1)"})
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

    for url in CANDIDATES:
        try:
            download(url, args.out / "mhc_ligand_full.zip")
            return 0
        except Exception as e:  # noqa: BLE001
            print(f"{url} failed: {e}", file=sys.stderr)

    print(
        "Automatic download failed. Manual route: https://www.iedb.org/database_export_v3.php "
        "(select MHC ligand -> CSV single file) and place the zip in "
        f"{args.out}/, then rerun the cleaning pipeline.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
