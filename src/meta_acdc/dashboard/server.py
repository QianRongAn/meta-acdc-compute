"""TCR-Safety-Radar dashboard server (KN-13, v0).

Zero-dependency web server (stdlib only) that renders the cross-reactivity
ranking results as an interactive page:
- per-TCR ranking tables (compatibility scores on AF3-predicted structures)
- TCR x peptide heatmap (traffic-light colors)
- gold-standard clinical cases panel (MAGE-A3/titin etc.)

Usage:
    .venv/bin/python src/meta_acdc/dashboard/server.py --port 8000
then open http://localhost:8000
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]  # repo root
STATIC = Path(__file__).resolve().parent / "static"
DEFAULT_SCORES = ROOT / "data/processed/prediction_scores.tsv"
DEFAULT_MAP = ROOT / "data/processed/structure_vdjdb_map.tsv"
DEFAULT_CLINICAL = ROOT / "data/processed/clinical_gold_standard.tsv"


def build_data(scores_path: Path, map_path: Path, clinical_path: Path) -> dict:
    scores = []
    if scores_path.exists():
        with open(scores_path, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["score"]:
                    scores.append({
                        "job": row["job_id"].rsplit("_model", 1)[0],
                        "score": float(row["score"]),
                        "std": float(row["std"]) if row.get("std") else None,
                    })
    # per-TCR groups
    groups: dict[str, list] = {}
    for s in scores:
        pdb, pep = s["job"].split("_", 1)
        groups.setdefault(pdb, []).append({"peptide": pep, **s})
    for v in groups.values():
        v.sort(key=lambda r: -r["score"])

    native = {}
    evidence = {}
    if map_path.exists():
        with open(map_path, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                native.setdefault(row["pdb"], row["pdb_peptide"])
                evidence.setdefault((row["pdb"], row["vdjdb_epitope"]), 0)
                evidence[(row["pdb"], row["vdjdb_epitope"])] += 1

    clinical = []
    if clinical_path.exists():
        with open(clinical_path, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                clinical.append({
                    "tcr": row["tcr"], "target": row["target"],
                    "off_targets": row["off_targets"].split(";"),
                    "fatal": row["fatal"], "evidence": row["evidence"],
                    "note": row["note"],
                })

    return {
        "groups": {pdb: {"native": native.get(pdb, ""), "rows": rows}
                   for pdb, rows in sorted(groups.items())},
        "clinical": clinical,
    }


class Handler(SimpleHTTPRequestHandler):
    data: dict = {}

    def do_GET(self):  # noqa: N802
        if self.path.startswith("/api/data"):
            body = json.dumps(self.data).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        super().do_GET()

    def log_message(self, *args):
        pass  # quiet


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--scores", type=Path, default=DEFAULT_SCORES)
    ap.add_argument("--map", type=Path, default=DEFAULT_MAP)
    ap.add_argument("--clinical", type=Path, default=DEFAULT_CLINICAL)
    args = ap.parse_args()

    Handler.data = build_data(args.scores, args.map, args.clinical)
    Handler.directory = str(STATIC)

    server = HTTPServer(("127.0.0.1", args.port), Handler)
    print(f"TCR-Safety-Radar: http://127.0.0.1:{args.port}")
    print(f"  scored TCRs: {len(Handler.data['groups'])}, "
          f"clinical cases: {len(Handler.data['clinical'])}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
