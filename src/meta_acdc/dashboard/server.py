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
from functools import partial
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]  # repo root
STATIC = Path(__file__).resolve().parent / "static"
# prefer the relabeled table (2026-09-23 chain-fingerprint correction)
DEFAULT_SCORES = ROOT / "data/processed/prediction_scores_ensemble.relabeled.tsv"
if not DEFAULT_SCORES.exists():
    DEFAULT_SCORES = ROOT / "data/processed/prediction_scores_ensemble.tsv"
DEFAULT_MAP = ROOT / "data/processed/structure_vdjdb_map.tsv"
DEFAULT_CLINICAL = ROOT / "data/processed/clinical_gold_standard.tsv"
DEFAULT_QC = ROOT / "data/processed/af3_qc.tsv"
DEFAULT_CLINICAL_SCAN = ROOT / "data/processed/kn11_clinical_scan.tsv"
NATIVE_TODO = ROOT / "data/processed/af3_native_todo.txt"
UPLOAD_TODO = ROOT / "data/processed/af3_upload_todo.txt"
NATIVE_MANIFEST = ROOT / "data/processed/af3_native/manifest.tsv"


def load_reliability(qc_path: Path) -> dict:
    """AF3 structural reliability gate (ipTM + chain completeness)."""
    if not qc_path.exists():
        return {}
    recs = []
    with open(qc_path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            recs.append(row)
    flagged = [r for r in recs if r.get("flag")]
    return {
        "n": len(recs),
        "n_flagged": len(flagged),
        "n_low_iptm": sum(1 for r in recs if "low_iptm" in r.get("flag", "")),
        "n_truncated": sum(1 for r in recs if "truncated_tcr" in r.get("flag", "")),
        "gate": "ipTM < 0.75 or TCR CA < 400",
        "examples": [{"instance": r["instance"], "job": r["job_id"],
                       "iptm": r["iptm_mean"], "tcr_ca": r["tcr_ca"],
                       "flag": r["flag"]} for r in flagged[:20]],
    }


def load_clinical_scan(path: Path) -> list[dict]:
    """No clinical scoring model is released; old scans never enter the UI."""
    return []


def load_submission_status() -> dict:
    """How many AF3 jobs the user still needs to submit (native / candidate)."""

    def lines(p: Path) -> set[str]:
        return {ln.strip() for ln in open(p)} if p.exists() else set()

    native_todo = lines(NATIVE_TODO)
    native_done = set()
    if NATIVE_MANIFEST.exists():
        with open(NATIVE_MANIFEST, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                native_done.add(r["job_id"])
    upload_todo = lines(UPLOAD_TODO)
    scored = set()
    for p in (ROOT / "data/processed").glob("da_scores_seed0_ensemble.tsv"):
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                scored.add(r["job_id"])
    return {
        "native_total": len(native_todo),
        "native_done": len(native_done & native_todo),
        "native_remaining": len(native_todo - native_done),
        "candidate_total": len(upload_todo),
        "candidate_done": len(scored),
        "candidate_remaining": len(upload_todo - scored),
    }


def build_data(scores_path: Path, map_path: Path, clinical_path: Path,
               qc_path: Path | None = None,
               clinical_scan_path: Path | None = None) -> dict:
    """Incident/reference view only; revoked scores never become rankings."""
    from meta_acdc.structure.artifact_status import catalog
    state=catalog()
    references=[]
    if clinical_path.exists():
        with clinical_path.open() as f:
            references=[r for r in csv.DictReader(f,delimiter='\t')
                        if r.get('record_status')=='REFERENCE_ONLY_NOT_GOLD_STANDARD']
    return dict(groups={}, clinical=[], clinical_scan=[], cross_evidence={},
        same_evidence={}, clinical_reference_facts=references,
        reliability=dict(load_reliability(qc_path) if qc_path else {},
                         status='historical_job_level_QC_not_release_authorization'),
        submission_status=dict(load_submission_status(),status='historical_counts_only'),
        artifact_status=state,
        domain_adaptation=dict(candidate_verdict='WITHDRAWN_UNVALIDATED',
                               note='No structural scoring or clinical model is released.'),
        caveats=['历史结构分数、临床危险标签及排名已撤回，当前不展示。',
                 'MAG-IC3 不是临床 a3a；临床参考事实不构成模型验收。',
                 '新亲和力研究基线及冻结验收集另行管理；本页不作风险判断。'])


class Handler(SimpleHTTPRequestHandler):
    data: dict = {}

    def do_GET(self):  # noqa: N802
        if self.path in ("/", ""):
            self.path = "/index.html"
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
    ap.add_argument("--qc", type=Path, default=DEFAULT_QC)
    ap.add_argument("--clinical-scan", type=Path,
                    default=DEFAULT_CLINICAL_SCAN)
    args = ap.parse_args()

    Handler.data = build_data(args.scores, args.map, args.clinical, args.qc,
                              args.clinical_scan)

    # Python 3.12+: SimpleHTTPRequestHandler.__init__ ignores the class
    # attribute and defaults to os.getcwd() — must pass directory= here.
    server = HTTPServer(("127.0.0.1", args.port),
                        partial(Handler, directory=str(STATIC)))
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
