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
    """KN-11 structural safety scan over clinical gold-standard cases."""
    if not path.exists():
        return []
    rows = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("role") == "no structure":
                continue
            rows.append({
                "case": r["case"], "family": r["family"],
                "peptide": r["peptide"], "role": r["role"],
                "score": float(r["score"]) if r.get("score") else None,
                "rank": r.get("rank", ""), "fatal": r["fatal"],
                "flag": r.get("flag", ""),
            })
    return rows


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
    # cross-reactivity evidence per PDB: VDJdb-validated epitopes that differ
    # from the crystallized peptide (KN-4+ map)
    cross: dict[str, list[str]] = {}
    same: dict[str, list[str]] = {}
    for (pdb, epi), n in evidence.items():
        if epi == native.get(pdb):
            same.setdefault(pdb, []).append(epi)
        else:
            cross.setdefault(pdb, []).append(epi)

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
        "cross_evidence": {pdb: sorted(set(v)) for pdb, v in cross.items()},
        "same_evidence": {pdb: sorted(set(v)) for pdb, v in same.items()},
        "reliability": load_reliability(qc_path) if qc_path else {},
        "clinical_scan": load_clinical_scan(clinical_scan_path)
        if clinical_scan_path else [],
        "submission_status": load_submission_status(),
        "domain_adaptation": {
            "n_native_positives": 144,
            "val_auroc": [0.959, 0.952, 0.964],
            "val_cross_seed_spearman": [0.943, 0.939, 0.947],
            "candidate_verdict": "STABLE",
            "candidate_mean_spearman": 0.729,
            "note": "144 个真实天然 AF3 正样本域适应;val AUROC 0.96,候选"
                    "判决 mean r=0.729(STABLE,过 0.5 门)——Rashomon 危机"
                    "解除,扩样单调(0.447→0.479→0.729),见 benchmarks.md",
        },
        "caveats": [
            "AF3 交叉反应排名已跨实例稳定(144 native 域适应,mean pairwise "
            "Spearman 0.729,过 0.5 门)——排名可报告,但须用固定多实例集成;"
            "单实例 EGNN 打分仍不可信。详见 benchmarks.md 2026-09-23",
            "打分须报告固定多实例集成;单实例 EGNN 打分不可信(非 AF3 不稳)。"
            "重提分析(31 任务):高置信(ipTM≥0.88)跨重提稳,低 ipTM(~0.5)"
            " AF3 自身会飘——见 reliability 面板",
            "pdb 标签可靠性警告(2026-09-23):kn5 提交清单 tcr_a/tcr_b 两列"
            "写反,历史脚本曾把 48/73 条结构标签张冠李戴;下表的 TCR 分组名"
            "请按链指纹家族理解,待结构映射修复",
        ],
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
