"""Freeze the candidate panel registry (review items #16/#33).

One canonical口径 per artifact set, computed from the files themselves:

- scored panel:  data/raw/af3_predictions/*_model_0.cif (188 files)
- pending panel: af3_pending_upload/**/*.fasta        (347 files)

For the scored panel, the receptor unit is the unordered {chain A, chain B}
sequence pair from the AF3 CIF chain labels. This reproduces the delivery
audit counts (188 files / 72 PDB prefixes / 59 exact A/B pairs) and must not
be read as biologically assigned roles: chain letters are AF3 output labels;
role assignment requires the local release metadata
(data/releases/*/structure_records.json, not in Git).

For the pending panel, roles come from the FASTA role headers
(>alpha/>beta/>peptide/>mhc/>b2m) — the same rule verified in
docs/forensic_followup_20261021, which reproduces the review counts
(347 files / 111 PDB prefixes / 91 peptides / 98 exact alpha-beta pairs /
45 MHC / 11 B2M sequences). Role headers are metadata, not sequence-truth:
NA/GS/MG prefixes are reported, never trimmed (review item #17).

Usage:
    .venv/bin/python scripts/build_panel_registry.py \
        --cifs data/raw/af3_predictions \
        --fastas af3_pending_upload \
        --out data/processed/candidate_panel_registry.tsv \
        --summary data/processed/candidate_panel_registry.summary.json
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from meta_acdc.data.release import sha256  # noqa: E402
from meta_acdc.structure.cif import parse_cif  # noqa: E402
from meta_acdc.structure.identity import chain_sequence  # noqa: E402
from meta_acdc.structure.graph import Residue  # noqa: E402


def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()


def chain_rows(residues: list[Residue]) -> dict[str, dict]:
    """Per-chain summary from parsed residues (sequence + sha256 + length)."""
    by_chain: dict[str, list[Residue]] = defaultdict(list)
    for r in residues:
        by_chain[r.chain].append(r)
    out = {}
    for c, rs in sorted(by_chain.items()):
        seq = chain_sequence(rs)
        out[c] = {"length": len(seq), "sequence": seq, "sha256": seq_hash(seq)}
    return out


def pct_identity(a: str, b: str) -> float:
    """Global sequence identity over the shorter length (no alignment)."""
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return sum(x == y for x, y in zip(a, b)) / n


def cluster_98(keys: list[str], seqs: dict[str, str]) -> dict[str, str]:
    """Single-linkage clustering of exact keys at >=98% sequence identity.

    Returns key -> cluster id (the lexicographically smallest member key).
    """
    parent = {k: k for k in keys}

    def find(k: str) -> str:
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            if pct_identity(seqs[a], seqs[b]) >= 0.98:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
    return {k: find(k) for k in keys}


def read_fasta(path: Path) -> list[str]:
    return [ln.strip() for ln in open(path, encoding="utf-8", errors="replace")
            if ln.strip() and not ln.startswith(">")]


_ROLES = {"alpha", "beta", "peptide", "mhc", "b2m"}


def read_fasta_roles(path: Path) -> dict[str, str]:
    """Parse a role-annotated FASTA (>alpha/>beta/>peptide/>mhc/>b2m).

    The pending panel FASTAs carry explicit role headers; this is the same
    rule verified in docs/forensic_followup_20261021 (98 exact alpha/beta
    pairs, 45 MHC, 11 B2M sequences). Role headers are metadata, not
    sequence-truth: NA/GS/MG prefixes are reported, never trimmed.
    """
    records: dict[str, str] = {}
    key = None
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            key = line[1:]
            if key in records:
                raise ValueError(f"{path}: duplicate role header {key}")
            records[key] = ""
        else:
            if key is None:
                raise ValueError(f"{path}: sequence before any role header")
            records[key] += line
    if set(records) != _ROLES:
        raise ValueError(
            f"{path}: expected role headers {sorted(_ROLES)}, "
            f"got {sorted(records)}")
    return records


def build_scored_registry(cifs_dir: Path, qc_path: Path | None) -> list[dict]:
    rows = []
    files = sorted(cifs_dir.glob("*_model_0.cif"))
    for cif in files:
        job = cif.name.rsplit("_model_0", 1)[0]
        pdb, _, peptide = job.partition("_")
        try:
            chains = chain_rows(parse_cif(cif))
        except Exception as exc:  # keep the row; mark the parse failure
            chains = {}
            parse_error = f"{type(exc).__name__}: {exc}"
        else:
            parse_error = ""
        ab = tuple(sorted(chains.get(c)["sha256"] for c in ("A", "B")
                          if c in chains))
        ab_seqs = tuple(sorted(chains.get(c)["sequence"] for c in ("A", "B")
                               if c in chains))
        rows.append({
            "panel": "scored",
            "job_id": job,
            "pdb": pdb,
            "peptide": peptide,
            "chain_lengths": ";".join(
                f"{c}:{chains[c]['length']}" for c in sorted(chains)),
            "receptor_key": "|".join(ab) if len(ab) == 2 else "",
            "_receptor_seqs": "|".join(ab_seqs) if len(ab_seqs) == 2 else "",
            "receptor_identity_note":
                "unordered {chain A, chain B} sequence pair from AF3 chain "
                "labels; biological roles require release metadata",
            "qc_flag": _historical_qc(qc_path, job),
            "parse_error": parse_error,
            "file_sha256": sha256(cif),
        })
    return rows


def _historical_qc(qc_path: Path | None, job: str) -> str:
    if qc_path is None or not qc_path.exists():
        return ""
    with open(qc_path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("job_id") == job:
                return row.get("flag", "")
    return ""


def build_pending_registry(fasta_dir: Path) -> list[dict]:
    rows = []
    for fasta in sorted(fasta_dir.rglob("*.fasta")):
        job = fasta.stem
        pdb, _, peptide = job.partition("_")
        try:
            seqs = read_fasta_roles(fasta)
        except ValueError as exc:
            # keep the row, mark the failure; roles stay empty
            seqs = {}
            parse_error = f"{type(exc).__name__}: {exc}"
        else:
            parse_error = ""
        pair = (seqs.get("alpha", ""), seqs.get("beta", ""))
        ok = all(pair) and set("".join(pair)) <= set("ACDEFGHIKLMNPQRSTVWY")
        receptor_key = ("|".join(sorted(seq_hash(s) for s in pair))
                        if ok else "")
        rows.append({
            "panel": "pending",
            "job_id": job,
            "pdb": pdb,
            "peptide": peptide,
            "chain_lengths": ";".join(
                f"{k}:{len(v)}" for k, v in sorted(seqs.items())),
            "chain_set_sha256": sha256_bytes(
                "\x00".join(sorted(seqs.values())).encode()) if seqs else "",
            "tcr_window_pair": "",
            "receptor_key": receptor_key,
            "_receptor_seqs": "|".join(sorted(pair)) if ok else "",
            "_receptor_key": receptor_key,
            "_mhc_seqs": seqs.get("mhc", ""),
            "_b2m_seqs": seqs.get("b2m", ""),
            "receptor_identity_note":
                "roles come from the FASTA role headers (alpha/beta/peptide/"
                "mhc/b2m), verified in docs/forensic_followup_20261021; "
                "NA/GS/MG prefixes are reported, never trimmed; role "
                "headers are metadata, not sequence-truth",
            "qc_flag": "",
            "parse_error": parse_error,
            "file_sha256": sha256(fasta),
        })
    return rows


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def summarize(rows: list[dict], panel: str) -> dict:
    sub = [r for r in rows if r["panel"] == panel]
    key_col = "receptor_key" if panel == "scored" else "_receptor_key"
    exact = Counter(r[key_col] for r in sub if r[key_col])
    seqs = {r[key_col]: r.get("_receptor_seqs", "") for r in sub
            if r[key_col] and r.get("_receptor_seqs")}
    clusters = cluster_98(list(seqs), seqs)
    by_cluster = Counter(clusters.get(k, k) for k in exact)
    extra = {}
    if panel == "pending":
        chain_sets = Counter(r["chain_set_sha256"] for r in sub
                             if r["chain_set_sha256"])
        mhc = Counter(r["_mhc_seqs"] for r in sub if r["_mhc_seqs"])
        b2m = Counter(r["_b2m_seqs"] for r in sub if r["_b2m_seqs"])
        extra["exact_chain_sets"] = len(chain_sets)
        extra["unique_mhc_sequences"] = len(mhc)
        extra["unique_b2m_sequences"] = len(b2m)
        extra["chain_set_identity_note"] = (
            "frozenset of all 5 role sequences per FASTA; roles from FASTA "
            "headers, verified in docs/forensic_followup_20261021")
        extra["percent98_cluster_identity_note"] = (
            ">=98% single-linkage over the concatenated alpha|beta pair is "
            "ONE descriptive rule, not a biological family; the forensic "
            "follow-up shows rule sensitivity (89 vs 57 groups under two "
            "both-chains rules), so no single cluster count is canonical")
    return {
        "panel": panel,
        "files": len(sub),
        "pdb_prefixes": len({r["pdb"] for r in sub}),
        "unique_peptides": len({r["peptide"] for r in sub}),
        "exact_receptor_units": len(exact),
        "receptor_unit": ("unordered {chain A, chain B} sequence pair"
                          if panel == "scored"
                          else "unordered (alpha, beta) pair from FASTA "
                               "role headers"),
        "receptor_units_with_multiple_candidates":
            sum(1 for v in exact.values() if v > 1),
        "largest_receptor_cluster": max(exact.values(), default=0),
        "percent98_clusters": len(by_cluster),
        "percent98_cluster_identity_note":
            ">=98% identity over the receptor unit; single-linkage; "
            "descriptive only",
        "parse_errors": sum(1 for r in sub if r["parse_error"]),
        **extra,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cifs", type=Path,
                    default=Path("data/raw/af3_predictions"))
    ap.add_argument("--fastas", type=Path,
                    default=Path("af3_pending_upload"))
    ap.add_argument("--qc", type=Path,
                    default=Path("data/processed/af3_qc.tsv"),
                    help="legacy job-level QC flags (historical; joined by "
                         "job_id for reference only)")
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/candidate_panel_registry.tsv"))
    ap.add_argument("--summary", type=Path,
                    default=Path(
                        "data/processed/candidate_panel_registry.summary.json"))
    args = ap.parse_args()
    if args.out.exists() or args.summary.exists():
        raise ValueError("refusing to overwrite an existing panel registry")
    scored = build_scored_registry(args.cifs, args.qc)
    pending = build_pending_registry(args.fastas)
    fieldnames = ["panel", "job_id", "pdb", "peptide", "chain_lengths",
                  "receptor_key", "chain_set_sha256", "tcr_window_pair",
                  "receptor_identity_note", "qc_flag", "parse_error",
                  "file_sha256"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for row in scored + pending:
            w.writerow(row)
    summary = {
        "generated_from": {
            "cifs": str(args.cifs),
            "fastas": str(args.fastas),
        },
        "canonical_scope": {
            "scored_panel": "188 model_0 CIFs; receptor unit = unordered "
                            "{chain A, chain B} pair (AF3 labels)",
            "pending_panel": "347 submission FASTAs; no role assignment; "
                             "receptor unit = full chain-set",
            "warning": "these counts describe corpus redundancy, not "
                       "independent biological observations; stability CIs "
                       "must cluster by true receptor, not PDB prefix",
        },
        "panels": [summarize(scored, "scored"),
                   summarize(pending, "pending")],
    }
    with open(args.summary, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    for panel in summary["panels"]:
        print(f"{panel['panel']}: files={panel['files']} "
              f"pdb={panel['pdb_prefixes']} "
              f"exact_receptors={panel['exact_receptor_units']} "
              f"clusters98={panel['percent98_clusters']}")
    print(f"wrote {args.out} and {args.summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
