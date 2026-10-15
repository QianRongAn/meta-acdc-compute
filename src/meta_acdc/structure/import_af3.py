"""Import AF3 web-download results into the pipeline (KN-5).

The AF3 web UI downloads one folder per job with model_0..4 CIFs and
confidence JSONs. This script:
1. Reads each job_request.json to recover the submitted sequences
2. Matches them against kn5_submission_list.tsv (peptide + TCR alpha)
3. Copies model_0 (top-ranked) CIF + summary_confidences_0.json into
   data/raw/af3_predictions/{job_id}_model_0.cif (+ _conf.json)

Usage:
    .venv/bin/python src/meta_acdc/structure/import_af3.py \
        --src "/mnt/c/Users/Administrator/Downloads/folds_2026_09_21_11_38" \
        --list data/processed/kn5_submission_list.tsv \
        --out data/raw/af3_predictions
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path


def read_job_sequences(request_path: Path) -> list[str]:
    with open(request_path) as fh:
        req = json.load(fh)
    if isinstance(req, list):
        req = req[0]  # web downloads wrap the job in a list
    seqs = []
    for item in req.get("sequences", []):
        for kind, spec in item.items():
            if isinstance(spec, dict) and "sequence" in spec:
                seqs.append(spec["sequence"])
    return sorted(seqs, key=len)  # peptide first


def tcr_pair(seqs: list[str], peptide: str) -> tuple[str | None, str | None]:
    """Alpha/beta by conserved pre-CDR3 motif (order-independent).

    Fallback: unusual TRAJ/TRBJ alleles miss the hardcoded motifs entirely
    (e.g. 3utt alpha has no YFCAx); then assign by chain length windows
    (class I: alpha ~170-225, beta ~225-265, MHC ~270+, b2m ~100).
    Motif path is tried first, so previously matched jobs are unchanged.
    """
    rest = [s for s in seqs if s != peptide]
    beta = next((s for s in rest if "YFCAS" in s or "YLCAS" in s), None)
    alpha = next((s for s in rest
                  if any(m in s for m in ("YFCAV", "YFCAI", "YFCAL", "YLCAV"))),
                 None)
    if alpha is None or beta is None:
        cands = [s for s in rest if s not in (alpha, beta)]
        if beta is None:
            beta = next((s for s in cands if 225 <= len(s) <= 265), None)
        if alpha is None:
            alpha = next((s for s in cands if 170 <= len(s) <= 225), None)
    return alpha, beta


def matches_chain_pair(alpha, beta, cand: dict) -> bool:
    """True iff the job's {alpha,beta} equals the candidate's {tcr_a,tcr_b}.

    kn5_submission_list.tsv stores tcr_a = beta and tcr_b = alpha (columns
    are swapped relative to the AF3 job FASTA), so an order-dependent
    alpha == tcr_a test never fires — hence the set comparison."""
    if not alpha or not beta:
        return False
    return {alpha, beta} == {cand.get("tcr_a"), cand.get("tcr_b")}


def _eq_upto_termini(s1: str, s2: str, max_trim: int = 2) -> bool:
    """Terminal-tolerant equality for AF3 pasted chains.

    AF3 web pastes sometimes shift chains by 1-2 terminal residues
    (observed: prepend 1 / append 1 / drop last residue), which breaks exact
    fingerprint matching. True iff s1/s2 agree on an identical core with at
    most one head-side trim and at most one tail-side trim TOTAL (max_trim
    residues combined) — so mid-sequence differences can never be masked.
    Peptides must still be compared exactly by the caller.
    """
    if s1 == s2:
        return True
    if abs(len(s1) - len(s2)) > max_trim:
        return False
    for h1 in range(max_trim + 1):              # head trim on s1
        for t1 in range(max_trim - h1 + 1):     # tail trim on s1
            core = s1[h1:len(s1) - t1 if t1 else None]
            h2max = max_trim - h1 - t1
            if h1 > 0:
                h2max = 0                       # head trim on one side only
            for h2 in range(h2max + 1):         # head trim on s2
                t2 = len(s2) - len(core) - h2   # tail trim on s2
                if t2 < 0 or t2 > max_trim - h1 - t1 - h2:
                    continue
                if t1 > 0 and t2 > 0:
                    continue                    # tail trim on one side only
                if core == s2[h2:len(s2) - t2 if t2 else None]:
                    return True
    return False


def matches_chain_pair_tol(alpha, beta, cand: dict, max_trim: int = 2) -> bool:
    """Tolerant version of matches_chain_pair (terminal edits allowed)."""
    if not alpha or not beta:
        return False
    a2, b2 = cand.get("tcr_a"), cand.get("tcr_b")
    if not a2 or not b2:
        return False
    return ((_eq_upto_termini(alpha, a2, max_trim)
             and _eq_upto_termini(beta, b2, max_trim))
            or (_eq_upto_termini(alpha, b2, max_trim)
                and _eq_upto_termini(beta, a2, max_trim)))


def load_manifest(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, required=True)
    ap.add_argument("--list", type=Path,
                    default=Path("data/processed/kn5_submission_list.tsv"))
    ap.add_argument("--out", type=Path, default=Path("data/raw/af3_predictions"))
    ap.add_argument("--unique-fallback", action="store_true",
                    help="when no exact chain-pair match exists, assign a "
                         "unique sequence-derived id instead of quarantining "
                         "(use for native submissions)")
    ap.add_argument("--allow-first-candidate-fallback", action="store_true",
                    help="DANGEROUS: keep the legacy behaviour of naming an "
                         "unmatched job after candidates[0]. This silently "
                         "mislabels structures; only use when you accept the "
                         "attribution risk. Off by default since 2026-10-15.")
    ap.add_argument("--report", type=Path, default=None,
                    help="write a TSV mapping job_id -> source folder, peptide, "
                         "evidence (for building a native manifest)")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(args.list)
    by_peptide: dict[str, list[dict]] = {}
    for row in manifest:
        by_peptide.setdefault(row["peptide"], []).append(row)

    imported = []
    unmatched = []
    quarantined = []
    for job_dir in sorted(args.src.iterdir()):
        if not job_dir.is_dir():
            continue
        reqs = sorted(job_dir.glob("*_job_request.json"))
        if not reqs:
            continue
        req = reqs[0]
        seqs = read_job_sequences(req)
        if not seqs:
            continue
        peptide = seqs[0]
        alpha, beta = tcr_pair(seqs, peptide)
        # clinical gold-standard set (MAG-IC3/5brz TCR): detect FIRST by the
        # exact alpha-chain sequence, before manifest matching
        CLINICAL_ALPHA = (
            "AQEVTQIPAALSVPEGENLVLNCSFTDSAIYNLQWFRQDPGKGLTSLLYVRPYQREQTSGRLNASLDKS"
            "SGRSTLYIAASQPGDSATYLCAVRPGGAGPFFVVFGKGTKLSVIPNIQNPDPAVYQLRDSKSSDKSVCL"
            "FTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIP")
        if alpha == CLINICAL_ALPHA:
            match = {"pdb": "5brz", "evidence": "clinical"}
        else:
            # match manifest: same peptide + same {alpha,beta} chain pair
            candidates = by_peptide.get(peptide, [])
            match = next((c for c in candidates
                          if matches_chain_pair(alpha, beta, c)), None)
            if match is None:
                # No exact chain pair in the submission list. Guessing a name
                # here is what produced the 2026-10-14 misattribution, so the
                # default is to quarantine the job and say so loudly.
                if args.unique_fallback or not candidates:
                    key = "|".join(sorted(seqs))  # full chain set: collision-safe
                    stem = "_".join([peptide,
                                     hashlib.sha1(key.encode()
                                                  ).hexdigest()[:8]])
                    match = {"pdb": f"nativetcr-{stem}", "evidence": "native"}
                elif args.allow_first_candidate_fallback:
                    match = candidates[0]
                else:
                    quarantined.append((job_dir.name, peptide, alpha[:12],
                                        beta[:12]))
                    continue
        job_id = f"{match['pdb']}_{peptide}"
        # copy ALL five ranked models (model_0..4) for structural ensembling
        n_copied = 0
        for cif_src in sorted(job_dir.glob("*_model_?.cif")):
            model_tag = cif_src.name.split("_model_")[1].split(".")[0]
            shutil.copy2(cif_src, args.out / f"{job_id}_model_{model_tag}.cif")
            n_copied += 1
        conf_srcs = sorted(job_dir.glob("*_summary_confidences_0.json"))
        if conf_srcs:
            shutil.copy2(conf_srcs[0], args.out / f"{job_id}_conf.json")
        imported.append((job_dir.name, job_id, peptide,
                         match["evidence"], n_copied))

    print(f"imported {len(imported)} jobs:")
    for src, job_id, pep, ev, n in imported:
        print(f"  {job_id:30s} <- {src}  [{ev}] {n} models")
    if unmatched:
        print(f"unmatched {len(unmatched)}: {unmatched}")
    if quarantined:
        print(f"QUARANTINED {len(quarantined)} jobs: no exact chain-pair "
              f"match in the submission list. No files were written for "
              f"them. Inspect manually before assigning a job_id.", flush=True)
        for src, pep, a, b in quarantined:
            print(f"  {src}: peptide={pep} alpha={a}... beta={b}...")
        qpath = args.out.parent / "af3_import_quarantine.tsv"
        with open(qpath, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["source", "peptide", "alpha_prefix", "beta_prefix"])
            w.writerows(quarantined)
        print(f"wrote quarantine list to {qpath}")
    if args.report is not None:
        with open(args.report, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["job_id", "source", "peptide", "evidence"])
            for src, job_id, pep, ev, _ in imported:
                w.writerow([job_id, src, pep, ev])
        print(f"wrote report to {args.report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
