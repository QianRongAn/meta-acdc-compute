"""Proteome-wide candidate selection for the first 50k oligo pool (Q3, KN-8).

Proposal: "select the most informative 0.05% (50,000) of candidates from a
filtered pool of 5 million peptides, prioritized by Expected Information
Gain (EIG)."

Pipeline (stages, resumable):
  prefilter: human Swiss-Prot FASTA -> all 9/10-mers -> MHCflurry
             presentation scores -> top 5M pool -> data/processed/acdc_pool5m.tsv.gz
  rank:      pool + TCR -> 20x bootstrap LR ensemble (VDJdb epitope-split
             train fold, same protocol as simulate_al.py) -> EIG (BALD)
             -> top 50k -> data/processed/kn8_top50k_<tcr>.tsv

The wet-lab ordering (KN-8) is the user's; this produces the ranked
candidate list. Model caveat: EIG from the weak sequence surrogate;
re-rank with the domain-adapted EGNN once TCRmodel2 structures land.

Usage:
    .venv/bin/python src/meta_acdc/active_learning/select_candidates.py --stage prefilter
    .venv/bin/python src/meta_acdc/active_learning/select_candidates.py --stage rank --tcr CASSLGRYNEQFF

Proteome source (data/raw/human_sprot.fasta, gitignored — fetch once):
    curl -sL 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(organism_id:9606)%20AND%20(reviewed:true)' \
        -o data/raw/human_sprot.fasta   # Swiss-Prot human, ~20.4k canonical
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
import time
from pathlib import Path

import numpy as np

ALLELE = "HLA-A*02:01"
POOL_SIZE = 5_000_000
TOP_K = 50_000
CHUNK = 200_000
AA = set("ACDEFGHIKLMNPQRSTVWY")

TCRS = {
    # Kimmtrak/gp100 case (clinical, approved, HLA-A*02:01): natural
    # precursor TCR, VDJdb score 2, PMID:12574392 — not the affinity-enhanced
    # drug construct (its engineered CDR3s are not public).
    "CASSLGRYNEQFF": "gp100/Kimmtrak (YLEPGPVTA) natural precursor, PMID:12574392",
    # A6 TCR (1ao7), well-validated Tax LLFGYPVYV-specific, structure-mapped
    "CAVTTDSWGKLQF": "A6 TCR (Tax LLFGYPVYV), PDB 1ao7, VDJdb-validated",
}

# per-TCR sanity peptides (target + validated cross-reactives), labeled
SANITY = {
    "CASSLGRYNEQFF": [("target YLEPGPVTA", "YLEPGPVTA"),
                      ("off-target YLEPGPVTV*", "YLEPGPVTV"),
                      ("off-target YLEPGPVTL*", "YLEPGPVTL")],
    # *not in Swiss-Prot canonical proteome (literature variants)
    "CAVTTDSWGKLQF": [("target LLFGYPVYV (Tax)", "LLFGYPVYV"),
                      ("cross LLFGYAVYV", "LLFGYAVYV"),
                      ("cross LGYGFVNYI", "LGYGFVNYI"),
                      ("cross LLFGPVYV", "LLFGPVYV")],
}


def iter_fasta(path: Path):
    acc, gene, seq = None, None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    yield acc, gene, "".join(seq)
                head = line[1:]
                acc = head.split("|")[1] if "|" in head else head.split()[0]
                gene = ""
                for tok in head.split():
                    if tok.startswith("GN="):
                        gene = tok[3:]
                seq = []
            else:
                seq.append(line.strip())
    if acc is not None:
        yield acc, gene, "".join(seq)


def iter_peptides(seq: str):
    """All 9/10-mers of a sequence, standard amino acids only."""
    for L in (9, 10):
        for i in range(len(seq) - L + 1):
            p = seq[i:i + L]
            if not (set(p) - AA):
                yield p


def stage_prefilter(fasta: Path, out: Path) -> int:
    """All 9/10-mers -> MHCflurry presentation -> top POOL_SIZE peptides.

    Two-pass memory-safe design: pass 1 scores a systematic ~1M sample to
    calibrate the keep-threshold; pass 2 streams all peptides in chunks,
    keeping only those >= threshold.
    """
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    t0 = time.time()
    proteins = list(iter_fasta(fasta))
    n_total = sum(sum(1 for _ in iter_peptides(s)) for _, _, s in proteins)
    print(f"proteins: {len(proteins)}, standard-AA 9/10-mers: "
          f"{n_total/1e6:.1f}M ({time.time()-t0:.0f}s)", flush=True)

    def predict(peps: list[str]) -> np.ndarray:
        df = predictor.predict(peptides=peps, alleles=[ALLELE], verbose=0)
        return df["presentation_score"].to_numpy()

    # ---- pass 1: calibrate threshold on a systematic sample
    k = max(n_total // 1_000_000, 1)
    sample = []
    c = 0
    for _, _, seq in proteins:
        for p in iter_peptides(seq):
            if c % k == 0:
                sample.append(p)
            c += 1
    svals = np.concatenate(
        [predict(sample[j:j + 500_000]) for j in range(0, len(sample), 500_000)])
    keep_frac = POOL_SIZE / n_total
    thresh = float(np.quantile(svals, 1 - keep_frac))
    print(f"sample {len(sample)} scored; keep threshold {thresh:.4f} "
          f"(keep {keep_frac*100:.1f}% of {n_total/1e6:.1f}M, "
          f"{time.time()-t0:.0f}s)", flush=True)

    # ---- pass 2: stream all peptides in chunks, keep >= thresh
    kept_pep: list[str] = []
    kept_score: list[float] = []
    buf: list[str] = []

    def flush() -> None:
        if not buf:
            return
        sc = predict(buf)
        m = sc >= thresh
        kept_pep.extend(p for p, k2 in zip(buf, m) if k2)
        kept_score.extend(sc[m].tolist())
        buf.clear()

    for pi, (_, _, seq) in enumerate(proteins):
        for p in iter_peptides(seq):
            buf.append(p)
            if len(buf) >= CHUNK:
                flush()
        if (pi + 1) % 4000 == 0:
            print(f"  {pi+1}/{len(proteins)} proteins, kept "
                  f"{len(kept_pep)/1e6:.2f}M ({time.time()-t0:.0f}s)",
                  flush=True)
    flush()
    print(f"kept after pass 2: {len(kept_pep)/1e6:.2f}M peptides "
          f"({time.time()-t0:.0f}s)", flush=True)

    # dedupe (same peptide from multiple proteins) and trim to POOL_SIZE
    order = np.argsort(-np.asarray(kept_score), kind="stable")
    seen: set[str] = set()
    final_pep, final_score = [], []
    for i in order:
        p = kept_pep[i]
        if p in seen:
            continue
        seen.add(p)
        final_pep.append(p)
        final_score.append(kept_score[i])
        if len(final_pep) >= POOL_SIZE:
            break
    print(f"final pool: {len(final_pep)} unique peptides "
          f"(min presentation {final_score[-1]:.4f})", flush=True)

    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["peptide", "presentation"])
        w.writerows(zip(final_pep, (f"{s:.4f}" for s in final_score)))
    print(f"saved -> {out} ({time.time()-t0:.0f}s total)", flush=True)
    return 0


def stage_rank(pool_path: Path, tcr: str, out: Path,
               vdjdb: Path, fasta: Path) -> int:
    from sklearn.linear_model import LogisticRegression
    from meta_acdc.active_learning.acquisition import expected_information_gain
    from meta_acdc.data.splits import group_split
    from meta_acdc.models.benchmark_baseline import featurize

    # ---- pool
    peps, pres = [], []
    with gzip.open(pool_path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            peps.append(row["peptide"])
            pres.append(float(row["presentation"]))
    pres = np.asarray(pres)
    pool_index = {p: i for i, p in enumerate(peps)}
    print(f"pool: {len(peps)} peptides", flush=True)

    # ---- train 20x bootstrap LR ensemble on the epitope-split train fold
    rows = []
    with open(vdjdb, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cdr3_beta"] and r["peptide"] and r["label"]:
                rows.append(r)
    split = group_split([r["cdr3_beta"] for r in rows],
                        [r["peptide"] for r in rows],
                        group_by="epitope", seed=42)
    tr = split.train_idx
    X_tr = np.vstack([featurize(rows[i]["cdr3_beta"], rows[i]["peptide"])
                      for i in tr])
    y_tr = np.array([int(rows[i]["label"]) for i in tr], dtype=float)
    rng = np.random.RandomState(7)
    clfs = []
    for k in range(20):
        idx = rng.randint(0, len(y_tr), len(y_tr))
        clf = LogisticRegression(max_iter=1000, class_weight="balanced")
        clf.fit(X_tr[idx], y_tr[idx])
        clfs.append(clf)
    print(f"ensemble trained on {len(y_tr)} pairs "
          f"({int(y_tr.sum())} positives)", flush=True)

    # ---- EIG over pool in chunks
    t0 = time.time()
    eig = np.zeros(len(peps), dtype=np.float32)
    for c0 in range(0, len(peps), CHUNK):
        c1 = min(c0 + CHUNK, len(peps))
        Xc = np.vstack([featurize(tcr, p) for p in peps[c0:c1]])
        probs = np.vstack([clf.predict_proba(Xc)[:, 1] for clf in clfs])
        eig[c0:c1] = expected_information_gain(probs)
        if (c0 // CHUNK) % 5 == 0:
            print(f"  scored {c1/1e6:.1f}M ({time.time()-t0:.0f}s)",
                  flush=True)

    rank_of = np.argsort(np.argsort(-eig)) + 1  # 1-based rank per pool entry
    top = np.argsort(-eig)[:TOP_K]
    top_set = {peps[i] for i in top}

    # ---- map top peptides back to source proteins
    src: dict[str, list[str]] = {}
    for acc, gene, seq in iter_fasta(fasta):
        for L in (9, 10):
            for i in range(len(seq) - L + 1):
                p = seq[i:i + L]
                if p in top_set:
                    src.setdefault(p, []).append(f"{acc}|{gene}|{i+1}-{i+L}")

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rank", "peptide", "eig", "presentation",
                    "n_sources", "sources"])
        for r, i in enumerate(top, 1):
            p = peps[i]
            s = ";".join(src.get(p, [])[:3])
            w.writerow([r, p, f"{eig[i]:.5f}", f"{pres[i]:.4f}",
                        len(src.get(p, [])), s])
    print(f"top {TOP_K} -> {out}", flush=True)

    # ---- sanity: known target/cross peptides & VDJdb A2 epitope enrichment
    for label, p in SANITY.get(tcr, [("target", "")]):
        if not p:
            continue
        i = pool_index.get(p)
        if i is not None:
            print(f"  {label}: pool rank {rank_of[i]} / {len(peps)} "
                  f"(eig {eig[i]:.4f}, pres {pres[i]:.3f})"
                  f"{'  <-- IN TOP 50K' if rank_of[i] <= TOP_K else ''}",
                  flush=True)
        else:
            print(f"  {label}: NOT IN POOL", flush=True)

    a2_pos = {r["peptide"] for r in rows
              if r["label"] == "1" and r["mhc_allele"] == ALLELE
              and len(r["peptide"]) in (9, 10)}
    in_pool = a2_pos & set(peps)
    in_top = a2_pos & top_set
    from scipy.stats import hypergeom
    n, K, N = len(peps), len(in_pool), TOP_K
    exp = N * K / n
    p_val = float(hypergeom.sf(len(in_top) - 1, n, K, N))
    print(f"  VDJdb {ALLELE} positives: {len(a2_pos)} known, {len(in_pool)} "
          f"in pool, {len(in_top)} in top50k "
          f"(random expectation {exp:.1f}, enrichment "
          f"{len(in_top)/max(exp,1e-9):.1f}x, hypergeom P={p_val:.2e})",
          flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["prefilter", "rank"], required=True)
    ap.add_argument("--fasta", type=Path,
                    default=Path("data/raw/human_sprot.fasta"))
    ap.add_argument("--pool", type=Path,
                    default=Path("data/processed/acdc_pool5m.tsv.gz"))
    ap.add_argument("--vdjdb", type=Path,
                    default=Path("data/processed/vdjdb.clean.tsv"))
    ap.add_argument("--tcr", default="CASSLGRYNEQFF",
                    help=f"TCR beta CDR3; known: {list(TCRS)}")
    args = ap.parse_args()

    if args.stage == "prefilter":
        return stage_prefilter(args.fasta, args.pool)
    name = args.tcr[:8]
    out = Path(f"data/processed/kn8_top50k_{name}.tsv")
    print(f"TCR: {args.tcr} ({TCRS.get(args.tcr, 'custom')})", flush=True)
    return stage_rank(args.pool, args.tcr, out, args.vdjdb, args.fasta)


if __name__ == "__main__":
    sys.exit(main())
