"""STRUCT-FINGERPRINT-20261021: receptor fingerprinting trainer/evaluator.

Pre-registered representation-quality benchmark (docs/fingerprint_retrieval_20261021/PLAN.md).
Trains an EGNN embedding with a supervised-contrastive objective (positives =
same receptor identity) and evaluates on family-grouped folds:

- Task A (primary): receptor retrieval (rank-1, MRR, same-vs-diff AUROC)
- Task B (secondary): family classification via k-NN
- Task C (diagnostic): crystal/predicted domain nearest-centroid

No binding labels, no clinical claims. Results are written verbatim,
whatever they are.

Usage:
    .venv/bin/python -m meta_acdc.models.train_fingerprint \
        --dataset data/processed/fingerprint_dataset.pt \
        --out-dir docs/fingerprint_retrieval_20261021 \
        --epochs 60 --batch-size 16
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

from meta_acdc.models.egnn import EGNN
from meta_acdc.models.train_egnn import normalize_node_features

SEED = 20261021
EMB_DIM = 64
TEMPERATURE = 0.1


def set_seed():
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


class Embedder(nn.Module):
    def __init__(self, node_dim: int, edge_dim: int):
        super().__init__()
        self.egnn = EGNN(node_dim=node_dim, edge_dim=edge_dim, depth=6,
                         hidden=128, update_coordinates=True)
        pooled_dim = 3 * 128 + 1
        self.proj = nn.Sequential(
            nn.Linear(pooled_dim, 128), nn.SiLU(),
            nn.Linear(128, EMB_DIM))
        self.ln = nn.LayerNorm(EMB_DIM)

    def forward(self, h, x, edge_index, edge_attr, batch, peptide_mask):
        pooled = self.egnn(h, x, edge_index, edge_attr, batch,
                           embed_only=True, peptide_mask=peptide_mask)
        return self.ln(self.proj(pooled))


def supcon_loss(embs: torch.Tensor, receptor_ids: list[str],
                temperature: float = TEMPERATURE) -> torch.Tensor:
    """Supervised contrastive loss; positives = same receptor identity.

    Anchors without any same-receptor partner in the batch contribute 0.
    """
    sim = embs @ embs.t() / temperature
    n = embs.shape[0]
    same = torch.zeros(n, n, device=embs.device)
    for i in range(n):
        for j in range(n):
            if i != j and receptor_ids[i] == receptor_ids[j]:
                same[i, j] = 1.0
    has_pos = same.sum(dim=1) > 0
    sim = sim.masked_fill(torch.eye(n, device=embs.device).bool(),
                          float("-inf"))
    pos = (sim.masked_fill(~(same > 0.5), float("-inf"))).logsumexp(dim=1)
    all_sum = sim.logsumexp(dim=1)
    per = -(pos - all_sum)
    per = torch.where(has_pos, per, torch.zeros_like(per))
    n_valid = has_pos.sum().clamp(min=1)
    return per.sum() / n_valid


def batch_graphs(graphs, records, idx, device, model):
    hs, xs, eis, eas, bts, pms, rids = [], [], [], [], [], [], []
    for i in idx:
        nf, coords, ei, ef = graphs[i]
        h = torch.tensor(nf, dtype=torch.float32)
        x = torch.tensor(coords, dtype=torch.float32)
        ei_t = torch.tensor(ei, dtype=torch.long).t()
        ea = torch.tensor(ef, dtype=torch.float32)
        pm = torch.tensor(records[i]["peptide_mask"], dtype=torch.bool)
        h = normalize_node_features(h)
        hs.append(h); xs.append(x); eis.append(ei_t); eas.append(ea)
        bts.append(torch.full((h.shape[0],), len(hs) - 1, dtype=torch.long))
        pms.append(pm)
        rids.append(records[i]["receptor_id"])
    return (torch.cat(hs).to(device), torch.cat(xs).to(device),
            torch.cat(eis, dim=1).to(device), torch.cat(eas).to(device),
            torch.cat(bts).to(device), torch.cat(pms).to(device), rids)


def embed_all(graphs, records, device, model):
    model.eval()
    embs = []
    with torch.no_grad():
        for i in range(len(graphs)):
            (h, x, ei, ea, bt, pm, _) = batch_graphs(
                graphs, records, [i], device, model)
            embs.append(model(h, x, ei, ea, bt, pm).squeeze(0).cpu())
    embs = torch.stack(embs)
    embs = embs / embs.norm(dim=1, keepdim=True).clamp(min=1e-8)
    return embs


def auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Mann-Whitney U AUC for binary labels."""
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    wins = sum(np.sum(p > neg) + 0.5 * np.sum(p == neg) for p in pos)
    return wins / (len(pos) * len(neg))


def evaluate(embs: torch.Tensor, records: list[dict], test_idx: list[int],
             family_of_receptor: dict[str, str]) -> dict:
    sim = embs @ embs.t()  # cosine (normalized)
    n = len(records)
    test = set(test_idx)

    # --- Task A: receptor retrieval -------------------------------------
    rank1s, rrs, chances = [], [], []
    pair_scores, pair_labels = [], []
    for q in test_idx:
        others = [j for j in range(n) if j != q]
        same = [j for j in others
                if records[j]["receptor_id"] == records[q]["receptor_id"]]
        order = sorted(others, key=lambda j: -float(sim[q, j]))
        rank = next((k + 1 for k, j in enumerate(order)
                     if records[j]["receptor_id"] == records[q]["receptor_id"]),
                    None)
        if rank is None:
            continue
        rank1s.append(1.0 if rank == 1 else 0.0)
        rrs.append(1.0 / rank)
        chances.append(len(same) / len(others))
        for j in others:
            pair_scores.append(float(sim[q, j]))
            pair_labels.append(1.0 if records[j]["receptor_id"]
                               == records[q]["receptor_id"] else 0.0)
    retrieval = {
        "rank1": float(np.mean(rank1s)) if rank1s else None,
        "mrr": float(np.mean(rrs)) if rrs else None,
        "chance_rank1": float(np.mean(chances)) if chances else None,
        "n_queries_with_same_receptor_in_gallery": len(rank1s),
        "auroc": auroc(np.array(pair_scores), np.array(pair_labels)),
        "n_pairs": len(pair_scores),
    }

    # --- Task B: family classification (k-NN) ----------------------------
    train_idx = [i for i in range(n) if i not in test]
    fam_counts = Counter(records[i]["family_id"] for i in train_idx)
    chance_top1 = max(fam_counts.values()) / len(train_idx) if train_idx else 0
    fam_acc = {1: [], 3: [], 5: []}
    for q in test_idx:
        others = [j for j in range(n) if j != q]
        order = sorted(others, key=lambda j: -float(sim[q, j]))
        for k in fam_acc:
            votes = Counter(records[j]["family_id"] for j in order[:k])
            top = votes.most_common(1)[0][0] if votes else None
            fam_acc[k].append(1.0 if top == records[q]["family_id"] else 0.0)
    family = {str(k): float(np.mean(v)) if v else None
              for k, v in fam_acc.items()}
    family["chance_top1"] = float(chance_top1)

    # --- Task C: domain diagnostic ---------------------------------------
    dom_cent = {}
    for d in ("crystal", "predicted"):
        idx = [i for i in train_idx if records[i]["domain"] == d]
        if idx:
            dom_cent[d] = embs[idx].mean(0)
    dom_acc = []
    for q in test_idx:
        if len(dom_cent) < 2:
            break
        pred = min(dom_cent, key=lambda d: float(
            1 - (embs[q] @ dom_cent[d])))
        dom_acc.append(1.0 if pred == records[q]["domain"] else 0.0)
    domain = {"accuracy": float(np.mean(dom_acc)) if dom_acc else None,
              "n": len(dom_acc)}

    return {"retrieval": retrieval, "family": family, "domain": domain}


def group_folds(families: list[str], n_splits: int = 5) -> list[list[str]]:
    fams = sorted(set(families))
    rng = random.Random(SEED)
    rng.shuffle(fams)
    folds = [[] for _ in range(n_splits)]
    for i, f in enumerate(fams):
        folds[i % n_splits].append(f)
    return folds


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", type=Path,
                    default=Path("data/processed/fingerprint_dataset.pt"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("docs/fingerprint_retrieval_20261021"))
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--device", default=None)
    ap.add_argument("--folds", type=int, default=5)
    args = ap.parse_args(argv)

    set_seed()
    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device: {device}", flush=True)

    data = torch.load(args.dataset, map_location="cpu", weights_only=False)
    graphs, records = data["graphs"], data["records"]
    # fold assignment is per FAMILY, persisted before training
    families = [r["family_id"] for r in records]
    folds = group_folds(families, args.folds)
    family_to_fold = {f: i for i, fs in enumerate(folds) for f in fs}
    fold_of = [family_to_fold[f] for f in families]

    run_dir = args.out_dir / "model_run"
    run_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "fold_assignment.json").write_text(json.dumps(
        {"seed": SEED, "n_splits": args.folds,
         "family_folds": folds,
         "per_graph_fold": fold_of,
         "graph_records": [r["complex_id"] for r in records]}, indent=2) + "\n")

    fold_results = []
    for fold in range(args.folds):
        test_idx = [i for i, f in enumerate(fold_of) if f == fold]
        train_idx = [i for i in range(len(graphs)) if i not in test_idx]
        print(f"fold {fold}: train={len(train_idx)} test={len(test_idx)} "
              f"test_families={len(set(families[i] for i in test_idx))}",
              flush=True)
        model = Embedder(26, 16).to(device)
        opt = torch.optim.Adam(model.parameters(), lr=3e-4, weight_decay=1e-4)
        model.train()
        for epoch in range(args.epochs):
            order = train_idx[:]
            random.shuffle(order)
            total, nb = 0.0, 0
            for s in range(0, len(order), args.batch_size):
                bidx = order[s:s + args.batch_size]
                h, x, ei, ea, bt, pm, rids = batch_graphs(
                    graphs, records, bidx, device, model)
                embs = model(h, x, ei, ea, bt, pm)
                loss = supcon_loss(embs, rids, TEMPERATURE)
                if not torch.isfinite(loss):
                    raise RuntimeError(
                        f"fold {fold} epoch {epoch} batch {s}: "
                        f"non-finite SupCon loss; training diverged")
                opt.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
                opt.step()
                total += float(loss.item())
                nb += 1
            if epoch % 15 == 0 or epoch == args.epochs - 1:
                print(f"  fold {fold} epoch {epoch}: loss {total / max(nb, 1):.4f}",
                      flush=True)
        torch.save({"model": model.state_dict(), "fold": fold, "seed": SEED},
                   run_dir / f"fold{fold}.pt")
        embs = embed_all(graphs, records, device, model)
        res = evaluate(embs, records, test_idx, None)
        res["fold"] = fold
        res["n_train"] = len(train_idx)
        res["n_test"] = len(test_idx)
        fold_results.append(res)
        partial = args.out_dir / "results_partial.json"
        partial.write_text(json.dumps(
            {"task_id": "STRUCT-FINGERPRINT-20261021",
             "fold_results": fold_results}, indent=2, allow_nan=False) + "\n")
        print(json.dumps({k: v for k, v in res.items() if k in
                          ("fold", "retrieval", "family")}, indent=2))

    agg = {
        "task_id": "STRUCT-FINGERPRINT-20261021",
        "n_graphs": len(graphs),
        "seed": SEED,
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "device": device,
        "fold_results": fold_results,
        "means": {
            "retrieval_rank1": float(np.nanmean(
                [r["retrieval"]["rank1"] for r in fold_results])),
            "retrieval_mrr": float(np.nanmean(
                [r["retrieval"]["mrr"] for r in fold_results])),
            "retrieval_chance_rank1": float(np.nanmean(
                [r["retrieval"]["chance_rank1"] for r in fold_results])),
            "retrieval_auroc": float(np.nanmean(
                [r["retrieval"]["auroc"] for r in fold_results])),
            "family_top1_k1": float(np.nanmean(
                [r["family"]["1"] for r in fold_results])),
            "family_chance_top1": float(np.nanmean(
                [r["family"]["chance_top1"] for r in fold_results])),
            "domain_accuracy": float(np.nanmean(
                [r["domain"]["accuracy"] for r in fold_results])),
        },
        "binding_labels": False,
        "interpretation_scope": (
            "representation-quality benchmark only; receptor identity is "
            "verified metadata, not a binding label; no binding prediction "
            "or clinical claim follows"),
    }
    out = args.out_dir / "results.json"
    if out.exists():
        raise ValueError("refusing to overwrite existing results.json")
    out.write_text(json.dumps(agg, indent=2, allow_nan=False) + "\n")
    print(json.dumps(agg["means"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
