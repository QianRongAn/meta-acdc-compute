"""FOLLOWUP-AUROC-RASHOMON-20261021: AUROC-tightening recipe (v2) + Rashomon
instances for cross-instance stability analysis.

Pre-registered in docs/fingerprint_retrieval_20261021/FOLLOWUP_AUROC_RASHOMON_20261021.md
and followup_plan.json. Parent benchmark: STRUCT-FINGERPRINT-20261021
(meta_acdc.models.train_fingerprint). Historical v2 used the following recipe. The shared batch assembler and
instance seed isolation were corrected on 2026-10-23; old artifacts remain archived:

- batch 16 -> 64
- projection head without the trailing LayerNorm (unit-sphere normalization
  before loss and evaluation, standard SupCon practice)
- temperature 0.1 -> 0.07
- emb dim 64 -> 128 (head 385->256->128)
- 60 -> 120 epochs, cosine decay (warmup 5, 3e-4 -> 3e-6)
- family-hard-negative weighting: same-family different-receptor negatives get
  weight (1 + lambda), lambda = 1.0, in the SupCon denominator

Everything else stays identical to v1: EGNN 6x128, feature schema 26/16,
invariant variant (coordinate updates off), family-grouped 5-fold assignment
(seed 20261021, v1 fold_assignment.json), Task A/B/C metric definitions,
dropout 0.1, grad clip 10. Deterministic algorithms are required (unsupported
CUDA operations fail); reproducibility must still be checked on the actual device.

Each (fold, seed) instance saves: weights, embeddings of all 331 graphs
(unit-normalized), the full similarity matrix, per-query retrieval rankings,
and its metrics fragment. These artifacts feed scripts/rashomon_fingerprint.py.

Usage:
    .venv/bin/python -m meta_acdc.models.train_fingerprint_v2 \
        --out-dir docs/fingerprint_retrieval_20261021/followup \
        --epochs 120 --batch-size 64 --seeds 20261021 1101 4242
"""
from __future__ import annotations

import argparse
import hashlib
import os
import json
import math
import random
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

from meta_acdc.models.egnn import EGNN
from meta_acdc.models.train_fingerprint import (
    auroc,
    batch_graphs,
    evaluate,
    group_folds,
    normalize_node_features,
)

SEED = 20261021
TEMPERATURE = 0.07
EMB_DIM = 128
EXECUTION_PROTOCOL = "isolated-fold-seed-and-offset-edges-v1"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_seed(seed: int):
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.use_deterministic_algorithms(True)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


class EmbedderV2(nn.Module):
    """EGNN (6x128, invariant) + projection head 385->256->128, no LayerNorm."""

    def __init__(self, node_dim: int = 26, edge_dim: int = 16,
                 emb_dim: int = EMB_DIM):
        super().__init__()
        self.egnn = EGNN(node_dim=node_dim, edge_dim=edge_dim, depth=6,
                         hidden=128, update_coordinates=False)
        pooled_dim = 3 * 128 + 1
        self.proj = nn.Sequential(
            nn.Linear(pooled_dim, 256), nn.SiLU(),
            nn.Linear(256, emb_dim))

    def forward(self, h, x, edge_index, edge_attr, batch, peptide_mask):
        pooled = self.egnn(h, x, edge_index, edge_attr, batch,
                           embed_only=True, peptide_mask=peptide_mask)
        return self.proj(pooled)


def family_neg_mask(receptor_ids: list[str], family_ids: list[str],
                    n: int) -> torch.Tensor:
    """(n, n) bool mask: (i, j) true when j is a same-family different-receptor
    negative for anchor i."""
    mask = torch.zeros(n, n, dtype=torch.bool)
    for i in range(n):
        for j in range(n):
            if i != j and receptor_ids[i] != receptor_ids[j] \
                    and family_ids[i] == family_ids[j]:
                mask[i, j] = True
    return mask


def supcon_loss_v2(embs: torch.Tensor, receptor_ids: list[str],
                   family_ids: list[str],
                   temperature: float = TEMPERATURE,
                   family_neg_weight: float = 1.0) -> torch.Tensor:
    """SupCon with family-hard-negative weighting.

    Anchors without any same-receptor partner in the batch contribute 0.
    """
    sim = embs @ embs.t() / temperature
    n = embs.shape[0]
    device = embs.device
    same = torch.zeros(n, n, device=device)
    for i in range(n):
        for j in range(n):
            if i != j and receptor_ids[i] == receptor_ids[j]:
                same[i, j] = 1.0
    hard = family_neg_mask(receptor_ids, family_ids, n).to(device)
    has_pos = same.sum(dim=1) > 0
    self_mask = torch.eye(n, device=device).bool()
    # denominator: weighted logsumexp over all j != i;
    # positives keep weight 1, hard negatives weight (1 + lambda).
    weights = torch.ones(n, n, device=device)
    weights[hard] = 1.0 + family_neg_weight
    logw = torch.zeros(n, n, device=device)
    logw[hard] = math.log(1.0 + family_neg_weight)
    denom = (sim + logw).masked_fill(self_mask, float("-inf")).logsumexp(dim=1)
    pos = sim.masked_fill(~(same > 0.5), float("-inf")).logsumexp(dim=1)
    per = -(pos - denom)
    per = torch.where(has_pos, per, torch.zeros_like(per))
    n_valid = has_pos.sum().clamp(min=1)
    return per.sum() / n_valid


def embed_all(graphs, records, device, model):
    """Embed every graph (single-graph batches), unit-normalize."""
    model.eval()
    embs = []
    with torch.no_grad():
        for i in range(len(graphs)):
            h, x, ei, ea, bt, pm, _ = batch_graphs(
                graphs, records, [i], device, model)
            embs.append(model(h, x, ei, ea, bt, pm).squeeze(0).cpu())
    embs = torch.stack(embs)
    return embs / embs.norm(dim=1, keepdim=True).clamp(min=1e-8)


def per_query_rankings(embs: torch.Tensor, test_idx: list[int]) -> dict:
    """Per-query gallery rankings (similarity-sorted indices + scores)."""
    sim = embs @ embs.t()
    n = embs.shape[0]
    out = {}
    for q in test_idx:
        others = [j for j in range(n) if j != q]
        order = sorted(others, key=lambda j: -float(sim[q, j]))
        out[str(q)] = {
            "query_graph": q,
            "gallery_order": order,
            "scores": [float(sim[q, j]) for j in order],
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", type=Path,
                    default=Path("data/processed/fingerprint_dataset.pt"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("docs/fingerprint_retrieval_20261021/followup"))
    ap.add_argument("--epochs", type=int, default=120)
    ap.add_argument("--warmup-epochs", type=int, default=5)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--temperature", type=float, default=TEMPERATURE)
    ap.add_argument("--emb-dim", type=int, default=EMB_DIM)
    ap.add_argument("--family-neg-weight", type=float, default=1.0)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--lr-min", type=float, default=3e-6)
    ap.add_argument("--wd", type=float, default=1e-4)
    ap.add_argument("--grad-clip", type=float, default=10.0)
    ap.add_argument("--seeds", type=int, nargs="+", default=[SEED])
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--device", default=None)
    ap.add_argument("--only-fold", type=int, default=None)
    ap.add_argument("--resume", action="store_true",
                    help="skip (fold, seed) instances whose full artifact set "
                         "already exists (weights + embeddings + similarity + "
                         "rankings + metrics); training math is unchanged")
    ap.add_argument("--max-epochs-probe", type=int, default=None,
                    help="internal: stop after N epochs (timing probe only)")
    args = ap.parse_args(argv)
    if (args.folds < 2 or args.epochs < 1 or args.batch_size < 1
            or len(set(args.seeds)) != len(args.seeds)
            or (args.only_fold is not None and not 0 <= args.only_fold < args.folds)):
        ap.error("invalid folds, epochs, batch size, duplicate seeds or fold index")

    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device: {device}", flush=True)

    data = torch.load(args.dataset, map_location="cpu", weights_only=False)
    graphs, records = data["graphs"], data["records"]
    families = [r["family_id"] for r in records]
    # same family-grouped fold assignment as v1 (seed fixed in group_folds)
    if len(set(families)) < args.folds:
        raise ValueError("fewer families than folds; empty test partitions are forbidden")
    folds = group_folds(families, args.folds)
    family_to_fold = {f: i for i, fs in enumerate(folds) for f in fs}
    fold_of = [family_to_fold[f] for f in families]

    run_dir = args.out_dir / "model_run_v2"

    config = {
        "task_id": "FOLLOWUP-AUROC-RASHOMON-20261021",
        "parent_task": "STRUCT-FINGERPRINT-20261021",
        "n_graphs": len(graphs),
        "execution_protocol": EXECUTION_PROTOCOL,
        "dataset_sha256": file_sha256(args.dataset),
        "probe_only": args.max_epochs_probe is not None,
        "source_sha256": {name: file_sha256(Path(__file__).with_name(name))
                          for name in ("train_fingerprint_v2.py", "train_fingerprint.py",
                                       "egnn.py", "train_egnn.py")},
        "torch_version": torch.__version__,
        "recipe_v2": {
            "batch_size": args.batch_size,
            "epochs": args.epochs,
            "warmup_epochs": args.warmup_epochs,
            "temperature": args.temperature,
            "emb_dim": args.emb_dim,
            "family_neg_weight": args.family_neg_weight,
            "lr_base": args.lr,
            "lr_min": args.lr_min,
            "wd": args.wd,
            "grad_clip": args.grad_clip,
            "egnn": "6 layers hidden 128, node/edge 26/16, invariant",
            "head": "385->256->128, no LayerNorm; unit-sphere before loss/eval",
            "fold_assignment": "v1 family-grouped, seed 20261021",
        },
        "seeds": args.seeds,
        "device": device,
        "binding_labels": False,
    }
    config_path = args.out_dir / "config.json"
    if config_path.exists():
        previous = json.loads(config_path.read_text())
        compatible = config
        old_compatible = previous
        if not args.resume or old_compatible != compatible:
            raise ValueError("refusing legacy/incompatible or non-resume output; use a fresh directory")
    else:
        if args.out_dir.exists() and any(args.out_dir.iterdir()):
            raise ValueError("output directory is nonempty without a compatible config")
        args.out_dir.mkdir(parents=True, exist_ok=True)
        config_path.write_text(json.dumps(config, indent=2, allow_nan=False) + "\n")
    run_dir.mkdir(parents=True, exist_ok=True)

    fold_list = range(args.folds) if args.only_fold is None \
        else [args.only_fold]

    for seed in args.seeds:
        for fold in fold_list:
            tag = f"fold{fold}_seed{seed}"
            artifacts = (run_dir / f"{tag}.pt", run_dir / f"emb_{tag}.npy",
                         run_dir / f"sim_{tag}.npy", args.out_dir / f"rankings_{tag}.json")
            metric_path = args.out_dir / f"metrics_{tag}.json"
            complete = all(p.exists() for p in (*artifacts, metric_path))
            if args.resume and complete:
                saved = json.loads(metric_path.read_text())
                hashes = {str(p.relative_to(args.out_dir)): file_sha256(p) for p in artifacts}
                if (saved.get("fold") != fold or saved.get("seed") != seed
                        or saved.get("execution_protocol") != EXECUTION_PROTOCOL
                        or saved.get("artifact_sha256") != hashes):
                    raise ValueError("resume artifacts are legacy, incomplete or modified; use a fresh directory")
                print(f"seed {seed} fold {fold}: verified complete, skipping (--resume)", flush=True)
                continue
            if any(p.exists() for p in (*artifacts, metric_path)):
                raise ValueError("refusing to overwrite existing/partial instance artifacts")
            # Each instance must be independent of earlier folds and skipped instances.
            set_seed(seed)
            test_idx = [i for i, f in enumerate(fold_of) if f == fold]
            train_idx = [i for i in range(len(graphs)) if i not in test_idx]
            print(f"seed {seed} fold {fold}: train={len(train_idx)} "
                  f"test={len(test_idx)}", flush=True)
            model = EmbedderV2(26, 16, emb_dim=args.emb_dim).to(device)
            opt = torch.optim.Adam(model.parameters(), lr=args.lr,
                                   weight_decay=args.wd)
            model.train()
            epochs_done = args.epochs
            if args.max_epochs_probe is not None:
                epochs_done = min(args.epochs, args.max_epochs_probe)
            for epoch in range(epochs_done):
                # cosine schedule with warmup
                if epoch < args.warmup_epochs:
                    lr = args.lr * (epoch + 1) / args.warmup_epochs
                else:
                    t = (epoch - args.warmup_epochs) / max(
                        1, args.epochs - args.warmup_epochs)
                    lr = args.lr_min + 0.5 * (args.lr - args.lr_min) \
                        * (1.0 + math.cos(math.pi * t))
                for g in opt.param_groups:
                    g["lr"] = lr
                order = train_idx[:]
                random.shuffle(order)
                total, nb = 0.0, 0
                for s in range(0, len(order), args.batch_size):
                    bidx = order[s:s + args.batch_size]
                    h, x, ei, ea, bt, pm, rids = batch_graphs(
                        graphs, records, bidx, device, model)
                    fams_b = [records[i]["family_id"] for i in bidx]
                    embs = model(h, x, ei, ea, bt, pm)
                    embs = embs / embs.norm(dim=1, keepdim=True).clamp(min=1e-8)
                    loss = supcon_loss_v2(embs, rids, fams_b,
                                          args.temperature,
                                          args.family_neg_weight)
                    if not torch.isfinite(loss):
                        raise RuntimeError(
                            f"seed {seed} fold {fold} epoch {epoch}: "
                            f"non-finite SupCon loss; training diverged")
                    opt.zero_grad()
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(
                        model.parameters(), args.grad_clip)
                    opt.step()
                    total += float(loss.item())
                    nb += 1
                if epoch % 15 == 0 or epoch == epochs_done - 1:
                    print(f"  seed {seed} fold {fold} epoch {epoch}: "
                          f"loss {total / max(nb, 1):.4f} lr {lr:.2e}",
                          flush=True)
            if args.max_epochs_probe is not None:
                print("probe mode: stopping without artifacts", flush=True)
                return 0
            torch.save({"model": model.state_dict(), "fold": fold,
                        "seed": seed, "config": config},
                       run_dir / f"{tag}.pt")
            embs = embed_all(graphs, records, device, model)
            sim = embs @ embs.t()
            sim_np = sim.numpy().astype(np.float64)
            med_pair = float(np.median(sim_np[np.triu_indices(len(graphs), 1)]))
            res = evaluate(embs, records, test_idx, None)
            res["fold"] = fold
            res["seed"] = seed
            res["n_train"] = len(train_idx)
            res["n_test"] = len(test_idx)
            res["median_pairwise_cosine"] = med_pair
            frag = args.out_dir / f"metrics_{tag}.json"
            np.save(run_dir / f"emb_{tag}.npy", embs.numpy().astype(np.float32))
            np.save(run_dir / f"sim_{tag}.npy", sim_np)
            rankings = per_query_rankings(embs, test_idx)
            (args.out_dir / f"rankings_{tag}.json").write_text(
                json.dumps(rankings) + "\n")
            res["execution_protocol"] = EXECUTION_PROTOCOL
            res["artifact_sha256"] = {str(p.relative_to(args.out_dir)): file_sha256(p)
                                      for p in artifacts}
            # Completion marker is written last; resume verifies all four bound artifacts.
            frag.write_text(json.dumps(res, indent=2, allow_nan=False) + "\n")
            print(json.dumps(
                {k: v for k, v in res.items()
                 if k in ("fold", "seed", "retrieval", "family", "domain",
                          "median_pairwise_cosine")}, indent=2), flush=True)

    # aggregate across completed (fold, seed) fragments into followup_results.json
    frags = sorted(args.out_dir.glob("metrics_fold*_seed*.json"))
    fold_results = []
    for fp in frags:
        d = json.loads(fp.read_text())
        if d.get("seed") == args.seeds[0]:
            fold_results.append(d)
    if frags:
        fold_results.sort(key=lambda d: d["fold"])
        means = {
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
            "domain_accuracy": float(np.nanmean(
                [r["domain"]["accuracy"] for r in fold_results])),
            "median_pairwise_cosine": float(np.nanmean(
                [r["median_pairwise_cosine"] for r in fold_results])),
        }
        agg = dict(config)
        agg["fold_results"] = fold_results
        agg["means"] = means
        out = args.out_dir / "followup_results.json"
        out.write_text(json.dumps(agg, indent=2, allow_nan=False) + "\n")
        print(json.dumps(means, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
