"""Domain-adaptation training: crystals + TCRmodel2-predicted structures.

Fixes the Rashomon-effect instability: the model must see "predicted-style"
inputs during training. TCRmodel2 PDBs (data/raw/tcrmodel/) are added as
extra POSITIVE examples alongside the crystal natives, with the SAME decoy
protocol (displaced + graft) applied to them.

pLDDT is zeroed in training AND scoring (the AF3 B-factor confound).

After training 3 seeds, reports pairwise Spearman rank correlation on the
held-out decoy set across seeds (the stability metric).

Usage:
    .venv/bin/python src/meta_acdc/models/train_domain_adapt.py \
        --base data/processed/egnn_dataset.pt \
        --tcrmodel data/raw/tcrmodel \
        --epochs 150
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score

from meta_acdc.models.dataset import (
    _peptide_chains,
    make_decoy,
    make_displaced_peptide_decoy,
)
from meta_acdc.models.egnn import EGNN
from meta_acdc.models.train_egnn import collate
from meta_acdc.structure.cif import parse_cif
from meta_acdc.structure.graph import (MHC, PEPTIDE, TCR,
                                       build_interface_graph_from_residues,
                                       classify_chains, parse_pdb)

SEEDS = [0, 1, 2]


def mask_plddt(gs):
    out = []
    for nf, coords, ei, ef in gs:
        nf = np.array(nf, dtype=np.float32).copy()
        nf[:, 3] = 0.0
        out.append((nf.tolist(), coords, ei, ef))
    return out


def _roles_ok(roles: dict) -> bool:
    return (sum(1 for v in roles.values() if v == TCR) >= 1
            and sum(1 for v in roles.values() if v == PEPTIDE) == 1
            and sum(1 for v in roles.values() if v == MHC) >= 1)


def add_complex_set(graphs, labels, pdbs, complexes: dict, tag: str):
    """Add predicted complexes as positives + displaced/graft decoys.

    Node-set pinned (forced_nodes) and graft donors sequence-deduplicated —
    the same leakage guard as build_dataset(). `complexes` maps stem ->
    residue list (from PDB or CIF, both feed the same graph builder).
    """
    rng = random.Random(43)

    def add(residues, label, name, forced_nodes=None):
        try:
            g = build_interface_graph_from_residues(residues,
                                                    forced_nodes=forced_nodes)
        except ValueError:
            return None
        if g.n_nodes < 20:
            return None
        graphs.append((g.node_features, g.node_coords,
                       g.edge_index, g.edge_features))
        labels.append(label)
        pdbs.append(name)
        return g

    native_keys = {}
    for stem, residues in complexes.items():
        g = add(residues, 1, f"{tag}{stem}")
        if g is not None:
            native_keys[stem] = set(g.node_keys)

    n_disp = n_graft = 0
    stems = list(complexes.keys())
    for stem, residues in complexes.items():
        if stem not in native_keys:
            continue
        add(make_displaced_peptide_decoy(residues), 0,
            f"{tag}{stem}-disp", forced_nodes=native_keys[stem])
        n_disp += 1
        peptides = _peptide_chains(residues)
        if not peptides:
            continue
        native_seq = "".join(r.resname[:1] for r in next(iter(peptides.values())))
        graft = None
        donor = stem
        for _ in range(20):
            if len(stems) < 2:
                break
            donor = rng.choice(stems)
            if donor == stem:
                continue
            donor_peps = _peptide_chains(complexes[donor])
            if not donor_peps:
                continue
            cand = rng.choice(list(donor_peps.values()))
            if "".join(r.resname[:1] for r in cand) != native_seq:
                graft = cand
                break
        if graft is not None:
            add(make_decoy(residues, graft), 0, f"{tag}{stem}<-{donor}",
                forced_nodes=native_keys[stem])
            n_graft += 1

    print(f"{tag}: {len(native_keys)} natives, {n_disp} displaced, "
          f"{n_graft} graft decoys", flush=True)
    return graphs, labels, pdbs


def add_tcrmodel_graphs(graphs, labels, pdbs, tcr_dir: Path):
    """TCRmodel2 PDBs as predicted-style positives (legacy route)."""
    complexes = {}
    for p in sorted(tcr_dir.glob("*.pdb")):
        try:
            residues = parse_pdb(p)
            roles = classify_chains(residues)
        except Exception:
            continue
        if _roles_ok(roles):
            complexes[p.stem] = residues
    return add_complex_set(graphs, labels, pdbs, complexes, "tcr-")


def add_af3_graphs(graphs, labels, pdbs, af3_dir: Path,
                   native_ids: set[str] | None = None):
    """AF3-predicted CIFs as predicted-style positives (current route).

    Only complexes whose job_id is in `native_ids` are used, so that
    cross-reactivity candidates are never mislabeled as cognate binders.
    """
    complexes = {}
    for cif in sorted(af3_dir.glob("*_model_0.cif")):
        job = cif.name.rsplit("_model_0", 1)[0]
        if native_ids is not None and job not in native_ids:
            continue
        try:
            residues = parse_cif(cif)
            roles = classify_chains(residues)
        except Exception:
            continue
        if _roles_ok(roles):
            complexes[job] = residues
    return add_complex_set(graphs, labels, pdbs, complexes, "af3-")


def pdb_group(name: str) -> str:
    return name.split("-disp")[0].split("<-")[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=Path,
                    default=Path("data/processed/egnn_dataset.pt"))
    ap.add_argument("--tcrmodel", type=Path, default=Path("data/raw/tcrmodel"),
                    help="legacy route: TCRmodel2 PDB positives")
    ap.add_argument("--af3", type=Path, default=Path("data/raw/af3_native"),
                    help="current route: AF3-predicted native CIFs "
                         "(style-matched positives)")
    ap.add_argument("--af3-manifest", type=Path, default=None,
                    help="optional: restrict AF3 positives to native job_ids "
                         "listed here (use only when --af3 mixes natives with "
                         "cross-reactivity candidates)")
    ap.add_argument("--epochs", type=int, default=150)
    ap.add_argument("--seeds", type=int, nargs="+", default=SEEDS,
                    help="random seeds to train (default 0 1 2); the stability "
                         "gate is a function of the ensemble size")
    ap.add_argument("--out-dir", type=Path, default=Path("data/processed"),
                    help="where da_seed{0,1,2}.model.pt are written")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()
    seeds = list(args.seeds)

    d = torch.load(args.base, weights_only=False)
    graphs = list(d["graphs"])
    labels = list(d["labels"])
    pdbs = list(d["pdbs"])

    native_ids = None
    if args.af3_manifest is not None and args.af3_manifest.exists():
        with open(args.af3_manifest) as fh:
            next(fh, None)  # header
            native_ids = {ln.split("\t")[0] for ln in fh if ln.strip()}
    if args.af3.is_dir() and any(args.af3.glob("*_model_0.cif")):
        print(f"adding AF3 predicted positives from {args.af3} "
              f"(native filter: {len(native_ids) if native_ids else 'off'})",
              flush=True)
        graphs, labels, pdbs = add_af3_graphs(
            graphs, labels, pdbs, args.af3, native_ids=native_ids)
    elif args.tcrmodel.is_dir():
        graphs, labels, pdbs = add_tcrmodel_graphs(
            graphs, labels, pdbs, args.tcrmodel)
    labels_t = torch.tensor(labels, dtype=torch.float32)
    graphs_masked = mask_plddt(graphs)
    print(f"DA dataset: {len(graphs)} graphs "
          f"({sum(labels)} pos / {len(labels) - sum(labels)} neg)", flush=True)

    groups: dict[str, list[int]] = {}
    for i, p in enumerate(pdbs):
        groups.setdefault(pdb_group(p), []).append(i)
    rng = random.Random(42)
    gids = list(groups.keys())
    rng.shuffle(gids)
    n_val = max(1, int(0.2 * len(gids)))
    val_g = set(gids[:n_val])
    tr = [i for g in gids if g not in val_g for i in groups[g]]
    va = [i for g in gids if g in val_g for i in groups[g]]

    device = args.device
    bs = 32
    all_probs = []
    for seed in seeds:
        torch.manual_seed(seed)
        net = EGNN(node_dim=26, edge_dim=16, depth=6, hidden=128).to(device)
        opt = torch.optim.Adam(net.parameters(), lr=3e-4, weight_decay=1e-4)
        loss_fn = nn.BCEWithLogitsLoss()
        rng2 = random.Random(seed)
        for epoch in range(args.epochs):
            net.train()
            order = tr[:]
            rng2.shuffle(order)
            for start in range(0, len(order), bs):
                idx = order[start:start + bs]
                b = {k: v.to(device) for k, v in
                     collate([graphs_masked[i] for i in idx]).items()}
                y = labels_t[idx].to(device)
                opt.zero_grad()
                loss = loss_fn(net(**b).squeeze(-1), y)
                loss.backward()
                opt.step()
        net.eval()
        vb = {k: v.to(device) for k, v in
              collate([graphs_masked[i] for i in va]).items()}
        with torch.no_grad():
            p = torch.sigmoid(net(**vb).squeeze(-1)).cpu().numpy()
        y = labels_t[va].numpy()
        auc = roc_auc_score(y, p)
        print(f"seed {seed}: val AUROC {auc:.3f}", flush=True)
        all_probs.append(p)
        torch.save({"model": net.state_dict(), "n_node_dim": 26,
                    "n_edge_dim": 16},
                   args.out_dir / f"da_seed{seed}.model.pt")

    # stability: pairwise Spearman of val scores across seeds
    from scipy.stats import spearmanr
    print("\npairwise Spearman (val set) across seeds:")
    for i in range(len(seeds)):
        for j in range(i + 1, len(seeds)):
            r, _ = spearmanr(all_probs[i], all_probs[j])
            print(f"  seed{i} vs seed{j}: r={r:.3f}")
    print("DA training complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
