"""Temperature calibration for EGNN compatibility scores.

The raw sigmoid scores are not calibrated binding probabilities (BCE training
only constrains within-group ordering). Temperature scaling fits a single
scalar T on a held-out labeled set (crystal decoy val split) so that
sigmoid(logit / T) is calibrated ON THAT TASK — a prerequisite for any
absolute score threshold in a safety screen. AUROC is invariant to T;
calibration affects the probability scale and expected calibration error (ECE).

WHAT THIS DOES NOT GIVE YOU (AUDIT.md #19): the calibration labels are the
decoy labels, i.e. "real crystal structure vs artificially constructed
decoy". Real off-target screening has a different class ratio, different
negative sources (peptides nobody tested) and a different structure
distribution (predicted, not crystallised). A calibrated score therefore may
not be read as "probability of binding", and certainly not as "probability
that the peptide is dangerous". Calibration must be re-estimated on a task
that matches its intended use, and T must not be fitted and reported on the
same set.

Usage:
    .venv/bin/python src/meta_acdc/models/calibrate.py \
        --model data/processed/da_seed0.model.pt \
        --data data/processed/egnn_dataset.pt \
        --out data/processed/da_seed0.temperature.json
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score

from meta_acdc.models.egnn import EGNN
from meta_acdc.models.train_egnn import collate


def mask_plddt(graphs):
    out = []
    for nf, coords, ei, ef in graphs:
        nf = np.array(nf, dtype=np.float32).copy()
        nf[:, 3] = 0.0
        out.append((nf.tolist(), coords, ei, ef))
    return out


def ece(probs: np.ndarray, labels: np.ndarray, n_bins: int = 10) -> float:
    """Expected calibration error with closed bins.

    Bin i covers (lo, hi] except the first, which is [0, hi]. The previous
    `probs > lo` left every prediction of exactly 0.0 outside all bins, so a
    model that always predicts 0 against all-positive labels scored ECE = 0.
    See AUDIT.md #18.
    """
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    total = 0.0
    for i, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
        lower_ok = probs >= lo if i == 0 else probs > lo
        m = lower_ok & (probs <= hi)
        if m.sum() == 0:
            continue
        conf = probs[m].mean()
        acc = labels[m].mean()
        total += m.mean() * abs(conf - acc)
    return float(total)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--data", type=Path,
                    default=None)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--device",
                    default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--release", type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():raise ValueError("refusing to overwrite calibration output")
    from meta_acdc.data.release import load_release, frozen_indices, sha256
    manifest, approved = load_release(args.release, data_path=args.data)
    args.data = approved

    data = torch.load(args.data, weights_only=False)
    graphs = list(data["graphs"])
    labels = np.asarray(data["labels"])
    pdbs = list(data["pdbs"])

    frozen = frozen_indices(manifest, data)
    va_idx, te_idx = frozen['calibration'], frozen['test']

    ck = torch.load(args.model, weights_only=False)
    if ck.get("release_sha256") != sha256(args.release):
        raise ValueError("checkpoint was not trained on this accepted release")
    if ck["mask_plddt"]:graphs=mask_plddt(graphs)
    net = EGNN(node_dim=ck["n_node_dim"], edge_dim=ck["n_edge_dim"],
               depth=6, hidden=128).to(args.device)
    net.load_state_dict(ck["model"])
    net.eval()

    def logits_for(idx):
        b = {k: v.to(args.device) for k, v in
             collate([graphs[i] for i in idx], ck["n_node_dim"], ck["n_edge_dim"]).items()}
        with torch.no_grad():
            return net(**b).squeeze(-1).cpu(), \
                torch.tensor(labels[idx], dtype=torch.float32)

    logits_cal, y_cal = logits_for(va_idx)
    logits_test, y_test = logits_for(te_idx)

    p0_test = torch.sigmoid(logits_test).numpy()
    print(f"cal n={len(y_cal)} pos={int(y_cal.sum())} | "
          f"held-out n={len(y_test)} pos={int(y_test.sum())}")
    print(f"before (held-out): AUROC={roc_auc_score(y_test, p0_test):.3f} "
          f"ECE={ece(p0_test, y_test.numpy()):.3f}")

    # fit temperature by minimizing BCE on logits/T — on the CALIBRATION split
    logT = torch.zeros(1, requires_grad=True)
    opt = torch.optim.LBFGS([logT], lr=0.1, max_iter=100)
    loss_fn = nn.BCEWithLogitsLoss()

    def closure():
        opt.zero_grad()
        loss = loss_fn(logits_cal / torch.exp(logT), y_cal)
        loss.backward()
        return loss

    opt.step(closure)
    T = float(torch.exp(logT).item())
    p1_test = torch.sigmoid(logits_test / T).detach().numpy()
    p1_cal = torch.sigmoid(logits_cal / T).detach().numpy()
    print(f"after : T={T:.3f} (fit on cal split)")
    print(f"  held-out : AUROC={roc_auc_score(y_test, p1_test):.3f} "
          f"ECE={ece(p1_test, y_test.numpy()):.3f}")
    print(f"  in-sample cal split (NOT independent): "
          f"ECE={ece(p1_cal, y_cal.numpy()):.3f}")

    with open(args.out, "w") as fh:
        json.dump({"release_sha256": sha256(args.release), "task": manifest["task"], "model": str(args.model), "temperature": T,
                   "n_cal": len(y_cal), "n_heldout": len(y_test),
                   "ece_heldout_before": ece(p0_test, y_test.numpy()),
                   "ece_heldout_after": ece(p1_test, y_test.numpy()),
                   "ece_insample_after": ece(p1_cal, y_cal.numpy())},
                  fh, indent=2)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())