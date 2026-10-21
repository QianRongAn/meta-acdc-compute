"""Diagnose the AUROC gap using the committed fold models (read-only analysis
of existing artifacts; no new training). Breaks the same/diff pair AUROC into:
  - same receptor vs different receptor (same family)   [hard confusables]
  - same receptor vs different receptor (cross family)
and times the embedding pass per model.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from meta_acdc.models.train_fingerprint import (  # noqa: E402
    Embedder, batch_graphs, evaluate, group_folds, auroc)

ROOT = Path(__file__).resolve().parents[1]
DS = ROOT / "data/processed/fingerprint_dataset.pt"
RUN = ROOT / "docs/fingerprint_retrieval_20261021/model_run"


def embed_one(model, graphs, records, device):
    embs = []
    model.eval()
    with torch.no_grad():
        for i in range(len(graphs)):
            h, x, ei, ea, bt, pm, _ = batch_graphs(graphs, records, [i], device, model)
            embs.append(model(h, x, ei, ea, bt, pm).squeeze(0).cpu())
    embs = torch.stack(embs)
    return embs / embs.norm(dim=1, keepdim=True).clamp(min=1e-8)


def main():
    device = "cuda" if torch.cuda.is_available() else "cpu"
    data = torch.load(DS, map_location="cpu", weights_only=False)
    graphs, records = data["graphs"], data["records"]
    fams = [r["family_id"] for r in records]
    folds = group_folds(fams, 5)
    f2f = {f: i for i, fs in enumerate(folds) for f in fs}
    fold_of = [f2f[f] for f in fams]

    out = {"per_fold": {}, "timing_s": {}}
    for fold in range(5):
        test_idx = [i for i, f in enumerate(fold_of) if f == fold]
        model = Embedder(26, 16, update_coordinates=False).to(device)
        sd = torch.load(RUN / f"fold{fold}.pt", map_location="cpu", weights_only=False)
        model.load_state_dict(sd["model"])
        t0 = time.time()
        embs = embed_one(model, graphs, records, device)
        dt = time.time() - t0
        out["timing_s"][f"fold{fold}"] = round(dt, 1)
        print(f"fold {fold}: embedding pass {dt:.1f}s", flush=True)

        sim = embs @ embs.t()
        n = len(records)
        # pairs anchored at test queries (same pair scheme as evaluate())
        sc_all, lb_all = [], []
        sc_samefam, lb_samefam = [], []
        sc_crossfam, lb_crossfam = [], []
        sc_samerec = []
        for q in test_idx:
            for j in range(n):
                if j == q:
                    continue
                s = float(sim[q, j])
                lab = records[j]["receptor_id"] == records[q]["receptor_id"]
                sc_all.append(s)
                lb_all.append(lab)
                if records[j]["family_id"] == records[q]["family_id"]:
                    sc_samefam.append(s)
                    lb_samefam.append(lab)
                else:
                    sc_crossfam.append(s)
                    lb_crossfam.append(lab)
                if lab:
                    sc_samerec.append(s)
        def au(x, y):
            return auroc(np.array(x), np.array(y))
        hard_au = au(sc_samefam, lb_samefam)
        cross_au = au(sc_crossfam, lb_crossfam)
        n_samefam_pos = sum(lb_samefam)
        n_samefam_neg = len(lb_samefam) - n_samefam_pos
        n_cross_neg = len(lb_crossfam)
        out["per_fold"][f"fold{fold}"] = {
            "auroc_all": au(sc_all, lb_all),
            "auroc_same_family_pairs": hard_au,
            "auroc_cross_family_pairs": cross_au,
            "n_samefam_pos": n_samefam_pos,
            "n_samefam_neg": n_samefam_neg,
            "n_crossfam_neg": n_cross_neg,
            "sim_same_receptor_mean": float(np.mean(sc_samerec)),
            "sim_diff_receptor_samefam_mean": float(np.mean([s for s, l in zip(sc_samefam, lb_samefam) if not l])),
        }
        print(json.dumps(out["per_fold"][f"fold{fold}"], indent=2), flush=True)

    path = ROOT / "docs/fingerprint_retrieval_20261021/auroc_gap_diagnosis.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", path)


if __name__ == "__main__":
    main()
