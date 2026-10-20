"""Explainability: which residue pairs drive a compatibility score?

Gradient-based attribution on the interface graph. The score is
differentiated with respect to the *side-chain contact* part of each edge
feature (the chemistry channels), giving a per-edge importance. Edges are
then aggregated per residue pair and reported with chain/residue identities
and roles (peptide / TCR / MHC). This turns a black-box score into a
human-readable "these contacts make the peptide look like a binder" report —
a first step toward mechanistic hypotheses for cross-reactivity.

Since the delivery review: chain identity is no longer guessed from chain
letters or lengths. The complex must be explicitly annotated from a release
record (`--reference-json`, the structure_records entry) and the checkpoint
must carry a release fingerprint. Legacy checkpoints and unannotated CIFs
are refused.

Usage:
    .venv/bin/python src/meta_acdc/structure/attribute.py \
        --cif data/raw/af3_predictions/5brz_ESDPIVAQY_model_0.cif \
        --model data/processed/da_seed0.model.pt --mask-plddt --top 15 \
        --reference-json <structure_record>.json --peptide ESDPIVAQY
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import torch

from meta_acdc.models.egnn import EGNN
from meta_acdc.structure.cif import parse_cif
from meta_acdc.structure.graph import (B2M, MHC, PEPTIDE, TCR,
                                       build_interface_graph_from_residues,
                                       classify_chains)
from meta_acdc.structure.quality import annotate_prediction
from meta_acdc.structure.score_predictions import graph_to_tensors

ROLE_NAME = {TCR: "TCR", PEPTIDE: "PEP", MHC: "MHC", B2M: "B2M"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", type=Path, required=True)
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--reference-json", type=Path, required=True,
                    help="structure_records entry with alpha/beta/"
                         "mhc_sequence/b2m_sequence; chain identity is "
                         "assigned from this, never guessed from chain "
                         "letters or lengths")
    ap.add_argument("--peptide", type=str, required=True,
                    help="expected peptide sequence of the complex")
    ap.add_argument("--mask-plddt", action="store_true")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--device",
                    default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    reference = json.loads(args.reference_json.read_text())
    for key in ("alpha", "beta", "mhc_sequence", "b2m_sequence"):
        if not reference.get(key):
            raise ValueError(
                f"--reference-json lacks {key}; provide a full "
                "structure_records entry")

    ck = torch.load(args.model, weights_only=False)
    if not ck.get("release_sha256"):
        raise ValueError(
            "legacy or mismatched checkpoint: attribution requires a "
            "release-bound model (see docs/delivery_audit_20261018/README.md); "
            "no structural retraining has been released")
    net = EGNN(node_dim=ck["n_node_dim"], edge_dim=ck["n_edge_dim"],
               depth=6, hidden=128).to(args.device)
    net.load_state_dict(ck["model"])
    net.eval()

    residues = parse_cif(args.cif, confidence_kind="af3_plddt")
    annotate_prediction(residues, reference, args.peptide)
    roles = classify_chains(residues)
    g = build_interface_graph_from_residues(residues)
    info = g.residue_info  # (chain, resname, resid) per node

    t = {k: v.to(args.device) for k, v in graph_to_tensors(g).items()}
    h = t["h"]
    if args.mask_plddt:
        h = h.clone()
        h[:, 3] = 0.0
    h = (h - h.mean(0, keepdim=True)) / (h.std(0, keepdim=True) + 1e-6)

    edge_attr = t["edge_attr"].clone().requires_grad_(True)
    logit = net(h, t["x"], t["edge_index"], edge_attr, t["batch"]).squeeze(-1)
    score = float(torch.sigmoid(logit).item())
    net.zero_grad()
    logit.backward()
    grad = edge_attr.grad  # (E, edge_dim)

    # importance = |grad| summed over the 4 contact-chemistry channels (12:16)
    contact_importance = grad[:, 12:].abs().sum(dim=-1)
    ei = t["edge_index"].cpu()
    imp = contact_importance.detach().cpu()

    # aggregate symmetric edges into residue pairs
    pair: dict[tuple[int, int], float] = {}
    for e in range(ei.shape[1]):
        i, j = int(ei[0, e]), int(ei[1, e])
        key = (min(i, j), max(i, j))
        pair[key] = pair.get(key, 0.0) + float(imp[e])

    print(f"{args.cif.name}  score={score:.3f}  "
          f"nodes={g.n_nodes} edges={ei.shape[1] // 2}")

    def fmt(i: int) -> str:
        c, r, n = info[i]
        return f"{c}:{r}{n}"

    def role(i: int) -> str:
        return ROLE_NAME.get(roles.get(info[i][0]), "?")

    def show(rows, title):
        print(f"\n{title}")
        print(f"  {'residue A':22s} {'residue B':22s} {'roles':10s} {'weight':>8s}")
        for (i, j), w in rows:
            print(f"  {fmt(i):22s} {fmt(j):22s} "
                  f"{role(i) + '-' + role(j):10s} {w:8.4f}")

    ranked = sorted(pair.items(), key=lambda kv: -kv[1])
    show(ranked[:args.top], f"top {args.top} contacts (all role pairs)")
    # mechanistically relevant: peptide-involving contacts only
    pep = [(k, w) for k, w in ranked
           if PEPTIDE in (roles.get(info[k[0]][0]), roles.get(info[k[1]][0]))]
    show(pep[:args.top], "top peptide-involving contacts (TCR-PEP / PEP-MHC)")
    return 0


if __name__ == "__main__":
    sys.exit(main())