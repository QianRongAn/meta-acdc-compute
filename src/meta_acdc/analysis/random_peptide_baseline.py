"""Random-peptide score distribution machinery (review item #20).

No valid structure-binding model is currently released, so this tool
REFUSES to run unless given an accepted release and a checkpoint whose
release fingerprint matches. When a future release passes the acceptance
gate, this produces the distribution of scores for construction-matched
random 9-mers on a chosen receptor scaffold — the reference needed to
interpret any absolute score.

What it is and is not:
- The random peptides are grafted onto the native peptide backbone with the
  SAME graft construction used for synthetic controls (make_decoy), so the
  comparison is construction-matched.
- The output is a DESCRIPTIVE score distribution only. Sigmoid outputs are
  not binding probabilities, and no clinical interpretation follows.

Usage:
    .venv/bin/python -m meta_acdc.analysis.random_peptide_baseline \
        --receptor-cif data/raw/af3_predictions/5brz_EVDPIGHLY_model_0.cif \
        --model data/processed/da_seed0.model.pt \
        --release data/releases/<accepted-release>/release_manifest.json \
        --n 1000 --seed 0 --out data/processed/random_peptide_baseline.tsv
"""
from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

from meta_acdc.data.release import load_release, sha256
from meta_acdc.models.dataset import make_decoy
from meta_acdc.models.egnn import EGNN
from meta_acdc.structure.artifact_qc import require_construction_screen
from meta_acdc.structure.cif import parse_cif
from meta_acdc.structure.graph import (
    PEPTIDE,
    build_interface_graph_from_residues,
    classify_chains,
)
from meta_acdc.structure.identity import chain_sequence
from meta_acdc.structure.quality import annotate_prediction
from meta_acdc.structure.score_predictions import score_graph

_AA = "ACDEFGHIKLMNPQRSTVWY"


def random_peptide(length: int, rng: random.Random) -> list:
    """Random peptide as (resname, resid, chain) residue shells.

    Only residue identity matters for grafting — geometry comes from the
    native backbone, so these carry no coordinates of their own.
    """
    return [type("_R", (), {"resname": rng.choice(_AA), "resid": i + 1,
                            "chain": "P", "atoms": [], "x": 0.0, "y": 0.0,
                            "z": 0.0, "plddt": None, "icode": "",
                            "role": PEPTIDE, "chain_kind": "peptide"})()
            for i in range(length)]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptor-cif", type=Path, required=True,
                    help="native complex CIF used as the receptor scaffold; "
                         "must pass the construction screen")
    ap.add_argument("--model", type=Path, required=True,
                    help="checkpoint whose release_sha256 must match --release")
    ap.add_argument("--release", type=Path, required=True,
                    help="accepted release manifest (release_manifest.json)")
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args(argv)
    if args.n <= 0:
        raise ValueError("--n must be positive")
    if args.out.exists():
        raise ValueError("refusing to overwrite an existing baseline table")

    import json
    import torch

    manifest, _ = load_release(args.release)
    records = {x["complex_id"]: x for x in json.loads(
        (args.release.parent / "structure_records.json").read_text())}

    ckpt = torch.load(args.model, map_location="cpu", weights_only=False)
    if ckpt.get("release_sha256") != sha256(args.release):
        raise ValueError(
            "legacy or mismatched checkpoint; a random-peptide baseline is "
            "only meaningful on an accepted release-bound model")
    net = EGNN(node_dim=ckpt["n_node_dim"], edge_dim=ckpt["n_edge_dim"],
               depth=6, hidden=128,
               update_coordinates=ckpt.get("update_coordinates", True))
    net.load_state_dict(ckpt["model"])
    net.eval()
    net.n_node_dim = ckpt["n_node_dim"]
    net.n_edge_dim = ckpt["n_edge_dim"]
    mask_plddt = bool(ckpt.get("mask_plddt", False))

    # Identify the scaffold against the release records.
    job = args.receptor_cif.stem.rsplit("_model_", 1)[0]
    pdb, _, expected_peptide = job.partition("_")
    residues = parse_cif(args.receptor_cif, confidence_kind="af3_plddt")
    candidates = [r for r in records.values()
                  if r.get("peptide") == expected_peptide]
    if not candidates:
        raise ValueError(f"no release record for peptide {expected_peptide}")
    ref = sorted(candidates, key=lambda x: x["complex_id"])[0]
    annotate_prediction(residues, ref, expected_peptide)
    require_construction_screen(residues)
    roles = classify_chains(residues)
    native_pep = [r for r in residues if roles[r.chain] == PEPTIDE]
    native_seq = chain_sequence(native_pep)
    g = build_interface_graph_from_residues(
        residues, electrostatics=ckpt["n_node_dim"] > 26)
    native_score = score_graph(net, g, args.device, mask_plddt=mask_plddt)

    rng = random.Random(args.seed)
    rows = [{"peptide": native_seq, "kind": "native",
             "score": native_score, "n_nodes": g.n_nodes}]
    rejected = 0
    for _ in range(args.n):
        pep = random_peptide(len(native_pep), rng)
        seq = "".join(r.resname for r in pep)
        if seq == native_seq:
            rejected += 1
            continue
        try:
            decoy = make_decoy(residues, pep)
            require_construction_screen(decoy)
            dg = build_interface_graph_from_residues(
                decoy, electrostatics=ckpt["n_node_dim"] > 26)
            s = score_graph(net, dg, args.device, mask_plddt=mask_plddt)
        except ValueError:
            rejected += 1
            continue
        rows.append({"peptide": seq, "kind": "random-graft",
                     "score": s, "n_nodes": dg.n_nodes})

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["peptide", "kind", "score",
                                           "n_nodes"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    scores = sorted(r["score"] for r in rows if r["kind"] == "random-graft")
    if scores:
        q = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
        import numpy as np
        quantiles = {f"q{int(p * 100)}": float(np.quantile(scores, p))
                     for p in q}
    else:
        quantiles = {}
    summary = {
        "receptor_job": job,
        "receptor_family": ref.get("family_id"),
        "native_score": native_score,
        "n_scored_random": len(scores),
        "n_rejected": rejected,
        "quantiles": quantiles,
        "notice": "descriptive score distribution on construction-matched "
                  "random grafts only; sigmoid is not a binding probability; "
                  "no clinical or binding interpretation follows",
    }
    sp = args.out.with_suffix(".summary.json")
    sp.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
