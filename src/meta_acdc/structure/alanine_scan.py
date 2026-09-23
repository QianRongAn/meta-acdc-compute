"""In-silico alanine scan to validate attribution (and probe mimicry).

For each peptide position, mutate that residue to alanine (keeping backbone
and CB) and re-score the interface. A large drop identifies positions the
model depends on — a direct, mutation-based cross-check of the gradient
attribution and a mechanistic probe of cross-reactivity.

Usage:
    .venv/bin/python src/meta_acdc/structure/alanine_scan.py \
        --cif data/raw/af3_predictions/5brz_ESDPIVAQY_model_0.cif \
        --model data/processed/da_seed0.model.pt --mask-plddt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import torch

from meta_acdc.models.egnn import EGNN
from meta_acdc.structure.cif import parse_cif
from meta_acdc.structure.graph import (AA3TO1, PEPTIDE, Residue,
                                       build_interface_graph_from_residues,
                                       classify_chains)
from meta_acdc.structure.score_predictions import graph_to_tensors

KEEP_ATOMS = ("N", "CA", "C", "O", "CB")


def to_mutant(res: Residue, target: str = "ALA") -> Residue:
    """Replace the residue identity, keeping backbone + CB coordinates.

    'ALA' is a permissive control (small, non-clashing); 'TRP'/'ARG' are
    clash probes (bulky / charged) that reveal positions the model depends
    on for compatibility.
    """
    atoms = [a for a in res.atoms if a[0] in KEEP_ATOMS]
    return Residue(res.chain, target, res.resid, res.x, res.y, res.z,
                   res.plddt, res.icode, atoms)


# backwards-compatible alias
def to_alanine(res: Residue) -> Residue:
    return to_mutant(res, "ALA")


def score(net, residues, device, mask_plddt):
    g = build_interface_graph_from_residues(residues)
    t = {k: v.to(device) for k, v in graph_to_tensors(g).items()}
    h = t["h"]
    if mask_plddt:
        h = h.clone()
        h[:, 3] = 0.0
    h = (h - h.mean(0, keepdim=True)) / (h.std(0, keepdim=True) + 1e-6)
    with torch.no_grad():
        logit = net(h, t["x"], t["edge_index"], t["edge_attr"],
                    t["batch"]).squeeze(-1)
    return float(torch.sigmoid(logit).item())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", type=Path, required=True)
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--mask-plddt", action="store_true")
    ap.add_argument("--mutant", default="ALA",
                    help="mutant residue 3-letter code; ALA = permissive "
                         "control, TRP/ARG = clash probe")
    ap.add_argument("--device",
                    default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    residues = parse_cif(args.cif)
    roles = classify_chains(residues)
    pep_chain = next((c for c, r in roles.items() if r == PEPTIDE), None)
    if pep_chain is None:
        print("no peptide chain found", file=sys.stderr)
        return 1
    pep_positions = [i for i, r in enumerate(residues)
                     if r.chain == pep_chain]

    ck = torch.load(args.model, weights_only=False)
    net = EGNN(node_dim=ck["n_node_dim"], edge_dim=ck["n_edge_dim"],
               depth=6, hidden=128).to(args.device)
    net.load_state_dict(ck["model"])
    net.eval()

    base = score(net, residues, args.device, args.mask_plddt)
    print(f"{args.cif.name}  baseline score={base:.3f}  "
          f"peptide={''.join(AA3TO1.get(r.resname[:3].upper(), 'X') for r in residues if r.chain == pep_chain)}")
    print(f"  {'pos':>3s} {'wt':>4s} {'mut':>4s} {'score':>7s} {'Δscore':>8s}")
    drops = []
    for n, i in enumerate(pep_positions, start=1):
        wt = AA3TO1.get(residues[i].resname[:3].upper(), "X")
        mut = list(residues)
        mut[i] = to_mutant(residues[i], args.mutant)
        s = score(net, mut, args.device, args.mask_plddt)
        drops.append((n, wt, s, base - s))
        print(f"  {n:>3d} {wt:>4s} {args.mutant:>4s} {s:7.3f} {base - s:8.3f}")
    top = sorted(drops, key=lambda x: -x[3])[:3]
    print("largest drops (model depends on these positions): "
          + ", ".join(f"{wt}{n} (−{d:.2f})" for n, wt, _, d in top))
    return 0


if __name__ == "__main__":
    sys.exit(main())