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

# Only truncations can be built correctly here, because the scan keeps the
# backbone coordinates and simply drops atoms. ALA is a real truncation; GLY
# additionally drops CB. Anything else needs a side-chain rotamer, which this
# code does not have — claiming a TRP or ARG scan by renaming the residue
# while keeping N/CA/C/O/CB is exactly what the 2026-10-15 audit found
# (AUDIT.md #15): the three "mutants" had identical atom sets, so the measured
# effect was a relabelling, not a steric or chemical change.
SUPPORTED_MUTANTS = {"ALA": ("N", "CA", "C", "O", "CB"),
                     "GLY": ("N", "CA", "C", "O")}


def to_mutant(res: Residue, target: str = "ALA") -> Residue:
    """Truncate a residue to the target's backbone-only representation.

    Supported: ALA, GLY. Any other target raises: building that side chain
    from scratch (rotamers, ideal geometry, clash relaxation) is out of scope,
    and silently renaming the residue produces a fake scan.
    """
    target = target.upper()
    if target not in SUPPORTED_MUTANTS:
        raise ValueError(
            f"cannot build {target} by truncation; supported targets are "
            f"{sorted(SUPPORTED_MUTANTS)}. A TRP/ARG 'mutation' here would only "
            f"rename the residue while keeping N/CA/C/O/CB, which the audit "
            f"identified as a virtual mutation, not a real one.")
    if target == "ALA" and not any(a[0] == "CB" for a in res.atoms):
        raise ValueError("ALA truncation requires an existing CB; cannot build alanine from glycine or incomplete atoms")
    keep = SUPPORTED_MUTANTS[target]
    atoms = [a for a in res.atoms if a[0] in keep]
    if target == "GLY":
        atoms = [a for a in atoms if a[0] != "CB"]
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
    from meta_acdc.models.train_egnn import normalize_node_features
    h = normalize_node_features(h)
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
                    choices=sorted(SUPPORTED_MUTANTS),
                    help="truncation target. ALA and GLY are real truncations; "
                         "ARG/TRP are refused because this scan cannot build "
                         "their side chains (see AUDIT.md #15)")
    ap.add_argument("--cumulative", action="store_true",
                    help="greedily accumulate the most damaging mutations and "
                         "report score vs number of mutations (distributed vs "
                         "hotspot: gradual drop = distributed)")
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

    if args.cumulative:
        # greedily accumulate the most damaging mutations; a gradual drop
        # means a distributed decision, a cliff means a hotspot
        current = list(residues)
        remaining = list(pep_positions)
        print(f"  {'n_mut':>5s} {'added':>6s} {'score':>7s} {'Δscore':>8s}")
        cur = base
        k = 0
        while remaining:
            best = None
            for i in remaining:
                trial = list(current)
                trial[i] = to_mutant(residues[i], args.mutant)
                s = score(net, trial, args.device, args.mask_plddt)
                if best is None or s < best[1]:
                    best = (i, s)
            i, s = best
            current[i] = to_mutant(residues[i], args.mutant)
            remaining.remove(i)
            k += 1
            wt = AA3TO1.get(residues[i].resname[:3].upper(), "X")
            pos = pep_positions.index(i) + 1
            print(f"  {k:>5d} {wt + str(pos):>6s} {s:7.3f} {cur - s:8.3f}")
            cur = s
        return 0

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