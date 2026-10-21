"""Build the STRUCT-FINGERPRINT-20261021 dataset (pre-registered).

Pool: v7 native inventory graphs (168, confidence-masked) + scored candidate
CIFs whose chain A/B sequence pair exactly matches an af3_pending_upload
role-header reference (alpha,beta). Identity-verified receptor labels only;
no binding labels anywhere.

Usage:
    .venv/bin/python scripts/build_fingerprint_dataset.py \
        --inventory data/releases/20261017-structural-v7-atomic-repair/native_inventory_graphs.pt \
        --cifs data/raw/af3_predictions \
        --refs af3_pending_upload \
        --out data/processed/fingerprint_dataset.pt
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import torch

from meta_acdc.data.release import sha256
from meta_acdc.structure.cif import parse_cif
from meta_acdc.structure.graph import build_interface_graph_from_residues
from meta_acdc.structure.identity import chain_sequence
from meta_acdc.structure.quality import annotate_prediction


def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()


def read_role_fastas(refs_dir: Path) -> dict[tuple[str, str], dict]:
    """(alpha,beta) -> reference dict (alpha/beta/mhc_sequence/b2m_sequence/peptide)."""
    out: dict[tuple[str, str], dict] = {}
    for f in sorted(refs_dir.rglob("*.fasta")):
        rec: dict[str, str] = {}
        key = None
        for line in f.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                key = line[1:]
                rec[key] = ""
            elif key:
                rec[key] += line
        if set(rec) != {"alpha", "beta", "peptide", "mhc", "b2m"}:
            continue
        pair = (rec["alpha"], rec["beta"])
        if pair in out and out[pair]["peptide"] != rec["peptide"]:
            continue  # ambiguous reference; keep the first
        out[pair] = {"alpha": rec["alpha"], "beta": rec["beta"],
                     "mhc_sequence": rec["mhc"], "b2m_sequence": rec["b2m"],
                     "peptide": rec["peptide"]}
    return out


def arrays(g):
    return (g.node_features, g.node_coords, g.edge_index, g.edge_features)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inventory", type=Path, required=True)
    ap.add_argument("--cifs", type=Path, default=Path("data/raw/af3_predictions"))
    ap.add_argument("--refs", type=Path, default=Path("af3_pending_upload"))
    ap.add_argument("--out", type=Path,
                    default=Path("data/processed/fingerprint_dataset.pt"))
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError("refusing to overwrite an existing dataset")

    refs = read_role_fastas(args.refs)
    print(f"reference receptors from pending panel: {len(refs)}", flush=True)

    graphs, records, labels = [], [], []
    seen_hashes: set[str] = set()
    n_dup = 0
    graph_hashes: list[str] = []

    inv = torch.load(args.inventory, map_location="cpu", weights_only=False)
    for a in inv["assets"]:
        meta = a["record"]
        g = a["graph"]
        digest = a["input_sha256"]
        if digest in seen_hashes:
            n_dup += 1
            continue
        seen_hashes.add(digest)
        graphs.append(g)
        records.append({
            "complex_id": meta["complex_id"],
            "domain": meta["domain"],
            "pdb": meta["pdb"],
            "peptide": meta["peptide"],
            "alpha": meta["alpha"],
            "beta": meta["beta"],
            "receptor_id": seq_hash(meta["alpha"] + "|" + meta["beta"]),
            "family_id": meta["family_id"],
            "source": "v7-inventory",
            "peptide_mask": a["peptide_mask"],
        })
        labels.append({
            "receptor_id": seq_hash(meta["alpha"] + "|" + meta["beta"]),
            "family_id": meta["family_id"],
            "domain": meta["domain"],
        })
        graph_hashes.append(digest)
    print(f"v7 inventory assets loaded: {len(inv['assets'])} "
          f"(kept {len(graphs)}, dup {n_dup})", flush=True)

    # v7 receptor -> family lookup for candidate labels
    family_of_receptor = {r["receptor_id"]: r["family_id"] for r in records}

    matched = skipped = failed = 0
    for cif in sorted(args.cifs.glob("*_model_0.cif")):
        job = cif.name.rsplit("_model_0", 1)[0]
        pdb, _, peptide = job.partition("_")
        try:
            residues = parse_cif(cif, confidence_kind="af3_plddt")
        except Exception as exc:
            skipped += 1
            continue
        chains: dict[str, list] = defaultdict(list)
        for r in residues:
            chains[r.chain].append(r)
        if "A" not in chains or "B" not in chains:
            skipped += 1
            continue
        a_seq = chain_sequence(chains["A"])
        b_seq = chain_sequence(chains["B"])
        ref = refs.get((a_seq, b_seq)) or refs.get((b_seq, a_seq))
        if ref is None:
            skipped += 1  # identity not verified -> excluded by design
            continue
        try:
            # mask confidence exactly like the v7 inventory path
            from dataclasses import replace
            residues = [replace(r, plddt=None) for r in residues]
            annotate_prediction(residues, ref, peptide)
            g = build_interface_graph_from_residues(residues)
            if g.n_nodes < 20:
                raise ValueError("small interface")
            if len(g.node_features[0]) != 26 or len(g.edge_features[0]) != 16:
                raise ValueError("feature schema mismatch vs v7 inventory")
        except ValueError as exc:
            failed += 1
            continue
        digest = sha256(cif)
        if digest in seen_hashes:
            n_dup += 1
            continue
        seen_hashes.add(digest)
        receptor_id = seq_hash(ref["alpha"] + "|" + ref["beta"])
        family_id = family_of_receptor.get(receptor_id, "pending:" + receptor_id)
        graphs.append(arrays(g))
        records.append({
            "complex_id": job,
            "domain": "predicted",
            "pdb": pdb,
            "peptide": peptide,
            "alpha": ref["alpha"],
            "beta": ref["beta"],
            "receptor_id": receptor_id,
            "family_id": family_id,
            "source": "scored-candidate",
            "peptide_mask": g.peptide_mask,
        })
        labels.append({
            "receptor_id": receptor_id,
            "family_id": family_id,
            "domain": "predicted",
        })
        graph_hashes.append(digest)
        matched += 1
    print(f"candidates: matched {matched}, identity-skipped {skipped}, "
          f"graph-failed {failed}", flush=True)

    torch.save({"graphs": graphs, "records": records, "labels": labels,
                "graph_hashes": graph_hashes,
                "task_id": "STRUCT-FINGERPRINT-20261021",
                "binding_labels": False}, args.out)
    manifest = {
        "task_id": "STRUCT-FINGERPRINT-20261021",
        "n_graphs": len(graphs),
        "n_v7_inventory": sum(1 for r in records if r["source"] == "v7-inventory"),
        "n_candidates": sum(1 for r in records if r["source"] == "scored-candidate"),
        "n_receptors": len({r["receptor_id"] for r in records}),
        "n_families": len({r["family_id"] for r in records}),
        "n_crystal": sum(1 for r in records if r["domain"] == "crystal"),
        "n_predicted": sum(1 for r in records if r["domain"] == "predicted"),
        "receptors_with_ge2_structures": sum(
            1 for v in __import__("collections").Counter(
                r["receptor_id"] for r in records).values() if v >= 2),
        "input_sha256_of_dataset_file": sha256(args.out),
        "binding_labels": False,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
