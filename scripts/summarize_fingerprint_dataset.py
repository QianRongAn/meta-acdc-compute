"""Dataset composition summary for the fingerprint benchmark (paper table)."""
import json
import sys
from collections import Counter
from pathlib import Path

import torch


def main() -> int:
    p = Path("data/processed/fingerprint_dataset.pt")
    d = torch.load(p, map_location="cpu", weights_only=False)
    records = d["records"]
    recs = Counter(r["receptor_id"] for r in records)
    fams = Counter(r["family_id"] for r in records)
    rows = {
        "n_graphs": len(records),
        "n_crystal": sum(r["domain"] == "crystal" for r in records),
        "n_predicted": sum(r["domain"] == "predicted" for r in records),
        "n_receptors": len(recs),
        "n_families": len(fams),
        "receptors_with_ge2": sum(v >= 2 for v in recs.values()),
        "receptors_with_ge5": sum(v >= 5 for v in recs.values()),
        "max_structures_per_receptor": max(recs.values()),
        "structures_from_v7_inventory": sum(r["source"] == "v7-inventory"
                                            for r in records),
        "structures_from_verified_candidates": sum(
            r["source"] == "scored-candidate" for r in records),
        "peptide_lengths": dict(Counter(len(r["peptide"]) for r in records)),
        "binding_labels": False,
    }
    out = Path("docs/fingerprint_retrieval_20261021/dataset_composition.json")
    out.write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps(rows, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
