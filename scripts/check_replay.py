"""Same-seed determinism check (success criterion 5 of the follow-up plan).

Compares the key evaluation fields of two runs of train_fingerprint_v2 with
the same seed/fold (main run vs replay) and reports whether they match exactly.
Exit code 0 when all compared fields are bit-identical.

Usage:
    .venv/bin/python scripts/check_replay.py <metrics_main.json> <metrics_replay.json>
"""
import json
import sys
from pathlib import Path

FIELDS = [
    ("retrieval", "rank1"),
    ("retrieval", "mrr"),
    ("retrieval", "auroc"),
    ("retrieval", "chance_rank1"),
    ("family", "1"),
    ("domain", "accuracy"),
]


def main(argv: list[str] | None = None) -> int:
    args = argv or sys.argv[1:]
    if len(args) != 2:
        print(__doc__)
        return 2
    a = json.loads(Path(args[0]).read_text())
    b = json.loads(Path(args[1]).read_text())
    rows = []
    ok = True
    for section, key in FIELDS:
        va = a.get(section, {}).get(key)
        vb = b.get(section, {}).get(key)
        same = va is not None and va == vb
        ok = ok and same
        rows.append((f"{section}.{key}", va, vb, same))
    for name, va, vb, same in rows:
        print(f"{name}: main={va} replay={vb} identical={same}")
    print("DETERMINISM:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
