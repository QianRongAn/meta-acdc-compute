"""Verify frozen inputs and replay fragment discovery and coordinate auditing."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from scripts.join_terminal_rescue import BASE, ROOT


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    plan = json.loads((BASE / "search_plan.json").read_text())
    assert sha(ROOT / plan["source_file"]) == plan["source_sha256"]
    for target in plan["targets"]:
        for source in target["evidence"]:
            assert sha(ROOT / source["annotation_file"]) == source["annotation_sha256"]
    old_root = ROOT / "docs/fingerprint_np105_recovery_20261029"
    old_index = json.loads((old_root / "frozen_integrity_after_storage_recovery.json").read_text())
    folders = [row["folder"] for row in old_index["previous_deliveries"]]
    folders += [old_root.name, "fingerprint_np105_context_20261030",
                "fingerprint_control_closure_20261101", "fingerprint_tk3_j_evidence_20261101"]
    old_files = 0
    for folder in folders:
        base = ROOT / "docs" / folder
        index = json.loads((base / "DELIVERY_MANIFEST.json").read_text())
        for name, receipt in index["files"].items():
            path = base / name
            assert sha(path) == receipt["sha256"], path
            old_files += 1
    assert len(folders) == 18 and old_files == 4795
    tests = subprocess.run([sys.executable, "-m", "unittest", "tests.test_paired_variable_identity",
                            "tests.test_cross_source_retrieval", "tests.test_j_reference_evidence",
                            "tests.test_control_closure"], cwd=ROOT, capture_output=True, text=True)
    assert tests.returncode == 0, tests.stderr
    (out / "tests.log").write_text(tests.stderr)
    for script, output in (("join_terminal_rescue", "search_results.json"),
                           ("audit_terminal_rescue", "audit.json")):
        subprocess.run([sys.executable, "-m", f"scripts.{script}", "--out", str(out / output)],
                       cwd=ROOT, check=True)
        assert (out / output).read_bytes() == (BASE / output).read_bytes(), output
    paths = [
        "src/meta_acdc/structure/graph.py",
        "src/meta_acdc/structure/identity.py",
        "src/meta_acdc/structure/rigid_controls.py",
        "src/meta_acdc/structure/paired_variable_identity.py",
        "docs/structure_repair_20261017/policy_v2.json",
    ]
    for path in paths:
        frozen = subprocess.check_output(["git", "show", f"6e59b478:{path}"], cwd=ROOT)
        assert frozen == (ROOT / path).read_bytes(), path
    manifest = json.loads((BASE / "DELIVERY_MANIFEST.json").read_text())
    for path, expected in manifest["files"].items():
        assert sha(ROOT / path) == expected, path
    result = dict(status="pass", test_summary=tests.stderr.splitlines()[-3:],
                  frozen_prior_deliveries=len(folders), frozen_prior_files=old_files,
                  receipt_and_audit_replays_exact=2, unchanged_policy_or_implementation_files=len(paths),
                  verified_delivery_files=len(manifest["files"]), independent_acceptance_accessed=False)
    (out / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
