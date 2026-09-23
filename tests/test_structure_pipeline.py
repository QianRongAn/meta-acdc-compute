"""Regression tests for the structure / AF3 pipeline (stdlib unittest).

Run:  .venv/bin/python -m unittest discover -s tests -v

Data-dependent tests skip automatically when the (gitignored) fixtures are
absent, so the suite still runs on a fresh checkout and inside the CPU image.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from meta_acdc.structure.graph import (B2M, MHC, PEPTIDE, TCR, Residue,
                                       classify_chains)
from meta_acdc.structure.import_af3 import matches_chain_pair, tcr_pair

ROOT = Path(__file__).resolve().parents[1]

AA3 = {"A": "ALA", "C": "CYS", "F": "PHE", "Y": "TYR", "V": "VAL",
       "S": "SER", "L": "LEU", "G": "GLY", "M": "MET", "Q": "GLN",
       "T": "THR", "K": "LYS", "E": "GLU", "R": "ARG", "P": "PRO",
       "D": "ASP", "N": "ASN", "I": "ILE", "H": "HIS", "W": "TRP"}


def chain_residues(chain: str, seq: str) -> list[Residue]:
    return [Residue(chain, AA3.get(a, "ALA"), i + 1, 0.0, 0.0, 0.0,
                    None, "", []) for i, a in enumerate(seq)]


class TestChainClassification(unittest.TestCase):
    def test_roles(self):
        residues = []
        residues += chain_residues("A", "Y" * 190 + "YLCAV")       # TCR alpha
        residues += chain_residues("B", "Y" * 238 + "YFCAS")       # TCR beta
        residues += chain_residues("C", "LGYGFVNYI")               # peptide
        residues += chain_residues("D", "G" * 274)                 # MHC-I
        residues += chain_residues("E", "M" * 100)                 # B2M
        roles = classify_chains(residues)
        self.assertEqual(roles["A"], TCR)
        self.assertEqual(roles["B"], TCR)
        self.assertEqual(roles["C"], PEPTIDE)
        self.assertEqual(roles["D"], MHC)
        self.assertEqual(roles["E"], B2M)


class TestChainPairMatching(unittest.TestCase):
    """Guards the tcr_a/tcr_b column-swap bug (2026-09-23)."""

    def setUp(self):
        self.peptide = "LGYGFVNYI"
        self.alpha = "K" * 150 + "YLCAV" + "T" * 45
        self.beta = "G" * 190 + "YFCAS" + "V" * 48
        self.seqs = sorted([self.peptide, self.alpha, self.beta,
                            "M" * 100, "G" * 274], key=len)

    def test_tcr_pair(self):
        a, b = tcr_pair(self.seqs, self.peptide)
        self.assertEqual(a, self.alpha)
        self.assertEqual(b, self.beta)

    def test_matches_swapped_columns(self):
        # kn5_submission_list stores tcr_a = beta, tcr_b = alpha
        cand = {"tcr_a": self.beta, "tcr_b": self.alpha}
        self.assertTrue(matches_chain_pair(self.alpha, self.beta, cand))
        # and an unrelated candidate must not match
        other = {"tcr_a": "X" * 100, "tcr_b": "Z" * 100}
        self.assertFalse(matches_chain_pair(self.alpha, self.beta, other))

    def test_matches_no_match_when_missing(self):
        self.assertFalse(matches_chain_pair(None, self.beta,
                                            {"tcr_a": self.beta,
                                             "tcr_b": self.alpha}))


class TestDataFixtures(unittest.TestCase):
    """Skip unless the local AF3 data tree is present."""

    def _require(self, p: Path):
        if not p.exists():
            self.skipTest(f"fixture missing: {p}")

    def test_cif_parse_and_identity(self):
        f = (ROOT / "data/raw/af3_instances/A"
             / "1ao7_LLFGYPVYV_model_0.cif")
        self._require(f)
        from meta_acdc.structure.cif import parse_cif
        from meta_acdc.structure.graph import PEPTIDE
        from meta_acdc.structure.relabel_instances_by_chains import (
            peptide_of, tcr_signature)
        res = parse_cif(f)
        self.assertEqual(peptide_of(res), "LLFGYPVYV")
        sig = tcr_signature(res)
        self.assertEqual(len(sig), 2)

    def test_relabel_map_present_and_sane(self):
        f = ROOT / "data/processed/af3_relabel_map.tsv"
        self._require(f)
        import csv
        with open(f) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertGreater(len(rows), 0)
        # the A6/B7 swap fix must be recorded
        self.assertTrue(any(r["old_job"].startswith("1ao7_")
                            and r["new_job"].startswith("1qrn_")
                            for r in rows))

    def test_qc_table(self):
        f = ROOT / "data/processed/af3_qc.tsv"
        self._require(f)
        import csv
        with open(f) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertIn("iptm_mean", rows[0])
        self.assertTrue(all("truncated_tcr" not in r["flag"] for r in rows),
                        "truncated_tcr flag was withdrawn (1ao7 is a real "
                        "truncated construct, not an AF3 defect)")


class TestImportUniqueFallback(unittest.TestCase):
    """Guards the job-collapse / mislabel bug in import_af3 (2026-09-23).

    Distinct native jobs sharing a peptide but with different TCRs must not
    collapse onto one pdb via the legacy candidates[0] fallback.
    """

    def setUp(self):
        import csv
        import json
        import tempfile
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.src = root / "folds"
        self.out = root / "out"
        self.list = root / "list.tsv"
        # two TCRs, same peptide; NEITHER chain pair matches the list
        self.pep = "SIYRYYGL"
        self.alpha1 = "Q" * 150 + "YFCAV" + "T" * 40
        self.alpha2 = "R" * 150 + "YFCAV" + "S" * 40
        self.beta = "G" * 190 + "YFCAS" + "V" * 48
        with open(self.list, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["pdb", "tcr_a", "tcr_b", "peptide", "mhc_hint",
                        "evidence", "pdb_peptide", "cdr3"])
            w.writerow(["9zzz", "OTHER" * 20, "ALSO" * 20, self.pep, "HLA-A",
                        "same", self.pep, "CXXX"])
        for tag, alpha in (("j1", self.alpha1), ("j2", self.alpha2)):
            d = self.src / tag
            d.mkdir(parents=True)
            with open(d / f"{tag}_job_request.json", "w") as fh:
                json.dump([{"sequences": [{"proteinChain": {"sequence": s}}
                                          for s in (self.pep, alpha,
                                                    self.beta,
                                                    "M" * 100, "H" * 274)]}],
                          fh)
            (d / f"{tag}_model_0.cif").write_text("data\n")

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, *extra):
        import subprocess
        import sys
        cmd = [sys.executable, "-m", "meta_acdc.structure.import_af3",
               "--src", str(self.src), "--list", str(self.list),
               "--out", str(self.out), *extra]
        subprocess.run(cmd, check=True, capture_output=True,
                       cwd=str(ROOT))

    def test_unique_fallback_no_collapse(self):
        self._run("--unique-fallback", "--report",
                  str(Path(self.tmp.name) / "rep.tsv"))
        cifs = sorted(p.name for p in self.out.glob("*_model_0.cif"))
        self.assertEqual(len(cifs), 2, f"jobs collapsed: {cifs}")
        self.assertTrue(all(c.startswith("nativetcr-") for c in cifs))

    def test_legacy_fallback_still_collapses_for_candidate_route(self):
        # default (no flag) keeps legacy behaviour (candidates[0]) — only
        # used for the cross-reactivity candidate route
        self._run()
        cifs = sorted(p.name for p in self.out.glob("*_model_0.cif"))
        self.assertEqual(len(cifs), 1)


class TestCalibration(unittest.TestCase):
    def test_ece_bounds_and_improvement(self):
        import numpy as np
        from meta_acdc.models.calibrate import ece
        # bin-calibrated: confidence equals the positive rate in the bin
        p = np.full(10, 0.5)
        y = np.array([1, 0, 1, 0, 1, 0, 1, 0, 1, 0])
        self.assertLess(ece(p, y), 1e-6)
        # systematically overconfident -> ECE > 0
        self.assertGreater(ece(np.full(100, 0.99), np.zeros(100)), 0.5)


class TestEnsembleMean(unittest.TestCase):
    def test_average_across_seeds(self):
        import csv
        import subprocess
        import sys
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            seeds = []
            for s, val in enumerate([0.2, 0.4, 0.6]):
                p = td / f"seed{s}.tsv"
                with open(p, "w", newline="") as fh:
                    w = csv.writer(fh, delimiter="\t")
                    w.writerow(["job_id", "score", "std"])
                    w.writerow(["x_PEP", str(val), "0.0"])
                seeds.append(str(p))
            out = td / "mean.tsv"
            subprocess.run([sys.executable, "-m",
                            "meta_acdc.structure.ensemble_mean",
                            "--scores", *seeds, "--out", str(out)],
                           check=True, cwd=str(ROOT), capture_output=True)
            row = next(csv.DictReader(open(out), delimiter="\t"))
            self.assertAlmostEqual(float(row["score"]), 0.4, places=3)

    def test_degenerate_seed_excluded(self):
        import csv
        import subprocess
        import sys
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            paths = []
            for name, rows in (
                ("good1", [("a", 0.2), ("b", 0.8)]),
                ("good2", [("a", 0.3), ("b", 0.7)]),
                ("bad", [("a", 0.0), ("b", 0.0)]),   # constant -> degenerate
            ):
                p = td / f"{name}.tsv"
                with open(p, "w", newline="") as fh:
                    w = csv.writer(fh, delimiter="\t")
                    w.writerow(["job_id", "score", "std"])
                    for j, v in rows:
                        w.writerow([j, str(v), "0.0"])
                paths.append(str(p))
            out = td / "mean.tsv"
            r = subprocess.run([sys.executable, "-m",
                                "meta_acdc.structure.ensemble_mean",
                                "--scores", *paths, "--out", str(out)],
                               check=True, cwd=str(ROOT),
                               capture_output=True, text=True)
            self.assertIn("degenerate", r.stdout)
            rows = {x["job_id"]: float(x["score"])
                    for x in csv.DictReader(open(out), delimiter="\t")}
            self.assertAlmostEqual(rows["a"], 0.25, places=3)
            self.assertAlmostEqual(rows["b"], 0.75, places=3)


class TestClinicalScan(unittest.TestCase):
    """KN-11: a fatal off-target scoring at/above the target must be flagged."""

    def test_fatal_mimicry_flagged(self):
        import csv
        import subprocess
        import sys
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            scores = td / "scores.tsv"
            with open(scores, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["job_id", "score"])
                w.writerow(["5brz_TARGETPEP", "0.9000"])
                w.writerow(["5brz_FATALPEPT", "0.9500"])
                w.writerow(["5brz_SAFEPEPXX", "0.1000"])
            clinical = td / "clinical.tsv"
            with open(clinical, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["tcr", "target", "off_targets", "fatal",
                            "evidence", "note"])
                w.writerow(["X", "TARGETPEP", "FATALPEPT;SAFEPEPXX", "yes",
                            "PMID", ""])
            out = td / "scan.tsv"
            subprocess.run([sys.executable, "-m",
                            "meta_acdc.structure.clinical_scan",
                            "--scores", str(scores), "--clinical", str(clinical),
                            "--out", str(out)], check=True, cwd=str(ROOT),
                           capture_output=True)
            rows = list(csv.DictReader(open(out), delimiter="\t"))
            by = {r["peptide"]: r for r in rows}
            self.assertEqual(by["FATALPEPT"]["flag"], "FATAL-mimicry flagged")
            self.assertEqual(by["TARGETPEP"]["flag"], "cognate")
            self.assertEqual(by["SAFEPEPXX"]["flag"], "")


class TestDashboard(unittest.TestCase):
    def test_build_data_serializable(self):
        from meta_acdc.dashboard.server import (DEFAULT_CLINICAL, DEFAULT_MAP,
                                                DEFAULT_QC, DEFAULT_SCORES,
                                                build_data)
        d = build_data(DEFAULT_SCORES, DEFAULT_MAP, DEFAULT_CLINICAL, DEFAULT_QC)
        json.dumps(d)  # must be JSON-serializable
        self.assertIn("groups", d)
        self.assertIn("caveats", d)


if __name__ == "__main__":
    unittest.main()
