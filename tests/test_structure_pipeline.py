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
