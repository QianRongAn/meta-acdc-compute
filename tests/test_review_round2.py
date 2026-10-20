"""Regression tests for the delivery-review round-2 remediation (2026-10-21).

Covers the remaining review items implemented in this round:
- #32 family-stratified Fisher pooling (rank_stats)
- #14 direction-aware displaced decoy (dataset)
- #34 manifest duplicate/legacy guards (import_af3.load_manifest)
- #16/#33 frozen panel registry reproduction
- #2  resolution-floor report machinery
- #20 random-peptide baseline machinery (refusal + construction helpers)
"""
import contextlib
import csv
import io
import json
import random
import tempfile
import unittest
from pathlib import Path

from scipy.stats import chi2

from meta_acdc.models.dataset import make_displaced_peptide_decoy
from meta_acdc.structure.graph import MHC, PEPTIDE, Residue
from meta_acdc.structure.import_af3 import load_manifest
from meta_acdc.structure.rank_stats import (
    family_stratified_fisher,
    fisher_pool,
)

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "data/processed/candidate_panel_registry.tsv"
REGISTRY_SUMMARY = ROOT / "data/processed/candidate_panel_registry.summary.json"
RECOVERED = ROOT / "docs/delivery_audit_20261018/recovered"


def _res(chain, resname, resid, x, y, z, role):
    return Residue(chain=chain, resname=resname, resid=resid, x=x, y=y, z=z,
                   atoms=[("N", x, y, z), ("CA", x, y, z),
                          ("C", x, y, z), ("O", x, y, z)], role=role)


class TestFisherPooling(unittest.TestCase):
    def test_exact_two_rank1_of_10(self):
        # chi2 = -2 * 2 * log(0.1) = 9.2103..., df = 4; exact P ~ 0.05605
        chi2_stat, p = fisher_pool([0.1, 0.1])
        self.assertAlmostEqual(chi2_stat, -4 * __import__("math").log(0.1),
                               places=6)
        self.assertAlmostEqual(p, chi2.sf(chi2_stat, 4), places=10)
        self.assertAlmostEqual(p, 0.05605, places=4)

    def test_single_pvalue_passthrough(self):
        chi2_stat, p = fisher_pool([0.25])
        self.assertAlmostEqual(p, chi2.sf(chi2_stat, 2), places=10)

    def test_rejects_invalid(self):
        with self.assertRaises(ValueError):
            fisher_pool([])
        with self.assertRaises(ValueError):
            fisher_pool([0.0])
        with self.assertRaises(ValueError):
            fisher_pool([1.5])


class TestFamilyStratifiedFisher(unittest.TestCase):
    def _rows(self):
        # two families, two structures each
        return [
            ("p1", "NAT1", 1, 10),
            ("p2", "NAT2", 1, 10),   # family A: both rank 1/10
            ("p3", "NAT3", 10, 10),
            ("p4", "NAT4", 10, 10),  # family B: both rank 10/10
        ]

    def test_family_level_pooling_matches_scipy(self):
        fam = {"p1": "A", "p2": "A", "p3": "B", "p4": "B"}
        chi2_stat, p, fam_ps = family_stratified_fisher(self._rows(), fam)
        # within-family A: fisher_pool([0.1, 0.1]); B: fisher_pool([1.0, 1.0])
        self.assertAlmostEqual(fam_ps["A"], chi2.sf(-4 * __import__("math").log(0.1), 4), places=10)
        self.assertAlmostEqual(fam_ps["B"], 1.0, places=10)
        # across families: fisher_pool([fam_ps["A"], fam_ps["B"]])
        expect_chi2 = -2 * (__import__("math").log(fam_ps["A"]) + __import__("math").log(fam_ps["B"]))
        self.assertAlmostEqual(chi2_stat, expect_chi2, places=6)
        self.assertAlmostEqual(p, chi2.sf(expect_chi2, 4), places=10)

    def test_single_member_family_passes_through(self):
        import math
        fam = {"p1": "A", "p2": "B", "p3": "B", "p4": "B"}
        _, _, fam_ps = family_stratified_fisher(self._rows(), fam)
        # A has one structure: p-value passes through unchanged
        self.assertAlmostEqual(fam_ps["A"], 0.1, places=10)
        # B holds p2 (rank 1/10) + p3, p4 (rank 10/10): within-family Fisher
        self.assertAlmostEqual(
            fam_ps["B"], chi2.sf(-2 * math.log(0.1), 6), places=10)

    def test_missing_family_key_raises(self):
        fam = {"p1": "A", "p2": "A", "p3": "B"}
        with self.assertRaises(KeyError):
            family_stratified_fisher(self._rows(), fam)


class TestDisplacedDecoyDirection(unittest.TestCase):
    def test_axial_case_moves_away_from_groove(self):
        mhc = [_res("M", "ALA", i + 1, x, 0.0, 0.0, MHC)
               for i, x in enumerate([-1.0, 0.0, 1.0])]
        pep = [_res("P", "GLY", i + 1, 0.0, 0.0, z, PEPTIDE)
               for i, z in enumerate([10.0, 12.0])]
        out = make_displaced_peptide_decoy(mhc + pep, shift=6.0)
        roles = {r.chain: r.role for r in out}
        moved = [r for r in out if roles[r.chain] == PEPTIDE]
        still = [r for r in out if roles[r.chain] == MHC]
        # direction = peptide centroid (0,0,11) - groove centroid (0,0,0) -> +z
        for r, z0 in zip(moved, [10.0, 12.0]):
            self.assertAlmostEqual(r.z, z0 + 6.0, places=9)
            self.assertTrue(all(abs(a[3] - (z0 + 6.0)) < 1e-9 for a in r.atoms))
        for r in still:
            self.assertAlmostEqual(r.z, 0.0, places=9)
        # centroid distance grows by exactly the shift
        def cz(rs, sel):
            pts = [a[3] for r in rs if r.chain == sel for a in r.atoms]
            return sum(pts) / len(pts)
        self.assertAlmostEqual(cz(out, "P") - cz(out, "M"), 11.0 + 6.0, places=9)

    def test_non_axial_direction_is_data_dependent(self):
        mhc = [_res("M", "ALA", 1, 0.0, 0.0, 0.0, MHC)]
        pep = [_res("P", "GLY", 2, 10.0, 0.0, 0.0, PEPTIDE)]
        out = make_displaced_peptide_decoy(mhc + pep, shift=5.0)
        moved = [r for r in out if r.role == PEPTIDE][0]
        # global +z shift would leave x untouched; direction-aware moves +x
        self.assertAlmostEqual(moved.x, 15.0, places=9)
        self.assertAlmostEqual(moved.z, 0.0, places=9)

    def test_missing_mhc_raises(self):
        pep = [_res("P", "GLY", 1, 0.0, 0.0, 10.0, PEPTIDE)]
        with self.assertRaisesRegex(ValueError, "MHC-role"):
            make_displaced_peptide_decoy(pep)

    def test_degenerate_centroids_raise(self):
        mhc = [_res("M", "ALA", 1, 0.0, 0.0, 0.0, MHC)]
        pep = [_res("P", "GLY", 2, 0.0, 0.0, 0.0, PEPTIDE)]
        with self.assertRaisesRegex(ValueError, "degenerate"):
            make_displaced_peptide_decoy(mhc + pep)

    def test_invalid_shift_raises(self):
        mhc = [_res("M", "ALA", 1, 0.0, 0.0, 0.0, MHC)]
        pep = [_res("P", "GLY", 2, 0.0, 0.0, 10.0, PEPTIDE)]
        with self.assertRaises(ValueError):
            make_displaced_peptide_decoy(mhc + pep, shift=0.0)
        with self.assertRaises(ValueError):
            make_displaced_peptide_decoy(mhc + pep, shift=-3.0)


class TestManifestGuards(unittest.TestCase):
    def _write(self, d: Path, name: str, text: str) -> Path:
        p = d / name
        p.write_text(text)
        return p

    def test_duplicate_rows_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p = self._write(Path(td), "list.tsv",
                            "pdb\ttcr_a\ttcr_b\tpeptide\n"
                            "5brz\tAAA\tBBB\tPEP\n"
                            "5brz\tAAA\tBBB\tPEP\n")
            with self.assertRaisesRegex(ValueError, "duplicate submission row"):
                load_manifest(p)
            rows = load_manifest(p, allow_duplicates=True)
            self.assertEqual(len(rows), 2)

    def test_legacy_swapped_list_announced(self):
        with tempfile.TemporaryDirectory() as td:
            p = self._write(Path(td), "kn5_submission_list.tsv",
                            "pdb\ttcr_a\ttcr_b\tpeptide\n"
                            "5brz\tAAA\tBBB\tPEP\n")
            buf = io.StringIO()
            with contextlib.redirect_stderr(buf):
                rows = load_manifest(p)
            self.assertEqual(len(rows), 1)
            self.assertIn("legacy submission list", buf.getvalue())
            self.assertIn("tcr_a=beta", buf.getvalue())

    def test_non_legacy_file_is_silent(self):
        with tempfile.TemporaryDirectory() as td:
            p = self._write(Path(td), "clean_manifest.tsv",
                            "pdb\ttcr_a\ttcr_b\tpeptide\n"
                            "5brz\tAAA\tBBB\tPEP\n")
            buf = io.StringIO()
            with contextlib.redirect_stderr(buf):
                load_manifest(p)
            self.assertEqual(buf.getvalue(), "")


class TestPanelRegistry(unittest.TestCase):
    def test_committed_registry_reproduces_review_counts(self):
        self.assertTrue(REGISTRY.exists(), "run scripts/build_panel_registry.py")
        by_panel = {"scored": 0, "pending": 0}
        pdbs = {"scored": set(), "pending": set()}
        peptides = {"scored": set(), "pending": set()}
        receptors = {"scored": set(), "pending": set()}
        with open(REGISTRY, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                by_panel[r["panel"]] += 1
                pdbs[r["panel"]].add(r["pdb"])
                peptides[r["panel"]].add(r["peptide"])
                if r["receptor_key"]:
                    receptors[r["panel"]].add(r["receptor_key"])
        # 188 scored CIFs / 72 pdb prefixes / 59 exact A+B pairs,
        # 347 pending FASTAs / 111 pdb prefixes / 91 unique peptides —
        # exactly the independent review counts.
        self.assertEqual(by_panel["scored"], 188)
        self.assertEqual(len(pdbs["scored"]), 72)
        self.assertEqual(len(receptors["scored"]), 59)
        self.assertEqual(by_panel["pending"], 347)
        self.assertEqual(len(pdbs["pending"]), 111)
        self.assertEqual(len(peptides["pending"]), 91)

    def test_summary_matches_rows(self):
        summary = json.loads(REGISTRY_SUMMARY.read_text())
        with open(REGISTRY, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        for panel in summary["panels"]:
            sub = [r for r in rows if r["panel"] == panel["panel"]]
            self.assertEqual(panel["files"], len(sub))
            self.assertEqual(
                panel["pdb_prefixes"], len({r["pdb"] for r in sub}))
            self.assertEqual(
                panel["unique_peptides"], len({r["peptide"] for r in sub}))

    def test_scored_rows_carry_identity_scope_note(self):
        with open(REGISTRY, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["panel"] == "scored":
                    self.assertIn("release metadata", r["receptor_identity_note"])
                    break


class TestResolutionFloor(unittest.TestCase):
    def test_report_written_from_recovered_tables(self):
        from meta_acdc.analysis.resolution_floor import main
        seed_tables = [RECOVERED / f"da_scores_seed{i}_ensemble.tsv"
                       for i in range(5)]
        model_table = (ROOT / "docs/delivery_audit_20261018/historical"
                       / "data/processed"
                       / "prediction_scores_ensemble.relabeled.tsv")
        self.assertTrue(all(p.exists() for p in seed_tables))
        self.assertTrue(model_table.exists())
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "resolution_floor.json"
            rc = main(["--seed-scores"] + [str(p) for p in seed_tables]
                      + ["--model-scores", str(model_table),
                         "--reference", "5brz_EVDPIGHLY",
                         "--out", str(out)])
            self.assertEqual(rc, 0)
            report = json.loads(out.read_text())
            self.assertEqual(report["reference"], "5brz_EVDPIGHLY")
            self.assertIn("DESCRIPTIVE ONLY", report["notice"])
            self.assertTrue(any(r["job_id"] == "5brz_ESDPIVAQY"
                                for r in report["rows"]))
            self.assertIsNotNone(report["median_abs_delta_vs_reference"])

    def test_refuses_overwrite_and_missing_reference(self):
        from meta_acdc.analysis.resolution_floor import main
        seed_tables = [RECOVERED / f"da_scores_seed{i}_ensemble.tsv"
                       for i in range(5)]
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "r.json"
            out.write_text("{}\n")
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                main(["--seed-scores"] + [str(p) for p in seed_tables]
                     + ["--reference", "5brz_EVDPIGHLY",
                        "--out", str(out)])
            with self.assertRaisesRegex(ValueError, "not present"):
                main(["--seed-scores"] + [str(p) for p in seed_tables]
                     + ["--reference", "no_such_job",
                        "--out", str(Path(td) / "r2.json")])


class TestRandomPeptideBaseline(unittest.TestCase):
    def test_random_peptide_shells(self):
        from meta_acdc.analysis.random_peptide_baseline import random_peptide
        rng = random.Random(0)
        pep = random_peptide(9, rng)
        self.assertEqual(len(pep), 9)
        for r in pep:
            self.assertEqual(len(r.resname), 1)
            self.assertEqual(r.chain, "P")
            self.assertEqual(r.role, PEPTIDE)

    def test_refuses_without_matching_release_bound_checkpoint(self):
        from meta_acdc.analysis.random_peptide_baseline import main
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            (td / "release_manifest.json").write_text("{}")
            with self.assertRaises(Exception):
                main(["--receptor-cif", str(td / "x.cif"),
                      "--model", str(td / "m.pt"),
                      "--release", str(td / "release_manifest.json"),
                      "--out", str(td / "out.tsv")])


if __name__ == "__main__":
    unittest.main()
