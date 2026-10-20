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
        mapping = {"A": TCR, "B": TCR, "C": PEPTIDE, "D": MHC, "E": B2M}
        for residue in residues: residue.role = mapping[residue.chain]
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
        from meta_acdc.structure.identity import chain_sequence
        from meta_acdc.structure.quality import annotate_prediction
        from meta_acdc.structure.relabel_instances_by_chains import (
            peptide_of, tcr_signature)
        res = parse_cif(f)
        # This historical instance stores the AF3 chain layout
        # A=tcr alpha, B=tcr beta, C=peptide, D=MHC, E=b2m. Since the
        # review (item #17) chain roles are only assigned from explicit
        # annotations — here via annotate_prediction with the release-style
        # reference, as the current scoring path does.
        seqs = {c: chain_sequence([r for r in res if r.chain == c])
                for c in {r.chain for r in res}}
        self.assertEqual(seqs["C"], "LLFGYPVYV")
        self.assertEqual(len(seqs), 5)
        reference = {"alpha": seqs["A"], "beta": seqs["B"],
                     "mhc_sequence": seqs["D"], "b2m_sequence": seqs["E"]}
        annotate_prediction(res, reference, "LLFGYPVYV")
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

    def test_legacy_fallback_removed_quarantines_by_default(self):
        # The old default silently named an unmatched job after candidates[0]
        # (the 2026-10-14 misattribution). Since the audit, the default
        # QUARANTINES the job: nothing is written, a quarantine list is.
        import csv
        self._run()
        cifs = sorted(p.name for p in self.out.glob("*_model_0.cif"))
        self.assertEqual(len(cifs), 0, f"unmatched jobs must not be named: {cifs}")
        qpath = Path(self.tmp.name) / "af3_import_quarantine.tsv"
        self.assertTrue(qpath.exists(), "quarantine list missing")
        with open(qpath) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(len(rows), 2, "both unmatched jobs quarantined")


class TestComplexityFilter(unittest.TestCase):
    def test_low_complexity_rejected(self):
        from meta_acdc.active_learning.select_candidates import (
            peptide_complexity_ok)
        self.assertFalse(peptide_complexity_ok("QQQEQQHQAW", 0.5))
        self.assertFalse(peptide_complexity_ok("AAAAAAAAA", 0.5))
        self.assertTrue(peptide_complexity_ok("GILGFVFTL", 0.5))
        self.assertTrue(peptide_complexity_ok("LLFGYPVYV", 0.5))
        # filter disabled
        self.assertTrue(peptide_complexity_ok("QQQEQQHQAW", 1.0))


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
        from meta_acdc.structure.ensemble_mean import combine
        rows,_=combine([{'x':.2,'y':.8},{'x':.4,'y':.6},{'x':.6,'y':.4}],['x','y'])
        self.assertAlmostEqual({x['job_id']:x['score'] for x in rows}['x'],.4)

    def test_degenerate_seed_retained_and_not_quality_eligible(self):
        from meta_acdc.structure.ensemble_mean import combine
        rows,diag=combine([{'a':.2,'b':.8},{'a':.3,'b':.7},{'a':0.,'b':0.}],['a','b'])
        self.assertEqual(diag['degenerate_seeds_retained'],[2])
        by={x['job_id']:x for x in rows}
        self.assertAlmostEqual(by['a']['score'],(.2+.3)/3)
        self.assertAlmostEqual(by['b']['score'],(.8+.7)/3)
        self.assertTrue(all(x['quality_eligible']=='no' for x in rows))

    def test_legacy_ensemble_cli_without_release_is_blocked(self):
        import subprocess,sys,tempfile
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'scores.tsv'
            p=subprocess.run([sys.executable,'-m','meta_acdc.structure.ensemble_mean',
                '--scores','old.tsv','--out',str(out)],capture_output=True,cwd=ROOT)
            self.assertNotEqual(p.returncode,0)
            self.assertFalse(out.exists())


class TestAttribution(unittest.TestCase):
    """Gradient attribution runs and reports peptide-involving contacts."""

    def test_to_mutant_truncation_only(self):
        # AUDIT.md #15: only backbone-truncating mutants are real. TRP/ARG
        # "mutants" used to keep N/CA/C/O/CB under a renamed residue — a
        # relabelling, not a mutation. They must now raise.
        from meta_acdc.structure.alanine_scan import to_mutant
        from meta_acdc.structure.graph import Residue
        r = Residue("C", "ILE", 5, 0.0, 0.0, 0.0, None, "",
                    [("N", 0, 0, 0), ("CA", 1, 0, 0), ("CB", 2, 0, 0),
                     ("CG1", 3, 0, 0), ("CD1", 4, 0, 0)])
        ala = to_mutant(r, "ALA")
        self.assertEqual(ala.resname, "ALA")
        self.assertEqual([a[0] for a in ala.atoms],
                         ["N", "CA", "CB"])  # backbone + CB kept
        gly = to_mutant(r, "GLY")
        self.assertEqual(gly.resname, "GLY")
        self.assertEqual([a[0] for a in gly.atoms], ["N", "CA"])  # no CB
        with self.assertRaises(ValueError):
            to_mutant(r, "TRP")
        with self.assertRaises(ValueError):
            to_mutant(r, "ARG")

    def test_attribution_refuses_legacy_checkpoint(self):
        import subprocess
        import sys
        import tempfile
        cif = ROOT / "data/raw/af3_predictions/5brz_ESDPIVAQY_model_0.cif"
        model = ROOT / "data/processed/da_seed0.model.pt"
        if not (cif.exists() and model.exists()):
            self.skipTest("attribution fixtures missing")
        # The legacy checkpoint has no release fingerprint; since the
        # delivery review, attribution must not run on it (and must not
        # guess chain identity). A valid run needs --reference-json +
        # --peptide and a release-bound checkpoint, none of which exist
        # while structural retraining is unreleased.
        with tempfile.TemporaryDirectory() as td:
            ref = Path(td) / "ref.json"
            ref.write_text(json.dumps({"alpha": "A" * 100,
                                       "beta": "B" * 100,
                                       "mhc_sequence": "M" * 270,
                                       "b2m_sequence": "K" * 100}))
            r = subprocess.run(
                [sys.executable, "-m", "meta_acdc.structure.attribute",
                 "--cif", str(cif), "--model", str(model), "--mask-plddt",
                 "--top", "5", "--reference-json", str(ref),
                 "--peptide", "ESDPIVAQY"],
                cwd=str(ROOT), capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("release-bound model", r.stderr)


class TestCandidateRouting(unittest.TestCase):
    """native_manifest routes non-native jobs to the candidate dir."""

    def _run_manifest(self, td: Path, seqs: list[str], pdb: str = "9zzz"):
        import csv
        import json
        import subprocess
        import sys
        lst = td / "list.tsv"
        # the real submission list stores tcr_a = beta, tcr_b = alpha
        with open(lst, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["pdb", "tcr_a", "tcr_b", "peptide", "mhc_hint",
                        "evidence", "pdb_peptide", "cdr3"])
            w.writerow([pdb, seqs[2], seqs[1], seqs[0], "HLA-A",
                        "different", "LLFGYPVYV", "CXXX"])
        job = td / "folds" / "j1"
        job.mkdir(parents=True, exist_ok=True)
        with open(job / "j1_job_request.json", "w") as fh:
            json.dump([{"sequences": [{"proteinChain": {"sequence": s}}
                                      for s in seqs]}], fh)
        (job / "j1_model_0.cif").write_text("data\n")
        subprocess.run([sys.executable, "-m",
                        "meta_acdc.structure.native_manifest",
                        "--src", str(td / "folds"), "--out", str(td / "native"),
                        "--manifest", str(td / "m.tsv"),
                        "--report", str(td / "r.tsv"),
                        "--candidate-out", str(td / "cand"),
                        "--list", str(lst)],
                       check=True, cwd=str(ROOT), capture_output=True)

    def test_candidate_routed_not_dropped(self):
        # a cross-reactivity candidate whose chain pair matches the list is
        # routed to the candidate dir (peptide-role TCR chains present)
        import tempfile
        with tempfile.TemporaryDirectory() as td2:
            td = Path(td2)
            seqs = ["LLFGPVYV",                          # peptide
                    "Q" * 150 + "YFCAV" + "T" * 40,      # alpha (motif)
                    "G" * 190 + "YFCAS" + "V" * 48,      # beta (motif)
                    "M" * 100, "H" * 274]
            self._run_manifest(td, seqs)
            self.assertEqual(len(list((td / "native").glob("*_model_0.cif"))), 0)
            routed = [p.name for p in (td / "cand").glob("*_model_0.cif")]
            self.assertEqual(routed, ["9zzz_LLFGPVYV_model_0.cif"])

    def test_unrecognizable_tcr_quarantined_not_crash(self):
        # a job with no recognizable TCR chain must be quarantined loudly,
        # not crash on a None chain prefix (regression for the 2026-10-15
        # quarantine-path fix)
        import tempfile
        with tempfile.TemporaryDirectory() as td2:
            td = Path(td2)
            seqs = ["LLFGPVYV", "A" * 190, "B" * 200, "M" * 100, "H" * 274]
            self._run_manifest(td, seqs)
            self.assertEqual(len(list((td / "native").glob("*_model_0.cif"))), 0)
            self.assertEqual(len(list((td / "cand").glob("*_model_0.cif"))), 0)


class TestClinicalScan(unittest.TestCase):
    """A high structural score must never fabricate a clinical hazard label."""
    def test_legacy_fatal_mimicry_rule_is_blocked(self):
        from meta_acdc.structure.clinical_scan import main
        with self.assertRaisesRegex(ValueError,'Clinical scoring blocked'):
            main(['--scores','does-not-need-to-exist.tsv','--clinical','legacy.tsv'])


class TestElectrostatics(unittest.TestCase):
    """v11 atom-level electrostatics features (ff14SB united-atom charges)."""

    def test_charge_tables_net_charge(self):
        from meta_acdc.structure.charges import unified_charges
        self.assertAlmostEqual(sum(unified_charges("LYS").values()), 1.0)
        self.assertAlmostEqual(sum(unified_charges("ARG").values()), 1.0)
        self.assertAlmostEqual(sum(unified_charges("ASP").values()), -1.0)
        self.assertAlmostEqual(sum(unified_charges("GLU").values()), -1.0)
        for aa in ("ALA", "HIS", "PRO", "SER", "TRP"):
            self.assertAlmostEqual(sum(unified_charges(aa).values()), 0.0)
        # HIS averages HID/HIE tautomers; ASH is protonated ASP alias
        self.assertEqual(unified_charges("ASH"), unified_charges("ASP"))
        # terminal variants carry the +1/-1 terminus
        self.assertAlmostEqual(
            sum(unified_charges("LYS", first_in_chain=True).values()), 2.0)
        self.assertAlmostEqual(
            sum(unified_charges("ASP", has_oxt=True).values()), -2.0)

    def test_coulomb_symmetry(self):
        import numpy as np
        from meta_acdc.structure.electrostatics import edge_coulomb
        a = np.array([[0., 0., 0., 0.5], [1., 0., 0., -0.5]])
        b = np.array([[0., 2., 0., -0.8]])
        c1, m1 = edge_coulomb(a, b)
        c2, m2 = edge_coulomb(b, a)
        self.assertAlmostEqual(c1, c2)
        self.assertAlmostEqual(m1, m2)
        # empty side -> zero
        c3, m3 = edge_coulomb(np.empty((0, 4)), b)
        self.assertEqual((c3, m3), (0.0, 0.0))

    def test_features_appended_and_respond(self):
        from pathlib import Path
        from meta_acdc.models.dataset import (_peptide_chains, make_decoy,
                                              make_displaced_peptide_decoy)
        from meta_acdc.structure.graph import (
            build_interface_graph_from_residues, parse_pdb)
        import numpy as np
        root = Path(__file__).resolve().parents[1]
        p1 = root / "data/raw/structures/1ao7.pdb"
        p2 = root / "data/raw/structures/5brz.pdb"
        if not (p1.exists() and p2.exists()):
            self.skipTest("structure fixtures missing")
        res = parse_pdb(p1)
        donor = next(iter(_peptide_chains(parse_pdb(p2)).values()))
        g0 = build_interface_graph_from_residues(res, electrostatics=False)
        g1 = build_interface_graph_from_residues(res, electrostatics=True)
        # dims: 26->28 node, 16->18 edge; base indices untouched
        self.assertEqual(len(g0.node_features[0]), 26)
        self.assertEqual(len(g0.edge_features[0]), 16)
        self.assertEqual(len(g1.node_features[0]), 28)
        self.assertEqual(len(g1.edge_features[0]), 18)
        for f0, f1 in zip(g0.node_features, g1.node_features):
            self.assertEqual(f0, f1[:26])
        for f0, f1 in zip(g0.edge_features, g1.edge_features):
            self.assertEqual(f0, f1[:16])
        # graft (chemistry swap) and displacement must move the new features
        keys = set(g1.node_keys)
        graft = build_interface_graph_from_residues(
            make_decoy(res, donor), electrostatics=True, forced_nodes=keys)
        disp = build_interface_graph_from_residues(
            make_displaced_peptide_decoy(res), electrostatics=True,
            forced_nodes=keys)
        nat_n = np.array(g1.node_features)[:, 26:]
        gr_n = np.array(graft.node_features)[:, 26:]
        di_n = np.array(disp.node_features)[:, 26:]
        self.assertGreater(np.abs(nat_n - gr_n).sum(), 1e-3)
        self.assertGreater(np.abs(nat_n - di_n).sum(), 1e-3)


class TestPackaging(unittest.TestCase):
    """Deliverable 3: console entry points resolve to importable callables."""

    def test_entry_points_importable(self):
        import importlib
        try:
            import tomllib
        except ModuleNotFoundError:  # py<3.11
            self.skipTest("tomllib unavailable")
        root = Path(__file__).resolve().parents[1]
        with open(root / "pyproject.toml", "rb") as fh:
            cfg = tomllib.load(fh)
        scripts = cfg.get("project", {}).get("scripts", {})
        self.assertGreaterEqual(len(scripts), 8)
        for name, target in scripts.items():
            mod, _, func = target.partition(":")
            m = importlib.import_module(mod)
            self.assertTrue(callable(getattr(m, func)), f"{name}: {target}")


class TestDashboard(unittest.TestCase):
    def test_build_data_serializable(self):
        from meta_acdc.dashboard.server import (DEFAULT_CLINICAL, DEFAULT_MAP,
                                                DEFAULT_QC, DEFAULT_SCORES,
                                                build_data)
        d = build_data(DEFAULT_SCORES, DEFAULT_MAP, DEFAULT_CLINICAL, DEFAULT_QC)
        json.dumps(d)  # must be JSON-serializable
        self.assertIn("groups", d)
        self.assertIn("caveats", d)
        self.assertEqual(d["groups"], {})
        self.assertEqual(d["clinical_scan"], [])
        self.assertEqual(d["clinical"], [])


if __name__ == "__main__":
    unittest.main()
