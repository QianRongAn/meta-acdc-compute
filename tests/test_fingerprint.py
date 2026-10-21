"""Regression tests for STRUCT-FINGERPRINT-20261021 machinery."""
import json
import math
import unittest
from pathlib import Path

import numpy as np
import torch

from meta_acdc.models.train_fingerprint import (
    auroc,
    evaluate,
    group_folds,
    set_seed,
    supcon_loss,
)

ROOT = Path(__file__).resolve().parents[1]


class TestGroupFolds(unittest.TestCase):
    def test_deterministic_and_family_disjoint(self):
        fams = [f"F{i}" for i in range(37)]
        a = group_folds(fams, 5)
        b = group_folds(fams, 5)
        self.assertEqual(a, b)
        flat = [f for fold in a for f in fold]
        self.assertEqual(sorted(flat), sorted(fams))
        for f in fams:
            self.assertEqual(sum(f in fold for fold in a), 1)

    def test_same_family_never_split(self):
        fams = ["A"] * 10 + ["B"] * 10 + ["C"] * 10
        folds = group_folds(fams, 3)
        for f in ("A", "B", "C"):
            self.assertEqual(sum(f in fold for fold in folds), 1)


class TestSupCon(unittest.TestCase):
    def test_loss_finite_and_lower_for_clustered(self):
        set_seed()
        embs = torch.randn(8, 16)
        embs = embs / embs.norm(dim=1, keepdim=True)
        rids = ["r1", "r1", "r2", "r2", "r3", "r3", "r4", "r4"]
        loss = supcon_loss(embs, rids)
        self.assertTrue(math.isfinite(float(loss)))
        # perfectly clustered embeddings (twin copies of 4 directions) give
        # a much lower loss; the twin of direction i sits at index i + 4
        base = torch.randn(4, 16)
        base = base / base.norm(dim=1, keepdim=True)
        clustered = torch.cat([base, base])
        clustered_rids = ["r1", "r2", "r3", "r4", "r1", "r2", "r3", "r4"]
        loss2 = supcon_loss(clustered, clustered_rids)
        self.assertLess(float(loss2), float(loss))


class TestAuroc(unittest.TestCase):
    def test_perfect_and_chance(self):
        s = np.array([0.1, 0.2, 0.3, 0.4])
        l = np.array([0, 0, 1, 1])
        self.assertAlmostEqual(auroc(s, l), 1.0)
        # all-tied scores with balanced labels give exactly 0.5
        self.assertAlmostEqual(auroc(np.ones(4), np.array([1, 0, 1, 0])), 0.5)
        # inverted scores give 0.0
        self.assertAlmostEqual(auroc(s, np.array([1, 1, 0, 0])), 0.0)
        self.assertTrue(math.isnan(auroc(s, np.array([0, 0, 0, 0]))))

    def test_inverted(self):
        s = np.array([0.4, 0.3, 0.2, 0.1])
        l = np.array([0, 0, 1, 1])
        self.assertAlmostEqual(auroc(s, l), 0.0)


class TestEvaluate(unittest.TestCase):
    def _records(self, n_rep=4, n_rec=8):
        recs = []
        for r in range(n_rec):
            for i in range(n_rep):
                recs.append({
                    "receptor_id": f"rec{r}",
                    "family_id": f"fam{r % 3}",
                    "domain": "crystal" if i % 2 == 0 else "predicted",
                    "complex_id": f"{r}_{i}",
                })
        return recs

    def test_perfect_embeddings_retrieve(self):
        n_rep, n_rec = 4, 8
        recs = self._records(n_rep, n_rec)
        # one-hot-ish orthogonal embeddings per receptor
        E = torch.zeros(len(recs), n_rec)
        for i, r in enumerate(recs):
            E[i, int(r["receptor_id"][3:])] = 1.0
        E = E / E.norm(dim=1, keepdim=True)
        test = list(range(len(recs)))
        res = evaluate(E, recs, test, None)
        self.assertEqual(res["retrieval"]["rank1"], 1.0)
        self.assertEqual(res["retrieval"]["mrr"], 1.0)
        self.assertEqual(res["retrieval"]["auroc"], 1.0)
        self.assertEqual(res["family"]["1"], 1.0)
        # chance rank-1 equals the share of same-receptor gallery items
        expected_chance = (n_rep - 1) / (len(recs) - 1)
        self.assertAlmostEqual(res["retrieval"]["chance_rank1"],
                               expected_chance, places=6)

    def test_random_embeddings_near_chance(self):
        recs = self._records()
        E = torch.randn(len(recs), 16)
        E = E / E.norm(dim=1, keepdim=True)
        res = evaluate(E, recs, list(range(len(recs))), None)
        # random embeddings: rank-1 close to the per-query chance
        self.assertAlmostEqual(
            res["retrieval"]["rank1"],
            res["retrieval"]["chance_rank1"], delta=0.15)
        self.assertLess(abs(res["retrieval"]["auroc"] - 0.5), 0.15)


class TestPlanArtifacts(unittest.TestCase):
    def test_plan_committed_and_consistent(self):
        plan_dir = ROOT / "docs/fingerprint_retrieval_20261021"
        plan = json.loads((plan_dir / "plan.json").read_text())
        self.assertEqual(plan["task_id"], "STRUCT-FINGERPRINT-20261021")
        self.assertFalse(plan["binding_prediction"])
        self.assertTrue(plan["success_criteria"]["A_rank1_ge_2x_chance"])


if __name__ == "__main__":
    unittest.main()
