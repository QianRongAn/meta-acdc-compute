"""Tests for FOLLOWUP-AUROC-RASHOMON-20261021 machinery
(train_fingerprint_v2 + rashomon_fingerprint)."""
import json
import math
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import torch

import scripts.rashomon_fingerprint as rashomon
from meta_acdc.models.train_fingerprint import set_seed
from meta_acdc.models.train_fingerprint_v2 import (
    EmbedderV2,
    family_neg_mask,
    supcon_loss_v2,
)

ROOT = Path(__file__).resolve().parents[1]


def _unit(embs):
    return embs / embs.norm(dim=1, keepdim=True).clamp(min=1e-8)


class TestFamilyNegMask(unittest.TestCase):
    def test_mask_shape_and_semantics(self):
        rids = ["r1", "r1", "r2"]
        fams = ["F0", "F0", "F0"]
        mask = family_neg_mask(rids, fams, 3)
        self.assertEqual(mask.shape, (3, 3))
        # (0,2) and (2,0): same family, different receptor -> hard negative
        self.assertTrue(mask[0, 2])
        self.assertTrue(mask[2, 0])
        # same receptor is a positive, never hard
        self.assertFalse(mask[0, 1])
        self.assertFalse(mask[1, 0])
        self.assertFalse(mask[0, 0])

    def test_cross_family_never_hard(self):
        rids = ["r1", "r2", "r3"]
        fams = ["F0", "F1", "F0"]
        mask = family_neg_mask(rids, fams, 3)
        self.assertTrue(mask[0, 2])   # same family
        self.assertFalse(mask[0, 1])  # cross family


class TestSupConV2(unittest.TestCase):
    def test_clustered_lower_than_random(self):
        set_seed()
        embs = torch.randn(8, 16)
        rids = ["r1", "r1", "r2", "r2", "r3", "r3", "r4", "r4"]
        fams = ["F0", "F0", "F0", "F0", "F1", "F1", "F2", "F2"]
        loss = supcon_loss_v2(_unit(embs), rids, fams)
        self.assertTrue(math.isfinite(float(loss)))
        base = torch.randn(4, 16)
        clustered = torch.cat([base, base])
        clustered_rids = ["r1", "r2", "r3", "r4", "r1", "r2", "r3", "r4"]
        clustered_fams = ["F0", "F0", "F0", "F0", "F0", "F0", "F0", "F0"]
        loss2 = supcon_loss_v2(_unit(clustered), clustered_rids, clustered_fams)
        self.assertLess(float(loss2), float(loss))

    def test_family_weighting_penalizes_hard_negatives(self):
        # anchors r1 (2 copies) plus one hard negative r2 of the same family
        # sitting exactly on the r1 cluster: weighting must raise the loss
        d = torch.tensor([[1.0, 0.0]])
        embs = torch.cat([d, d, d])  # all three identical
        rids = ["r1", "r1", "r2"]
        fams = ["F0", "F0", "F0"]
        l0 = float(supcon_loss_v2(embs, rids, fams, temperature=0.07,
                                  family_neg_weight=0.0))
        l1 = float(supcon_loss_v2(embs, rids, fams, temperature=0.07,
                                  family_neg_weight=1.0))
        self.assertGreater(l1, l0)

    def test_anchor_without_positives_contributes_zero(self):
        # single anchor, no positive partner: loss is exactly 0
        embs = torch.randn(1, 8)
        l = float(supcon_loss_v2(_unit(embs), ["r1"], ["F0"]))
        self.assertEqual(l, 0.0)
        # two same-receptor anchors: positive contribution, finite loss
        embs2 = torch.randn(3, 8)
        l2 = float(supcon_loss_v2(_unit(embs2), ["r1", "r1", "r2"],
                                  ["F0", "F0", "F1"]))
        self.assertTrue(math.isfinite(l2) and l2 > 0.0)


class TestEmbedderV2(unittest.TestCase):
    def test_forward_shape_and_no_layernorm(self):
        set_seed()
        m = EmbedderV2(26, 16, emb_dim=128)
        self.assertFalse(any(isinstance(x, torch.nn.LayerNorm)
                             for x in m.modules()))
        h = torch.randn(12, 26)
        x = torch.randn(12, 3)
        ei = torch.randint(0, 12, (2, 20))
        ea = torch.randn(20, 16)
        bt = torch.zeros(12, dtype=torch.long)
        pm = torch.zeros(12, dtype=torch.bool)
        pm[4:8] = True
        out = m(h, x, ei, ea, bt, pm)
        # graph-level pooled embedding: 12 nodes in one graph -> (1, 128)
        self.assertEqual(out.shape, (1, 128))


class TestRashomonHelpers(unittest.TestCase):
    def test_spearman(self):
        a = np.array([3.0, 1.0, 2.0, 4.0])
        self.assertAlmostEqual(rashomon.spearman_rho(a, a), 1.0)
        mono = np.array([1.0, 2.0, 3.0, 4.0])
        self.assertAlmostEqual(rashomon.spearman_rho(mono, mono[::-1]), -1.0)

    def test_jaccard(self):
        self.assertEqual(rashomon.jaccard([1, 2, 3], [2, 3, 4]), 0.5)
        self.assertEqual(rashomon.jaccard([1, 2], [1, 2]), 1.0)
        self.assertEqual(rashomon.jaccard([1], [2]), 0.0)

    def test_flip_rate(self):
        # identical matrices -> no flips
        S = np.array([[0.9, 0.1], [0.2, 0.8]])
        pop = [(0, 1), (1, 0)]
        pairs = [(1, 2)]
        sims = {1: S, 2: S.copy()}
        self.assertAlmostEqual(
            rashomon.flip_rate(sims, pop, pairs)["mean"], 0.0)
        # relative ordering of the two pairs reverses -> everything flips
        # (instance-median thresholds are per-instance, so a monotone shift
        # would NOT flip; only a rank change does)
        S1 = np.array([[0.9, 0.8], [0.1, 0.5]])
        S2 = np.array([[0.5, 0.1], [0.8, 0.9]])
        sims = {1: S1, 2: S2}
        self.assertAlmostEqual(
            rashomon.flip_rate(sims, pop, pairs)["mean"], 1.0)

    def test_bootstrap_ci(self):
        vals = np.arange(100, dtype=float)
        ci = rashomon.bootstrap_ci(vals, [str(i) for i in range(100)], n=200)
        self.assertEqual(ci["n"], 100)
        self.assertLess(ci["ci95"][0], ci["mean"])
        self.assertLess(ci["mean"], ci["ci95"][1])

    def test_fold_rashomon_perfect_instances(self):
        # 4 graphs, 2 receptors (r0: graphs 0/1, r1: graphs 2/3), 1 family
        recs = [{"receptor_id": f"r{i // 2}", "family_id": "F0",
                 "complex_id": str(i), "domain": "crystal"}
                for i in range(4)]
        E = np.zeros((4, 2))
        E[0], E[1], E[2], E[3] = [1, 0], [1, 0], [0, 1], [0, 1]
        E = E / np.linalg.norm(E, axis=1, keepdims=True)
        S = E @ E.T
        sims = {s: S for s in (20261021, 1101, 4242)}
        embs = {s: E for s in (20261021, 1101, 4242)}
        rankings = {}
        for s in (20261021, 1101, 4242):
            rankings[s] = {}
            for q in range(4):
                g = [j for j in range(4) if j != q]
                order = sorted(g, key=lambda j: -S[q, j])
                rankings[s][str(q)] = {"gallery_order": order,
                                       "scores": [float(S[q, j]) for j in order]}
        metrics = {}
        for s in (20261021, 1101, 4242):
            metrics[s] = {"retrieval": {"auroc": 1.0, "rank1": 1.0},
                          "family": {"1": 1.0}, "domain": {"accuracy": 0.5}}
        with patch.object(rashomon, "SEEDS", [20261021, 1101, 4242]), \
                patch.object(rashomon, "SEED_PAIRS", [(20261021, 1101)]):
            fr = rashomon.fold_rashomon(0, [0, 1], recs, sims, embs,
                                        rankings, metrics)
        self.assertEqual(fr["ranking_stability"]["pooled"]["mean"], 1.0)
        self.assertEqual(fr["retrieval_agreement"]["rank1_agreement"]["mean"],
                         1.0)
        self.assertEqual(fr["pair_flip_rate"]["mean"], 0.0)
        self.assertEqual(fr["pooled"]["auroc"], 1.0)
        self.assertEqual(fr["embedding_agreement"]
                         ["same_graph_across_instances"]["mean"], 1.0)


class TestPlanArtifacts(unittest.TestCase):
    def test_followup_plan_committed_and_consistent(self):
        plan_dir = ROOT / "docs/fingerprint_retrieval_20261021"
        plan = json.loads((plan_dir / "followup_plan.json").read_text())
        self.assertEqual(plan["task_id"], "FOLLOWUP-AUROC-RASHOMON-20261021")
        self.assertEqual(plan["parent_task"], "STRUCT-FINGERPRINT-20261021")
        self.assertIn("pooled same/diff pair AUROC >= 0.9",
                      plan["success_criteria"])
        self.assertEqual(len(plan["rashomon"]["seeds"]), 3)
        md = (plan_dir / "FOLLOWUP_AUROC_RASHOMON_20261021.md").read_text()
        self.assertIn("预注册", md)
        self.assertIn("Rashomon", md)


if __name__ == "__main__":
    unittest.main()
