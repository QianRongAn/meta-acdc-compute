"""Regression checks for the 2026-10-23 fingerprint review."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import torch
from torch import nn
from meta_acdc.models.train_fingerprint import batch_graphs
from meta_acdc.models import train_fingerprint_v2 as v2
from scripts import rashomon_fingerprint as rash
from scripts.fingerprint_shortcut_audit import query_auc

torch.set_num_threads(1)

class TinyEmbedder(nn.Module):
    def __init__(self, node_dim, edge_dim, emb_dim):
        super().__init__()
        self.proj = nn.Linear(node_dim, emb_dim)
    def forward(self, h, x, ei, ea, batch, pm):
        z = self.proj(h)
        return torch.stack([z[batch == i].mean(0) for i in range(int(batch.max()) + 1)])

def data_fixture():
    rng = np.random.default_rng(7)
    graphs, records = [], []
    for i in range(12):
        h = rng.normal(size=(4, 26)).astype('float32')
        h[:, 24] = [0, 0, 1, 1]
        graphs.append((h, rng.normal(size=(4, 3)).astype('float32'),
                       np.array([[0, 1], [1, 0], [2, 3], [3, 2]]),
                       rng.normal(size=(4, 16)).astype('float32')))
        records.append(dict(family_id=f'F{i//4}', receptor_id=f'R{i//2}',
                            complex_id=str(i), domain='crystal' if i % 2 else 'predicted',
                            peptide_mask=[False, False, True, True]))
    return dict(graphs=graphs, records=records)

class BatchChecks(unittest.TestCase):
    def test_offsets_and_batch_invariance(self):
        d = data_fixture()
        v2.set_seed(10)
        m = v2.EmbedderV2(26, 16, emb_dim=8).eval()
        joint = batch_graphs(d['graphs'], d['records'], [0, 1], 'cpu', m)
        ei = joint[2]
        self.assertTrue(torch.equal(ei[:, 4:], torch.tensor([[4,5,6,7],[5,4,7,6]])))
        self.assertTrue(torch.equal(joint[4][ei[0]], joint[4][ei[1]]))
        with torch.no_grad():
            both = m(*joint[:6])
            for i in range(2):
                single = batch_graphs(d['graphs'], d['records'], [i], 'cpu', m)
                torch.testing.assert_close(both[i], m(*single[:6])[0], atol=2e-6, rtol=2e-5)
    def test_invalid_local_edge_rejected(self):
        d = data_fixture()
        d['graphs'][0][2][0, 0] = 4
        with self.assertRaisesRegex(ValueError, 'own node'):
            batch_graphs(d['graphs'], d['records'], [0], 'cpu', None)

class InstanceChecks(unittest.TestCase):
    def test_full_only_fold_and_resume_are_equal_and_bound(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            torch.save(data_fixture(), p/'data.pt')
            common = ['--dataset',str(p/'data.pt'),'--epochs','2','--batch-size','8',
                      '--folds','3','--seeds','7','--device','cpu','--emb-dim','8']
            def run(name, extra=()):
                with patch.object(v2, 'EmbedderV2', TinyEmbedder), contextlib.redirect_stdout(io.StringIO()):
                    return v2.main(common + ['--out-dir',str(p/name)] + list(extra))
            run('full')
            run('isolated', ['--only-fold','1'])
            run('resumed', ['--only-fold','0'])
            run('resumed', ['--resume'])
            for name in ['isolated', 'resumed']:
                a = torch.load(p/'full/model_run_v2/fold1_seed7.pt', weights_only=False)['model']
                b = torch.load(p/name/'model_run_v2/fold1_seed7.pt', weights_only=False)['model']
                for k in a:
                    self.assertTrue(torch.equal(a[k], b[k]), k)
                np.testing.assert_array_equal(np.load(p/'full/model_run_v2/emb_fold1_seed7.npy'),
                                              np.load(p/name/'model_run_v2/emb_fold1_seed7.npy'))
            with self.assertRaisesRegex(ValueError, 'non-resume'):
                run('full')
            with (p/'resumed/model_run_v2/sim_fold1_seed7.npy').open('ab') as f:
                f.write(b'tampered')
            with self.assertRaisesRegex(ValueError, 'modified'):
                run('resumed', ['--resume'])
            (p/'legacy').mkdir()
            (p/'legacy/config.json').write_text('{}')
            with self.assertRaisesRegex(ValueError, 'legacy'):
                run('legacy', ['--resume'])

class ShortcutChecks(unittest.TestCase):
    def test_same_family_keeps_positives_and_reports_missing_classes(self):
        records = [dict(receptor_id='R',family_id='A'),
                   dict(receptor_id='R',family_id='B'),
                   dict(receptor_id='X',family_id='A'),
                   dict(receptor_id='Y',family_id='C')]
        scores = np.array([[1,0.9,0.1,0.2],[0.9,1,0.2,0.3],
                           [0.1,0.2,1,0.4],[0.2,0.3,0.4,1]])
        result, rows = query_auc(scores,records,[0,1,2,3],'same_family')
        self.assertEqual(rows[0],1.0)
        self.assertTrue(np.isnan(rows[1:]).all())
        self.assertEqual(result['n_undefined_queries'],3)
        self.assertEqual(result['n_pairs_defined_queries'],2)

class ClusterChecks(unittest.TestCase):
    def test_duplicate_seed_pairs_do_not_narrow_interval(self):
        a = np.array([[0.1], [0.3], [0.8], [0.9]])
        x = rash.bootstrap_ci(a, ['A','A','B','C'])
        y = rash.bootstrap_ci(np.repeat(a, 12, axis=1), ['A','A','B','C'])
        np.testing.assert_allclose(x['ci95'], y['ci95'], atol=1e-15, rtol=0)
        self.assertAlmostEqual(x['mean'], (0.2+0.8+0.9)/3)
    def test_one_family_and_undefined_are_explicit(self):
        x = rash.bootstrap_ci(np.array([[np.nan, np.nan], [0.5, np.nan]]), ['A','A'])
        self.assertEqual(x['ci95'], [None,None])
        self.assertEqual(x['n_undefined_values'], 3)
        self.assertEqual(x['n'], 1)
    def test_alignment_uses_only_training_graphs(self):
        recs = [dict(receptor_id=str(i//2), family_id=str(i//2)) for i in range(6)]
        rng = np.random.default_rng(8)
        base = rng.normal(size=(6,3)); base /= np.linalg.norm(base,axis=1,keepdims=True)
        embs = {s:base.copy() for s in rash.SEEDS}
        sims = {s:e@e.T for s,e in embs.items()}
        ranks = {s:{str(q):{'gallery_order':sorted([i for i in range(6) if i!=q],
                              key=lambda i:-sims[s][q,i])} for q in [0,1]} for s in rash.SEEDS}
        metrics = {s:{'retrieval':{'auroc':0.5,'rank1':0.0}} for s in rash.SEEDS}
        svd = np.linalg.svd
        with patch.object(np.linalg,'svd',wraps=svd) as spy:
            result = rash.fold_rashomon(0,[0,1],recs,sims,embs,ranks,metrics)
        for call in spy.call_args_list:
            np.testing.assert_allclose(call.args[0], base[2:].T@base[2:])
        self.assertEqual(result['embedding_agreement_aligned']['fit_graph_indices'],[2,3,4,5])
        self.assertEqual(result['embedding_agreement_aligned']['evaluation_graph_indices'],[0,1])
        self.assertAlmostEqual(result['pooled']['auroc'],result['embedding_agreement_aligned']['pooled_auroc'])

if __name__ == '__main__':
    unittest.main()
