"""Frozen experiment orchestration: tiny end-to-end exercise, not scientific results."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import torch
from tests.test_fingerprint_review import TinyEmbedder, data_fixture
from meta_acdc.models.train_fingerprint import group_folds
from meta_acdc.models import train_fingerprint_v2 as trainer
from scripts import run_fingerprint_retrain as runner

class FrozenExperimentChecks(unittest.TestCase):
    def test_complete_tiny_run_includes_control_replay_and_paired_evaluation(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);data=data_fixture();torch.save(data,root/'data.pt')
            families=[r['family_id'] for r in data['records']]
            folds=group_folds(families,3)
            fa=dict(per_graph_fold=[next(i for i,fs in enumerate(folds) if f in fs) for f in families],
                    graph_records=[r['complex_id'] for r in data['records']])
            (root/'folds.json').write_text(json.dumps(fa))
            sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
            protocol=dict(dataset=str(root/'data.pt'),dataset_sha256=sha(root/'data.pt'),
                          fold_assignment=str(root/'folds.json'),fold_assignment_sha256=sha(root/'folds.json'),
                          source_sha256={},n_folds=3,seeds=[20261021,1101,4242],arms={'true':'true','shuffled':'train-shuffle'},
                          recipe_cli={'epochs':2,'folds':3,'batch-size':8,'emb-dim':8,'cpu-threads':1})
            (root/'protocol.json').write_text(json.dumps(protocol))
            (root/'protocol.sha256').write_text(sha(root/'protocol.json'))
            argv=['run','--protocol',str(root/'protocol.json'),'--out-dir',str(root/'run')]
            with patch('sys.argv',argv), patch.object(trainer,'EmbedderV2',TinyEmbedder), contextlib.redirect_stdout(io.StringIO()):
                runner.main()
            status=json.loads((root/'run/status.json').read_text())
            self.assertEqual(status['state'],'complete')
            self.assertEqual(status['n_main_instances'],18)
            result=json.loads((root/'run/evaluation.json').read_text())
            self.assertTrue(result['replay']['weights_bitwise_equal'])
            self.assertTrue(result['replay']['embeddings_bitwise_equal'])
            self.assertIn('same_family',result['model_minus_baselines']['true'])
            self.assertFalse(result['practical_or_biological_gain_claim_allowed'])
            changes=[x['n_changed_receptor_labels'] for x in result['permutation_checks'] if x['arm']=='shuffled']
            self.assertTrue(all(x>0 for x in changes))
            self.assertEqual(sha(root/'data.pt'),protocol['dataset_sha256'])
            (root/'protocol.json').write_text('{}')
            with self.assertRaisesRegex(ValueError,'lock mismatch'):
                runner.validate_protocol(root/'protocol.json')

if __name__ == '__main__':
    unittest.main()
