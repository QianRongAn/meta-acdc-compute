"""Release/fit/acceptance isolation and statistical unit safeguards."""
import copy,json,tempfile,unittest,math
from pathlib import Path
from unittest.mock import patch
from meta_acdc.data.affinity_training_release import write,verify
from meta_acdc.models.frozen_affinity_baseline import fit,verify_lock,weighted_metrics,block_weights
from meta_acdc.models.affinity_acceptance import hierarchical_weights,metrics,compare,acceptance_rule,claim_test

class FrozenAffinityPipeline(unittest.TestCase):
    def test_source_and_lineage_weights_not_observation_count(self):
        rows=[dict(research_source='a',lineage='x')]*7+[dict(research_source='a',lineage='y')]+[dict(research_source='b',lineage='z')]*2
        w=hierarchical_weights(rows);self.assertAlmostEqual(sum(w[:7]),.25);self.assertAlmostEqual(w[7],.25);self.assertAlmostEqual(sum(w[8:]),.5)
    def test_repeated_identical_measurement_does_not_inflate_metric(self):
        a=[dict(research_source='s',lineage='a',pkd=1),dict(research_source='t',lineage='b',pkd=5)]
        first=compare(a,[2,4],3);second=compare([a[0]]*8+[a[1]],[2]*8+[4],3)
        self.assertEqual(first['model'],second['model']);self.assertEqual(first['train_mean'],second['train_mean'])
    def test_independent_math_metrics_match_hand_calculation(self):
        self.assertEqual(metrics([0,2],[1,1],[1,1]),{'rmse_pkd':1.,'mae_pkd':1.})
        self.assertAlmostEqual(metrics([0,4],[0,0],[3,1])['rmse_pkd'],2)
    def test_fit_and_independent_metric_implementations_agree(self):
        a=metrics([0,4,2],[1,2,3],[.5,.25,.25]);b=weighted_metrics([0,4,2],[1,2,3],[.5,.25,.25])
        for k in a:self.assertAlmostEqual(a[k],b[k],places=12)
    def test_invalid_predictions_rejected(self):
        with self.assertRaises(ValueError):metrics([0],[math.nan],[1])
        with self.assertRaises(ValueError):metrics([0],[1],[0])
    def test_source_sensitivity_failure_blocks_performance_pass(self):
        main=dict(rmse_improvement=.1,mae_improvement=.1)
        self.assertFalse(acceptance_rule(main,{'a':dict(rmse_improvement=-.01)}))
        self.assertTrue(acceptance_rule(main,{'a':dict(rmse_improvement=.01)}))
    def test_primary_mae_must_also_improve(self):
        self.assertFalse(acceptance_rule(dict(rmse_improvement=1,mae_improvement=-.01),{'a':dict(rmse_improvement=1)}))
    def test_exclusive_test_claim_preserves_first_access(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'claim.json';claim_test(p,{'model':'first'});before=p.read_bytes()
            with self.assertRaisesRegex(ValueError,'already claimed'):claim_test(p,{'model':'second'})
            self.assertEqual(before,p.read_bytes())
    def test_old_blocked_release_still_cannot_fit(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);write(p/'manifest.json',dict(task='experimental_SPR_pKD_regression',training_allowed=False))
            with self.assertRaisesRegex(ValueError,'not released'):verify(p,for_fit=True)
    def test_post_test_fit_is_blocked_before_reading_targets(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);ledger=p/'test.json';write(ledger,{})
            write(p/'manifest.json',dict(task='experimental_SPR_pKD_regression',training_allowed=True,test_access_ledger=str(ledger)))
            with self.assertRaisesRegex(ValueError,'test already accessed'):verify(p,for_fit=True)
    def test_duplicate_fit_dispatch_is_blocked(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);ledger=p/'fit.json';write(ledger,{})
            write(p/'manifest.json',dict(task='experimental_SPR_pKD_regression',training_allowed=True,test_access_ledger=str(p/'test.json'),fit_access_ledger=str(ledger)))
            with self.assertRaisesRegex(ValueError,'duplicate dispatch'):verify(p,for_fit=True)
    def test_training_and_lock_read_development_only(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t);release=r/'release';release.mkdir();run=r/'run';write(release/'manifest.json',{})
            write(release/'plan.json',dict(alpha_grid=[1.,10.,100.]))
            rows=[]
            for fold in ['train','validation']:
                for i in range(6):
                    rows.append(dict(evidence_id=fold+str(i),fold=fold,split_block=fold+str(i//2),alpha_variable='ACDEFGHIKLMNPQRSTVWY'*5,beta_variable='QWERTYIPASDFGHKLCVNM'*5,peptide='AAAAAAA' if fold=='train' else 'WWYYYYY',mhc_allele='HLA-A*02:01',primary_assay={'instrument':'Biacore T200'},temperature_K=298.15,pkd=4+i/10))
            write(release/'development.json',rows);m=dict(task='experimental_SPR_pKD_regression',code_files={},fit_access_ledger=str(r/'fit_claim.json'))
            from meta_acdc.data.affinity_training_release import read as actual_read
            reads=[]
            def monitored(path):
                reads.append(Path(path).name)
                if 'test' in Path(path).name or 'frozen' in Path(path).name:raise AssertionError('heldout read during fit')
                return actual_read(path)
            with patch('meta_acdc.models.frozen_affinity_baseline.verify',return_value=m),patch('meta_acdc.models.frozen_affinity_baseline.read',side_effect=monitored):
                report=fit(release,run);verify_lock(release,run)
                model=actual_read(run/'model.json');self.assertEqual(model['fit_samples'],6);self.assertNotIn('p2:YY',model['feature_names']);self.assertFalse(report['test_labels_loaded'])
                self.assertEqual(report['selected_alpha'],min(report['candidate_selection'],key=lambda x:(x['validation']['rmse_pkd'],x['alpha']))['alpha'])
                with self.assertRaisesRegex(ValueError,'duplicate training'):fit(release,run)
                model['intercept']+=1;(run/'model.json').write_text(json.dumps(model))
                with self.assertRaisesRegex(ValueError,'locked artifact changed'):verify_lock(release,run)
            self.assertNotIn('frozen_test.json',reads)

if __name__=='__main__':unittest.main()
