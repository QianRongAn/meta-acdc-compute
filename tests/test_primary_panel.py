"""Cross-study curation and whole-family/donor reservation regression tests."""
import copy,json,unittest
from unittest.mock import patch
from pathlib import Path
from meta_acdc.data.primary_panel import validate_measurements,reservation_exclusions
ROOT=Path(__file__).resolve().parents[1]

class PrimaryPanel(unittest.TestCase):
    def setUp(self):self.rows=json.loads((ROOT/'docs/evidence/primary_measurements.json').read_text())
    def check_bad(self,mutate,pattern):
        rows=copy.deepcopy(self.rows);mutate(rows)
        with self.assertRaisesRegex(ValueError,pattern):validate_measurements(rows,ROOT)
    def test_real_source_snapshots(self):
        # Network-free; downloaded raw sources are intentionally untracked.
        if not (ROOT/self.rows[0]['sources'][0]['path']).exists():self.skipTest('local primary source snapshots not available')
        result=validate_measurements(self.rows,ROOT);self.assertEqual(result['studies'],4);self.assertEqual(result['endpoints'],{'tcr_binding_spr':21,'tcr_cellular_tnf_alpha':2})
    def test_parent_and_engineered_wt_are_different_observations(self):
        wt=[x for x in self.rows if x['peptide']=='VVVGAGGVGK'];self.assertEqual({x['receptor_name']:x['observed_outcome'] for x in wt},{'JDI':'not_detected','JDIa41b1':'measured_binding'})
        self.assertNotEqual(wt[0]['alpha_variable'],wt[1]['alpha_variable']);self.assertEqual(len({x['lineage'] for x in wt}),1)
    def test_positive_only_receptors_do_not_gain_fabricated_wt_negatives(self):
        spr=[x for x in self.rows if x['pmid']=='35484264' and x['endpoint']=='tcr_binding_spr'];self.assertEqual(len(spr),2);self.assertTrue(all(x['observed_outcome']=='measured_binding' for x in spr))
    def test_all_external_training_labels_are_empty(self):self.assertTrue(all(x['training_label'] is None for x in self.rows))
    def test_whole_family_and_graft_donors_are_removed(self):
        x=self.rows[0];reg=[x];structures=[{'complex_id':'a','family_id':'F','alpha':x['alpha_variable'],'beta':x['beta_variable']},{'complex_id':'b','family_id':'F','alpha':'A'*80,'beta':'G'*80},{'complex_id':'c','family_id':'G','alpha':'A'*80,'beta':'G'*80}]
        records=[{'graph_id':'a','family_id':'F','donor':None},{'graph_id':'b','family_id':'F','donor':None},{'graph_id':'c.graft','family_id':'G','donor':'b'},{'graph_id':'c','family_id':'G','donor':None}]
        keep,reject,complexes=reservation_exclusions(structures,records,reg);self.assertEqual(keep,[3]);self.assertEqual(complexes,{'a','b'});self.assertEqual(reject[-1]['reason'],'reserved_graft_donor')
    def test_censored_value_is_not_exact_affinity(self):
        if not (ROOT/self.rows[0]['sources'][0]['path']).exists():self.skipTest('local primary source snapshots not available')
        self.check_bad(lambda rows:rows[0].update(qualifier='>'),'invalid exact affinity')
    def test_source_tampering_is_rejected(self):
        with patch('meta_acdc.data.primary_panel.sha256',return_value='1'*64):
            self.check_bad(lambda rows:rows[0]['sources'][0].update(sha256='0'*64),'source changed')
    def test_training_label_is_rejected(self):self.check_bad(lambda rows:rows[0].update(training_label=1),'training labels')

if __name__=='__main__':unittest.main()
