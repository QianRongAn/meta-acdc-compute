"""Negative cases for independent SPR source/identity/fold contracts."""
import copy,json,unittest
from meta_acdc.data.independent_spr import FREEZE,ROOT,author_keys,source_components,validate_partition,verify_freeze,connected_test_blocks
from meta_acdc.data.experimental_evidence import assert_reserved_absent

class IndependentSPRTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        read=lambda n:json.loads((FREEZE/n).read_text())
        cls.base=[read('frozen_test.json'),read('development_records.json'),read('test_study_sources.json'),read('development_study_sources.json'),json.loads((ROOT/'docs/sequence_followup_20261018/review_ledger.json').read_text())]
    def args(self):return copy.deepcopy(self.base)
    def test_valid_partition(self):
        x=validate_partition(*self.args());self.assertEqual((x['observations'],x['primary'],x['lineages'],x['research_source_components']),(15,14,7,3));self.assertEqual(x['development'],{'train':10,'validation':10})
    def test_original_connected_block_rule_retained_without_reassigning_folds(self):
        self.assertEqual(len(connected_test_blocks(self.base[0])),4)
        self.assertTrue(all(x['fold']=='external_test' for x in self.base[0]))
    def test_shared_coauthor_sources_merge_transitively(self):
        f=lambda p,n:dict(pmid=p,authors={'author':[{'fullName':x} for x in n]})
        self.assertEqual(source_components([f('a',['One A']),f('b',['Two B']),f('c',['One A','Two B'])]),[['a','b','c']])
    def test_surnames_alone_do_not_merge(self):
        self.assertNotEqual(author_keys({'authors':{'author':[{'fullName':'Chen G'}]}}),author_keys({'authors':{'author':[{'fullName':'Chen J'}]}}))
    def test_diacritics_and_orcid_aliases_do_not_hide_source_overlap(self):
        f=lambda name:dict(authors={'author':[dict(fullName=name,authorId={'type':'ORCID','value':'0000-0001-1234-5678'})]})
        self.assertEqual(author_keys(f('Dündar F')),author_keys(f('Dundar F')))
        self.assertIn('orcid:0000-0001-1234-5678',author_keys(f('Different Alias')))
    def test_shared_source_cannot_reenter_development(self):
        a=self.args();p=next(x for x in a[4] if x['evidence_id']=='ATLAS:717');a[1].append(copy.deepcopy(p['original_record']))
        with self.assertRaisesRegex(ValueError,'source overlap'):validate_partition(*a)
    def test_old_fold_cannot_be_reassigned(self):
        a=self.args();a[1][0]['fold']='validation' if a[1][0]['fold']=='train' else 'train'
        with self.assertRaisesRegex(ValueError,'fold or block'):validate_partition(*a)
    def test_outcome_dependent_development_drop_rejected(self):
        a=self.args();a[1].pop()
        with self.assertRaisesRegex(ValueError,'filter coverage'):validate_partition(*a)
    def test_test_training_flag_rejected(self):
        a=self.args();a[0][0]['training_allowed']=True
        with self.assertRaisesRegex(ValueError,'endpoint promotion'):validate_partition(*a)
    def test_no_clinical_label_promotion(self):
        a=self.args();a[0][0]['clinical_risk_label']=0
        with self.assertRaisesRegex(ValueError,'endpoint promotion'):validate_partition(*a)
    def test_censored_not_exact_regression_label(self):
        a=self.args();a[0][0]['qualifier']='>'
        with self.assertRaisesRegex(ValueError,'censored'):validate_partition(*a)
    def test_range_extrapolation_not_primary_metric(self):
        a=self.args();next(x for x in a[0] if x['metric_role']=='sensitivity_only')['metric_role']='primary'
        with self.assertRaisesRegex(ValueError,'extrapolation'):validate_partition(*a)
    def test_wrong_chain_rejected(self):
        a=self.args();a[0][-1]['cdr3_alpha']='CAVWRONGF'
        with self.assertRaisesRegex(ValueError,'paired variable'):validate_partition(*a)
    def test_unit_factor_mismatch_rejected(self):
        a=self.args();a[0][-1]['raw_unit']='nM'
        with self.assertRaisesRegex(ValueError,'KD conversion'):validate_partition(*a)
    def test_shared_V_does_not_inflate_lineages(self):
        a=self.args();new=[x for x in a[0] if x['pmid']=='39500929'];self.assertEqual(len(new),8);self.assertEqual({x['lineage'] for x in new},{'S2','S8'});self.assertEqual({x['trav'] for x in new},{'TRAV21'})
    def test_new_family_reserved_in_global_training_guard(self):
        row=next(x for x in self.base[0] if x['receptor_name']=='S8-WT')
        with self.assertRaisesRegex(ValueError,'reserved external receptor family'):assert_reserved_absent([row])
    def test_complete_freeze(self):
        seal=json.loads((FREEZE/'manifest.json').read_text())
        required=list(seal['reference_files'])+list(json.loads((FREEZE/'required_sources.json').read_text()))
        if any(not (ROOT/name).exists() for name in required):
            self.skipTest('local primary sources/historical snapshots absent; full integration validation not claimed')
        self.assertFalse(verify_freeze()['training_allowed'])

if __name__=='__main__':unittest.main()
