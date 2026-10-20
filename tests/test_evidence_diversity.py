import copy,hashlib,json,tempfile,unittest
from pathlib import Path
from meta_acdc.data.evidence_completion import CURATION
from meta_acdc.data.evidence_diversity import validate_kras_reference,validate_partial_spr
from meta_acdc.data.experimental_evidence import assert_reserved_absent
ROOT=Path(__file__).resolve().parents[1]
class EvidenceDiversity(unittest.TestCase):
    def load(self,n):return json.loads((CURATION/n).read_text())
    def kras(self):return self.load('KRAS10mer_published_protein_reference.json')
    def partial(self):return self.load('external_differentV_candidates.json')
    def test_full_protein_reference_is599_not_DNA(self):
        x=self.kras();self.assertTrue(validate_kras_reference(x,check_sources=False));self.assertEqual((len(x['alpha']),len(x['beta']),len(x['linker'])),(269,303,27));self.assertIsNone(x['full_DNA'])
    def test_backtranslation_cannot_be_called_published_insert(self):
        x=self.kras();x['full_DNA']='ATG'*599
        with self.assertRaisesRegex(ValueError,'DNA proof'):validate_kras_reference(x,check_sources=False)
    def test_reference_cannot_be_called_clinical_lot(self):
        x=self.kras();x['clinical_exact_insert_verified']=True
        with self.assertRaisesRegex(ValueError,'clinical'):validate_kras_reference(x,check_sources=False)
    def test_wrong_order_rejected(self):
        x=self.kras();x['chain_order']='beta-SEQ58-alpha'
        with self.assertRaisesRegex(ValueError,'order'):validate_kras_reference(x,check_sources=False)
    def test_transcription_typo_rejected(self):
        x=self.kras();x['alpha']='A'+x['alpha'][1:]
        with self.assertRaisesRegex(ValueError,'transcription'):validate_kras_reference(x,check_sources=False)
    def test_differentV_candidates_are_not_full_sequence_core(self):
        self.assertEqual(validate_partial_spr(self.partial(),check_sources=False),dict(reserved_partial_identity=2,historically_exposed=1,model_ready=0));self.assertEqual({x['trbv'] for x in self.partial()[:2]},{'TRBV20-1'})
    def test_germline_imputation_rejected(self):
        x=self.partial();x[0]['alpha_variable']='FAKE'
        with self.assertRaisesRegex(ValueError,'imputation'):validate_partial_spr(x,check_sources=False)
    def test_KS1_history_cannot_be_hidden(self):
        x=self.partial();x[2]['scope']='reserved_partial_identity_not_model_ready'
        with self.assertRaisesRegex(ValueError,'historical'):validate_partial_spr(x,check_sources=False)
    def test_NB_cannot_be_binary_negative_or_infinite_KD(self):
        for key,val in [('binary_binding_label',0),('kd_um',float('inf')),('max_tested_um',500)]:
            x=self.partial();x[0]['counterligand'][key]=val
            with self.assertRaisesRegex(ValueError,'NB promoted'):validate_partial_spr(x,check_sources=False)
    def test_pending_CDR3_reservation_not_bypassable_by_empty_registry(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'empty.json';p.write_text('{"receptors":[]}')
            for x in self.partial()[:2]:
                with self.assertRaisesRegex(ValueError,'reserved'):assert_reserved_absent([x],p)
    def test_primary_author_temperature_and_pair_link_explicit(self):
        rows={x['evidence_id']:x for x in self.load('review_ledger.json')}
        for k in ['ATLAS:747','ATLAS:751']:
            x=rows[k];self.assertEqual(x['temperature_K'],298.15);self.assertIn('primary-author thesis',x['temperature_evidence_kind']);self.assertIn('115',x['method_locator']);self.assertTrue(x['training_allowed'])
        self.assertEqual(rows['ATLAS:747']['parent_block'],rows['ATLAS:751']['parent_block'])
    def test_D30_source_values_and_folds(self):
        rows={x['evidence_id']:x for x in self.load('review_ledger.json')}
        for key,val in [('ATLAS:727',136.),('ATLAS:728',142.5)]:
            x=rows[key];self.assertEqual(x['primary_kd_um'],val);self.assertEqual(x['parent_fold'],'validation');self.assertEqual(x['temperature_K'],298.15)
    def test_core_not_enlarged_with_partial_identity(self):
        old=json.loads((ROOT/'docs/spr_followup_20261018/external_classical_panel.json').read_text());self.assertEqual(self.load('external_classical_panel.json'),old);rows=self.load('review_ledger.json');test=[x for x in rows if x['training_allowed'] and x['parent_fold']=='test'];self.assertEqual(len(test),3);self.assertEqual(len({x['parent_block'] for x in test}),2)
    def test_LC13_direct_KD_not_kinetic_derived_column(self):
        rows={x['evidence_id']:x for x in self.load('review_ledger.json')}
        self.assertEqual(rows['ATLAS:145']['primary_kd_um'],1.54);self.assertEqual(rows['ATLAS:148']['primary_kd_um'],49.)
        for key in ['ATLAS:145','ATLAS:148']:self.assertTrue(rows[key]['training_allowed']);self.assertEqual(rows[key]['temperature_K'],298.15)
    def test_numeric_only_supplements_do_not_release_missing_conditions(self):
        rows={x['evidence_id']:x for x in self.load('review_ledger.json')}
        for key,val in [('ATLAS:514',4.58),('ATLAS:515',1.1),('ATLAS:532',8.89),('ATLAS:537',15.07)]:
            x=rows[key];self.assertEqual(x['primary_kd_um'],val);self.assertTrue(x['primary_supplement_retrieved']);self.assertFalse(x['training_allowed']);self.assertTrue(x['numeric_source_sha256'])
    def test_local_primary_checksums_and_reference(self):
        if not (ROOT/self.kras()['sources'][0]['path']).exists():self.skipTest('third-party raw sources local-only')
        validate_kras_reference(self.kras());validate_partial_spr(self.partial())
if __name__=='__main__':unittest.main()
