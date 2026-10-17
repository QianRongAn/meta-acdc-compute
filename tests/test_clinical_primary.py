import copy,json,tempfile,unittest
from pathlib import Path
from meta_acdc.data.clinical_primary import translate_fragment,foy_receptors,validate_bundle
from meta_acdc.data.experimental_evidence import reserved_hits,assert_reserved_absent
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'docs/clinical_primary_20261018'

class ClinicalPrimary(unittest.TestCase):
    def bundle(self):return json.loads((D/'clinical_bundle.json').read_text())
    def test_real_bundle_has_no_receptor_clinical_labels(self):self.assertTrue(validate_bundle(self.bundle()))
    def test_sequence_conflict_is_quarantined(self):
        b=self.bundle();x=next(x for x in b['receptors'] if x['receptor_name']=='Foy:TCR919')
        self.assertFalse(x['alpha_fragment']['cdr3_matches']);self.assertFalse(x['beta_fragment']['cdr3_matches'])
        x['identity_status']='variable_pair_verified'
        with self.assertRaisesRegex(ValueError,'quarantine'):validate_bundle(b)
    def test_cdr3_is_checked_after_dna_translation(self):
        self.assertTrue(translate_fragment('TGTGCTTCTTCT','CASS')['cdr3_matches'])
        self.assertFalse(translate_fragment('TGTGCTTCTTCT','CASR')['cdr3_matches'])
    def test_terminal_partial_codon_is_retained(self):
        x=translate_fragment('TGTGCTTCTTCTA','CASS');self.assertEqual(x['trailing_nt'],'A');self.assertEqual(x['amino_acids'],'CASS')
    def test_internal_stop_rejected(self):
        with self.assertRaisesRegex(ValueError,'stop'):translate_fragment('TGTTAGGCT','CASS')
    def test_unknown_base_rejected(self):
        with self.assertRaises(ValueError):translate_fragment('TGTNCT','CASS')
    def test_sequence_poisoning_is_detected(self):
        b=self.bundle();b['receptors'][0]['alpha_variable']='A'*110
        with self.assertRaisesRegex(ValueError,'changed'):validate_bundle(b)
    def test_no_promoting_patient_event_to_receptor_label(self):
        b=self.bundle();b['receptors'][0]['clinical_risk_label']=0
        with self.assertRaisesRegex(ValueError,'risk labels'):validate_bundle(b)
    def test_multireceptor_join_poisoning_rejected(self):
        b=self.bundle();b['patients'][0]['receptors'].pop()
        with self.assertRaises(ValueError):validate_bundle(b)
    def test_missing_assessment_is_preserved(self):
        p=next(x for x in self.bundle()['patients'] if x['patient_id']=='0417');self.assertEqual(p['response_table3'],'No post-baseline assessment')
    def test_no_per_receptor_causal_labels(self):
        b=self.bundle();b['patients'][0]['per_receptor_causality_established']=True
        with self.assertRaisesRegex(ValueError,'causal'):validate_bundle(b)
    def test_e7_published_insert_has_both_cdrs_and_two_terminal_stops(self):
        x=next(x for x in self.bundle()['receptors'] if x['receptor_name']=='HPV16-E7-optimized')
        self.assertEqual(len(x['dna']),1800);self.assertEqual(len(x['protein']),598);self.assertEqual(x['terminal_stop_codons'],2)
        self.assertIn(x['cdr3_alpha'],x['protein']);self.assertIn(x['cdr3_beta'],x['protein']);self.assertFalse(x['clinical_lot_sequence_verified'])
    def test_partial_clinical_cdr3_reservation_is_enforced(self):
        ref=dict(receptor_name='clinical',alpha_variable='',beta_variable='',cdr3_alpha='CASSLGWRGGRYNEQF')
        self.assertEqual(reserved_hits(dict(cdr3_alpha='CASSLGWRGGRYNEQF'),[ref]),['clinical'])
    def test_explicit_empty_registry_cannot_bypass_clinical_guard(self):
        x=self.bundle()['receptors'][0]
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'registry.json';p.write_text('{"receptors": []}')
            with self.assertRaisesRegex(ValueError,'reserved external'):assert_reserved_absent([x],p)
    def test_new_affinity_dataset_keeps_old_fold_membership(self):
        old={x['evidence_id']:x for x in json.loads((ROOT/'docs/database_binding/affinity_dataset.json').read_text())}
        for x in json.loads((D/'affinity_dataset.json').read_text()):
            self.assertEqual((x['fold'],x['split_block']),(old[x['evidence_id']]['fold'],old[x['evidence_id']]['split_block']))
    def test_removed_parent_connected_blocks_cannot_bridge_holdout(self):
        old=json.loads((ROOT/'docs/database_binding/affinity_dataset.json').read_text())
        new=json.loads((D/'affinity_dataset.json').read_text());ids={x['evidence_id'] for x in new}
        blocked={x['split_block'] for x in old if x['evidence_id'] not in ids}
        self.assertFalse(blocked & {x['split_block'] for x in new})
    def test_affinity_does_not_include_clinical_candidates(self):
        refs=json.loads((D/'clinical_reservations.json').read_text())['receptors']
        for x in json.loads((D/'affinity_dataset.json').read_text()):self.assertEqual(reserved_hits(x,refs),[])

if __name__=='__main__':unittest.main()
