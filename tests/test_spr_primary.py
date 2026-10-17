import copy,json,tempfile,unittest
from pathlib import Path
from meta_acdc.data.spr_primary import validate_review
from meta_acdc.models.affinity_baseline import run
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'docs/spr_primary_20261018'
class PrimarySPR(unittest.TestCase):
    def rows(self):return json.loads((D/'review_ledger.json').read_text())
    def parent(self):return json.loads((ROOT/'docs/clinical_primary_20261018/affinity_dataset.json').read_text())
    def row(self,rows,key):return next(x for x in rows if x['evidence_id']=='ATLAS:'+key)
    def test_all_51_have_individual_dispositions(self):self.assertTrue(validate_review(self.rows(),self.parent(),check_sources=False))
    def test_missing_row_rejected(self):
        with self.assertRaisesRegex(ValueError,'coverage'):validate_review(self.rows()[:-1],self.parent(),check_sources=False)
    def test_duplicate_row_rejected(self):
        x=self.rows();x[-1]=copy.deepcopy(x[0])
        with self.assertRaises(ValueError):validate_review(x,self.parent(),check_sources=False)
    def test_unknown_status_rejected(self):
        x=self.rows();x[0]['status']='approved_somehow'
        with self.assertRaisesRegex(ValueError,'status'):validate_review(x,self.parent(),check_sources=False)
    def test_raw_database_record_preserved(self):
        x=self.rows();x[0]['original_record']['pkd']=7
        with self.assertRaisesRegex(ValueError,'original'):validate_review(x,self.parent(),check_sources=False)
    def test_no_test_to_train_reassignment(self):
        x=self.rows();x[0]['parent_fold']='test' if x[0]['parent_fold']!='test' else 'train'
        with self.assertRaisesRegex(ValueError,'split'):validate_review(x,self.parent(),check_sources=False)
    def test_pn_wrong_peptide_number_is_quarantined(self):
        x=self.rows();v=self.row(x,'1159');self.assertEqual(v['original_record']['measurement']['kd_um'],420);self.assertEqual(v['primary_kd_um'],6.1)
        v['training_allowed']=True
        with self.assertRaisesRegex(ValueError,'quarantine'):validate_review(x,self.parent(),check_sources=False)
    def test_octet_cannot_be_relabelled_verified_spr(self):
        x=self.rows();v=self.row(x,'772');v['status']='primary_verified_soluble_SPR';v['training_allowed']=True
        with self.assertRaisesRegex(ValueError,'SPR'):validate_review(x,self.parent(),check_sources=False)
    def test_unverified_temperature_cannot_train(self):
        x=self.rows();v=self.row(x,'747');v['status']='primary_verified_soluble_SPR';v['training_allowed']=True
        with self.assertRaisesRegex(ValueError,'temperature'):validate_review(x,self.parent(),check_sources=False)
    def test_47picomolar_is_not_47nanomolar(self):self.assertAlmostEqual(self.row(self.rows(),'747')['primary_kd_um'],.000047)
    def test_equilibrium_is_not_replaced_with_kinetic_kd(self):
        self.assertEqual(self.row(self.rows(),'1130')['primary_value'],.43);self.assertEqual(self.row(self.rows(),'1134')['primary_value'],.41)
    def test_binary_label_promotion_rejected(self):
        x=self.rows();x[0]['binary_binding_label']=1
        with self.assertRaisesRegex(ValueError,'endpoint'):validate_review(x,self.parent(),check_sources=False)
    def test_primary_value_corruption_rejected(self):
        x=self.rows();self.row(x,'1127')['primary_value']=6530
        with self.assertRaisesRegex(ValueError,'value'):validate_review(x,self.parent(),check_sources=False)
    def test_p53_all39_fragments_concordant_only_one_engineered_outcome(self):
        x=json.loads((D/'p53_paired_sequence_ledger.json').read_text());self.assertEqual(len(x),39);self.assertTrue(all(a['cdr3_pair_verified'] for a in x))
        linked=[a for a in x if a['clinical_patient_id']];self.assertEqual(len(linked),1);self.assertEqual(linked[0]['tcr_name'],'4196-AV6/BV11');self.assertEqual(linked[0]['clinical_patient_id'],'4349')
    def test_no_detection_is_not_exact_kd_or_binary_label(self):
        for x in json.loads((D/'conditional_nondetection_controls.json').read_text()):self.assertIsNone(x['kd_um']);self.assertIsNone(x['binary_binding_label'])
    def test_frozen_clinical18_artifact_unchanged(self):
        import hashlib
        rel='docs/clinical_primary_20261018/reserved_acceptance_candidates.json'
        old=json.loads((ROOT/'docs/clinical_primary_20261018/review_manifest.json').read_text())['artifacts'][rel]
        self.assertEqual(hashlib.sha256((ROOT/rel).read_bytes()).hexdigest(),old)
    def test_c4_comparator_is_not_receptor_negative(self):
        row=next(x for x in json.loads((D/'clinical_link_followup.json').read_text()) if x['receptor']=='WT1-TCR-C4')
        self.assertEqual(row['treated_patients'],12);self.assertEqual(row['n_infusions'],21);self.assertFalse(row['complete_final_sequence']);self.assertIsNone(row['clinical_risk_label'])
    def test_full_length_pn_corrected_affinities_kept_separate(self):
        rows=json.loads((D/'full_length_TCR_CD3_affinities.json').read_text());v={x['evidence_id']:x for x in rows}
        self.assertEqual(v['PRIMARY:37100770:PN45545:MAGEA4']['kd_um'],6.1);self.assertEqual(v['PRIMARY:37100770:PN45545:MAGEA8']['kd_um'],420)
        self.assertIsNone(v['PRIMARY:37100770:PN45545:MAGEA8']['structure_complex_pdb']);self.assertTrue(all(not x['training_allowed'] for x in rows))
    def test_inspection_release_cannot_run_optimizer(self):
        release=ROOT/'data/external_evidence/20261018-primary-reviewed-affinity-v3'
        if not release.exists():self.skipTest('local raw snapshot/release not in Git')
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'model'
            with self.assertRaisesRegex(ValueError,'insufficient independent'):run(release,output)
            self.assertFalse(output.exists())
if __name__=='__main__':unittest.main()
