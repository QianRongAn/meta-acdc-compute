"""Regression tests for assay identity, conditional labels and external holdouts."""
import copy,json,tempfile,unittest
from pathlib import Path
from meta_acdc.data.experimental_evidence import (assay_ids,endpoint,finite_nm,make_ledger,condition_conflicts,verify_pilot,assert_reserved_absent)

ROOT=Path(__file__).resolve().parents[1]


def fixture(aid=11,outcome='Negative'):
    r={'receptor__type':'alphabeta','assay__type':'T cell','assay__iedb_ids':str(aid),'receptor__iedb_receptor_id':23,'receptor_group_id':7,'reference__iedb_iri':'https://www.iedb.org/reference/1','epitope__iedb_iri':'https://www.iedb.org/epitope/8','epitope__name':'SIINFEKL'}
    for i,kind in [(1,'alpha'),(2,'beta')]:r.update({f'chain_{i}__type':kind,f'chain_{i}__organism_iri':'http://purl.obolibrary.org/obo/NCBITaxon_9606',f'chain_{i}__protein_sequence':'A'*80+'CASSFGAG',f'chain_{i}__cdr3_calculated':'ASS'})
    a={'assay_id':aid,'epitope__name':'SIINFEKL','mhc_restriction__name':'HLA-A*02:01','assay__method':'surface plasmon resonance (SPR)','assay__response_measured':'dissociation constant KD','assay__qualitative_measurement':outcome,'reference__pmid':'1'}
    s={'tcell_id':aid,'linear_sequence':'SIINFEKL','structure_id':8,'reference_id':1,'receptor_ids':[23],'tcr_receptor_group_ids':[7],'mhc_allele_name':'HLA-A*02:01'}
    return r,a,s


class ExperimentalEvidence(unittest.TestCase):
    def test_database_negative_does_not_become_universal_training_label(self):
        r,a,s=fixture();x=make_ledger([r],[a],[s],[])[0]
        self.assertEqual(x['observed_outcome'],'conditional_negative');self.assertEqual(x['label'],'');self.assertEqual(x['use'],'review_only')
    def test_missing_alpha_does_not_inherit_group_sequence(self):
        r,a,s=fixture();r.pop('chain_1__protein_sequence');x=make_ledger([r],[a],[s],[])[0]
        self.assertEqual(x['alpha_variable'],'');self.assertIn('complete canonical paired',x['screening_reasons'])
    def test_group_match_cannot_substitute_for_exact_receptor_assay_link(self):
        r,a,s=fixture();s['receptor_ids']=[24];x=make_ledger([r],[a],[s],[])[0]
        self.assertIn('assay_receptor_link_mismatch',x['screening_reasons'])
    def test_epitope_mismatch_is_quarantined(self):
        r,a,s=fixture();s['linear_sequence']='GILGFVFTL';self.assertIn('mismatched_peptide',make_ledger([r],[a],[s],[])[0]['screening_reasons'])
    def test_condition_conflicts_are_retained(self):
        r,a,s=fixture();r2,a2,s2=fixture(12,'Positive');a['effector_cell__culture_condition']='condition A';a2['effector_cell__culture_condition']='condition B'
        rows=make_ledger([r,r2],[a,a2],[s,s2],[]);self.assertEqual(len(rows),2);self.assertEqual(len(condition_conflicts(rows)),1)
        self.assertEqual({x['effector_culture'] for x in rows},{'condition A','condition B'});self.assertTrue(all(x['label']=='' for x in rows))
    def test_censored_kd_retains_inequality(self):
        r,a,s=fixture();a.update(assay__quantitative_measurement=200000,assay__measurement_inequality='>',assay__units='nM');x=make_ledger([r],[a],[s],[])[0]
        self.assertEqual(x['qualifier'],'>');self.assertEqual(x['value_nm'],200000);self.assertIn('censored_measurement',x['screening_reasons']);self.assertEqual(x['label'],'')
    def test_endpoints_are_separate(self):
        self.assertNotEqual(endpoint('ELISA','IFNg release'),endpoint('surface plasmon resonance (SPR)','dissociation constant KD'))
        self.assertEqual(endpoint('x-ray crystallography','3D structure'),'other_assay')
    def test_invalid_assay_id_and_unit_are_refused(self):
        with self.assertRaises(ValueError):assay_ids('12.3')
        self.assertEqual(assay_ids('12; 13,12'),[12,13]);self.assertIsNone(finite_nm(2,'nmol'));self.assertIsNone(finite_nm(float('inf'),'nM'))
    def test_ambiguous_receptor_exports_are_preserved_and_quarantined(self):
        r,a,s=fixture();r2=copy.deepcopy(r);r2['chain_1__protein_sequence']='G'+r2['chain_1__protein_sequence'][1:]
        rows=make_ledger([r,r2],[a],[s],[]);self.assertEqual(len(rows),2);self.assertTrue(all('ambiguous_receptor_export' in x['screening_reasons'] for x in rows))
    def test_duplicate_assay_source_ids_fail_closed(self):
        r,a,s=fixture()
        with self.assertRaises(ValueError):make_ledger([r],[a,a],[s],[])


class ExternalPanel(unittest.TestCase):
    def setUp(self):self.pilot=json.loads((ROOT/'docs/evidence/p53_primary_curation.json').read_text())
    def test_panel_has_matched_positive_and_conditional_negative(self):
        result=verify_pilot(self.pilot,[]);self.assertEqual(result['n_receptors'],3);self.assertEqual(result['n_positive'],3);self.assertEqual(result['n_conditional_negative'],3);self.assertFalse(result['benchmark_allowed'])
    def test_undetected_is_not_zero_affinity(self):
        rows=copy.deepcopy(self.pilot);next(x for x in rows if x['label']==0)['kd_nm']=0
        with self.assertRaisesRegex(ValueError,'non-detection'):verify_pilot(rows,[])
    def test_known_family_is_refused_even_if_only_one_chain_matches(self):
        ref={'alpha':self.pilot[0]['alpha_variable'],'beta':'A'*100,'family_id':'reserved'}
        with self.assertRaisesRegex(ValueError,'overlap'):verify_pilot(self.pilot,[ref])
    def test_future_training_reservation_is_enforced(self):
        with self.assertRaisesRegex(ValueError,'reserved external'):assert_reserved_absent([{'alpha':self.pilot[0]['alpha_variable'],'beta':'A'*100}])
    def test_single_chain_associations_are_also_reserved(self):
        with self.assertRaisesRegex(ValueError,'reserved external'):assert_reserved_absent([{'cdr3_alpha':self.pilot[0]['cdr3_alpha'],'cdr3_beta':''}])
    def test_missing_registry_is_not_silently_ignored(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'missing'):assert_reserved_absent([],Path(d)/'absent.json')

if __name__=='__main__':unittest.main()
