import copy,hashlib,json,tempfile,unittest,zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from meta_acdc.data.evidence_completion import CURATION,validate_c4,validate_external,structure_registry
from meta_acdc.data.spr_primary import validate_review
from meta_acdc.data.experimental_evidence import reserved_hits,assert_reserved_absent
ROOT=Path(__file__).resolve().parents[1]
class EvidenceCompletion(unittest.TestCase):
    def load(self,name):return json.loads((CURATION/name).read_text())
    def test_primary_51_covered_and_31_eligible(self):
        rows=self.load('review_ledger.json');parent=json.loads((ROOT/'docs/clinical_primary_20261018/affinity_dataset.json').read_text())
        self.assertTrue(validate_review(rows,parent,check_sources=False));self.assertEqual(sum(x['training_allowed'] for x in rows),31)
    def test_tk3_q5_value_is_not_wt(self):
        row=next(x for x in self.load('review_ledger.json') if x['evidence_id']=='ATLAS:914')
        self.assertEqual(row['peptide'],'HPVGQADYFEY');self.assertEqual(row['status'],'quarantine_wrong_ligand_value');self.assertFalse(row['training_allowed'])
    def test_two_original_temperatures_stay_missing(self):
        rows={x['evidence_id']:x for x in self.load('review_ledger.json')}
        for key in ['714','716']:self.assertIsNone(rows['ATLAS:'+key]['temperature_K']);self.assertFalse(rows['ATLAS:'+key]['training_allowed'])
    def test_unverified_database_error_bars_removed(self):
        rows={x['evidence_id']:x for x in self.load('review_ledger.json')}
        for key in ['502','505','485']:self.assertIsNone(rows['ATLAS:'+key]['primary_error'])
    def test_c25_c7_structure_mhc_mismatch_explicit(self):
        for row in self.load('review_ledger.json'):
            if row['evidence_id'] in ['ATLAS:644','ATLAS:645']:self.assertFalse(row['structure_assay_mhc_match']);self.assertFalse(row['structure_based_affinity_training_allowed'])
    def test_c4_published_reference_translation(self):self.assertTrue(validate_c4(self.load('C4_published_reference_cassette.json')))
    def test_c4_wrong_nucleotide_rejected(self):
        x=self.load('C4_published_reference_cassette.json');x['dna']='TAA'+x['dna'][3:]
        with self.assertRaisesRegex(ValueError,'translation'):validate_c4(x)
    def test_c4_lot_claim_rejected(self):
        x=self.load('C4_published_reference_cassette.json');x['clinical_lot_sequence_certificate']=True
        with self.assertRaisesRegex(ValueError,'clinical'):validate_c4(x)
    def test_c4_wrong_order_rejected(self):
        x=self.load('C4_published_reference_cassette.json');x['chain_order']='alpha-P2A-beta'
        with self.assertRaisesRegex(ValueError,'order'):validate_c4(x)
    def test_six_external_pairs_are_one_lineage(self):
        x=validate_external(self.load('external_HLAE_challenge.json'));self.assertEqual(x['measurements'],6);self.assertEqual(x['biological_receptors'],1);self.assertEqual(x['studies'],1)
    def test_external_molar_conversion_rejects_millionfold_error(self):
        x=self.load('external_HLAE_challenge.json');x[0]['kd_um']/=1e6
        with self.assertRaisesRegex(ValueError,'conversion'):validate_external(x)
    def test_external_not_training_or_binary_negative(self):
        for key,val in [('training_allowed',True),('binary_binding_label',0),('model_scored',True)]:
            x=self.load('external_HLAE_challenge.json');x[0][key]=val
            with self.assertRaisesRegex(ValueError,'promotion'):validate_external(x)
    def test_historical_alpha_beta_fields_actually_compared(self):
        x=self.load('external_HLAE_challenge.json')[0];ref=structure_registry([dict(complex_id='historical',alpha=x['alpha_variable'],beta=x['beta_variable'])]);self.assertEqual(reserved_hits(x,ref),['historical'])
    def test_missing_historical_identity_fails_closed(self):
        with self.assertRaisesRegex(ValueError,'identity'):structure_registry([{'complex_id':'broken'}])
    def test_external_core_is_previous_reservation_not_reshuffle(self):
        core=self.load('external_classical_core.json');old={x['evidence_id']:x for x in json.loads((ROOT/'docs/evidence/primary_measurements.json').read_text())}
        self.assertEqual(len(core),4)
        for x in core:
            self.assertTrue(x['not_newly_discovered']);self.assertFalse(x['model_scored']);self.assertEqual(x['kd_nm'],old[x['evidence_id']]['kd_nm']);self.assertEqual(x['alpha_variable'],old[x['evidence_id']]['alpha_variable'])
    def test_new_HBV_reservation_cannot_be_bypassed(self):
        with self.assertRaisesRegex(ValueError,'reserved'):assert_reserved_absent([self.load('external_HLAE_challenge.json')[0]])
    def test_new_C4_reservation_protects_exact_pair(self):
        c=self.load('C4_published_reference_cassette.json')
        with self.assertRaisesRegex(ValueError,'reserved'):assert_reserved_absent([{'alpha':c['alpha_chain'],'beta':c['beta_chain'],'cdr3_alpha':c['cdr3_alpha'],'cdr3_beta':c['cdr3_beta']}])
    def test_earlier_clinical_and_spr_freezes_unchanged(self):
        for n,v in self.load('previous_freezes.json').items():self.assertEqual(hashlib.sha256((ROOT/n).read_bytes()).hexdigest(),v,n)
    def test_external_raw_source_cells_match_recorded_replicates(self):
        p=ROOT/'data/external_evidence/20261018-evidence-completion-v1/HBV_7.xlsx'
        if not p.exists():self.skipTest('third-party raw workbook is local-only')
        ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        with zipfile.ZipFile(p) as z:
            workbook=ET.fromstring(z.read('xl/workbook.xml'));names=[a.get('name') for a in workbook.findall('s:sheets/s:sheet',ns)];index=names.index('Figure 2c')+1
            sheet=ET.fromstring(z.read('xl/worksheets/sheet'+str(index)+'.xml'));cells={c.get('r'):float(c.find('s:v',ns).text) for c in sheet.findall('.//s:c',ns) if c.get('t') not in ['s','inlineStr'] and c.find('s:v',ns) is not None}
        for i,x in enumerate(self.load('external_HLAE_challenge.json'),12):self.assertEqual(x['raw_kd_replicates'],[cells[col+str(i)] for col in 'DEF'])
if __name__=='__main__':unittest.main()
