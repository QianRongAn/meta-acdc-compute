import copy, hashlib, json, unittest
from pathlib import Path
from meta_acdc.data.evidence_completion import CURATION
from meta_acdc.data.spr_followup import paired_receptor_fasta, validate_classical_additions
from meta_acdc.data.experimental_evidence import assert_reserved_absent, reserved_hits
ROOT=Path(__file__).resolve().parents[1]
class SPRFollowup(unittest.TestCase):
    def load(self,n):return json.loads((CURATION/n).read_text())
    def rows(self):return self.load('external_classical_additions.json')
    def test_three_measurements_not_three_receptors(self):
        self.assertEqual(validate_classical_additions(self.rows(),check_sources=False),dict(measurements=3,lineages=2,studies=2,research_groups=1))
    def test_paired_selector_does_not_take_MHC_alpha(self):
        s='>one|Clone T cell receptor alpha chain\nCAVAAAAAAF\n>two|Clone T cell receptor beta chain\nCASSGGGGF\n>three|MHC alpha chain\nAAAAAAAAAA\n'
        self.assertEqual(paired_receptor_fasta(s,'Clone')['alpha'],'CAVAAAAAAF')
    def test_missing_named_beta_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'paired'):paired_receptor_fasta('>x|MHC alpha chain\nAAAA','Clone')
    def test_ambiguous_named_chain_rejected(self):
        with self.assertRaisesRegex(ValueError,'ambiguous'):paired_receptor_fasta('>1|Clone T cell receptor alpha chain\nAAAA\n>2|Clone T cell receptor alpha chain\nAAAA','Clone')
    def test_units_error_rejected(self):
        x=self.rows();x[0]['kd_nm']/=1000
        with self.assertRaisesRegex(ValueError,'conversion'):validate_classical_additions(x,check_sources=False)
    def test_unsupported_endpoint_promotions(self):
        for key,value in [('model_scored',True),('training_allowed',True),('clinical_risk_label',0),('binary_binding_label',0)]:
            x=self.rows();x[0][key]=value
            with self.assertRaisesRegex(ValueError,'promotion'):validate_classical_additions(x,check_sources=False)
    def test_noncanonical_norleucine_never_silently_coded(self):
        x=self.rows();x[0]['peptide']='I[Nle]DQVPFSV'
        with self.assertRaisesRegex(ValueError,'noncanonical'):validate_classical_additions(x,check_sources=False)
    def test_missing_temperature_fails(self):
        x=self.rows();x[0]['temperature_K']=None
        with self.assertRaisesRegex(ValueError,'temperature'):validate_classical_additions(x,check_sources=False)
    def test_chain_CDR_mismatch_fails(self):
        x=self.rows();x[0]['cdr3_beta']='CASSWRONGSEQF'
        with self.assertRaisesRegex(ValueError,'identity'):validate_classical_additions(x,check_sources=False)
    def test_source_specific_T4H2_values_not_averaged(self):
        x=[a for a in self.rows() if a['lineage']=='T4H2'];self.assertEqual({a['kd_nm'] for a in x},{1800.,6800.});self.assertEqual(len({a['study'] for a in x}),2)
    def test_new_lineages_cannot_enter_training(self):
        for x in self.rows()[:2]:
            with self.assertRaisesRegex(ValueError,'reserved'):assert_reserved_absent([x])
    def test_prior_panel_and_folds_unchanged(self):
        panel=self.load('external_classical_panel.json');old=json.loads((ROOT/'docs/evidence_completion_20261018/external_classical_core.json').read_text());self.assertEqual(panel[:4],old);self.assertEqual(len(panel),7)
        rows=self.load('review_ledger.json');self.assertEqual(sum(a['training_allowed'] for a in rows),25)
        for a in rows:self.assertEqual(a['parent_fold'],a['original_record']['fold'])
    def test_recovered_fulltext_is_not_exact_numeric_confirmation(self):
        x={a['evidence_id']:a for a in self.load('review_ledger.json')}
        for key in ['ATLAS:149','ATLAS:154']:
            self.assertTrue(x[key]['primary_text_retrieved']);self.assertEqual(x[key]['temperature_K'],298.15);self.assertEqual(x[key]['status'],'quarantine_primary_value_unverified');self.assertFalse(x[key]['training_allowed'])
    def test_historical_P1_15_not_counted_fresh(self):
        x=self.load('28IL_history_check.json');self.assertTrue(any(h['hits'] for h in x['checks']['history']));self.assertTrue(all(a['receptor_name']!='P1-15' for a in self.rows()))
    def test_primary_sources_and_named_chains_when_local(self):
        if not (ROOT/self.rows()[0]['sources'][0]['path']).exists():self.skipTest('third-party raw primary source is local-only')
        validate_classical_additions(self.rows())
    def test_previous_frozen_artifacts_unchanged(self):
        for n,digest in self.load('previous_freezes.json').items():self.assertEqual(hashlib.sha256((ROOT/n).read_bytes()).hexdigest(),digest,n)
if __name__=='__main__':unittest.main()
