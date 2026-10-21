"""Mandatory primary-cell, duplicate-scale, identity-block and no-promotion checks."""
import hashlib,json,math,unittest
from pathlib import Path
from test_six_blockers import cells
from scripts.audit_measured_intake_connections import components,author_key,audit
from scripts.curate_measured_endpoint_intake import require_unique
R=Path(__file__).resolve().parents[1];D=R/'docs/measured_intake_20261021';S=R/'docs/six_blockers_20261021/sources'
def load(n):return json.loads((D/n).read_text())
def col(n):
    s=''
    while n:n,k=divmod(n-1,26);s=chr(65+k)+s
    return s
class MeasuredIntakeTests(unittest.TestCase):
    def test_required_sources_and_identity_only_inputs_unchanged(self):
        for document in [load('primary_observations.json'),load('connections.json')]:
            for p,h in document['source_sha256'].items():self.assertEqual(hashlib.sha256((R/p).read_bytes()).hexdigest(),h,p)
        for m in [load('kras_primary_retrieval.json')]+load('kras_supplement_retrieval.json'):
            self.assertEqual(m['status'],'retrieved');self.assertEqual(hashlib.sha256((R/m['path']).read_bytes()).hexdigest(),m['sha256'])
        self.assertEqual(len(load('source_metadata_retrieval.json')),16)
        for m in load('source_metadata_retrieval.json'):
            self.assertEqual(m['status'],'retrieved');self.assertEqual(hashlib.sha256((R/m['path']).read_bytes()).hexdigest(),m['sha256'])
    def test_every_primary_response_and_replicate_reconstructed(self):
        primary=load('primary_observations.json');cached={}
        for x in primary['observations']:
            key=(x['source_file'],x['source_sheet'])
            if key not in cached:cached[key]=cells(S/key[0],key[1])
            self.assertAlmostEqual(x['value'],cached[key][col(x['source_excel_column'])+str(x['source_excel_row'])],places=10)
            expected=sum(x['normalized_replicates'])/2 if x['normalized_replicates'] else x['positive_control']['unnormalized_activation']*x['positive_control']['normalization_factor']
            self.assertAlmostEqual(x['value'],expected,places=10)
            if x['positive_control']:self.assertAlmostEqual(x['positive_control']['mean_PMA_ionomycin']*x['positive_control']['normalization_factor'],100,places=10)
        self.assertEqual(len({x['evidence_id'] for x in primary['observations']}),4378)
    def test_secondary_is_duplicate_and_retains_different_scales(self):
        p={x['evidence_id']:x for x in load('primary_observations.json')['observations']};s=load('secondary_reconciliation.json');raw=cells(S/'batcave_mhci.xlsx','Sheet1')
        self.assertEqual(len(s['rows']),4378)
        for x in s['rows']:
            obs=p[x['primary_evidence_id']];row=str(x['secondary_excel_row'])
            self.assertEqual(raw['Q'+row],obs['peptide']);self.assertAlmostEqual(raw['R'+row],x['secondary_value'],places=10)
            self.assertFalse(x['admit_as_new_measurement']);self.assertTrue(x['assay_name_conflict'])
        self.assertEqual(sum(x['matches_primary_normalized'] for x in s['rows']),3440)
        self.assertEqual(sum(x['matches_primary_unnormalized'] for x in s['rows']),938)
        self.assertEqual(s['new_independent_measurements'],0)
    def test_wild_type_denominator_and_sequence_quarantine(self):
        p=load('primary_observations.json');valid=[x for x in p['observations'] if x['sequence_review_pass']]
        self.assertEqual(len(valid),4072);self.assertEqual(sum(x['is_cognate_parent'] for x in valid),25)
        self.assertEqual(sum(x['source_file']=='pteam_rnf43.xlsx' and not x['is_cognate_parent'] for x in valid),798)
        bad={x['receptor_name'] for x in p['receptors'] if x['sequence_reasons']}
        self.assertEqual(bad,{'TCR1-4','R27'});self.assertEqual(len(p['authors']),11)
        self.assertEqual({x['normalization'] for x in valid},{'TCR_transgenic_fraction','PMA_ionomycin_positive_control'})
    def test_transitive_source_and_family_blocking_not_relabelled_families(self):
        self.assertEqual(components(['A','B','C','D'],[dict(a='A',b='B'),dict(a='B',b='C')]),[['A','B','C'],['D']])
        self.assertEqual(author_key('Drost','Felix'),author_key('Drost','F'))
        self.assertIsNone(author_key('Drost',''))
        with self.assertRaises(ValueError):require_unique([dict(id=1),dict(id=1)],'id')
        a=load('connections.json');self.assertEqual(a['summary']['primary_joint_blocks'],1)
        block=next(g for g in a['block_review'] if 'PTEAM:cmv:TCR2-4' in g['members'])
        self.assertTrue(block['prior_family_exposure']);self.assertFalse(block['independent_acceptance_allowed'])
        self.assertEqual(a['summary']['primary_receptors_with_prior_reservation_hits'],7)
        self.assertEqual(a['summary']['unresolved_secondary_author_pmids'],[])
    def test_entire_graph_reproduces_and_no_labels_or_model_promotion(self):
        expected=load('connections.json');self.assertEqual(audit(R),expected)
        self.assertIsNone(expected['partition_assignment']);self.assertFalse(expected['frozen_test_labels_used_for_tuning'])
        self.assertFalse(expected['summary']['new_frozen_test_created'])
        for x in load('primary_observations.json')['observations']:
            for k in ['molecular_binding_label','clinical_label','binary_activation_label','peptide_pulse_concentration']:self.assertIsNone(x[k])
            self.assertFalse(x['training_allowed']);self.assertFalse(x['model_scored'])
if __name__=='__main__':unittest.main()
