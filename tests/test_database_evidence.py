"""Experimental endpoint, identity, split and clinical-label safeguards."""
import copy,json,tempfile,unittest
from pathlib import Path
from meta_acdc.data.atlas_binding import kd_measurement,fasta_records,paired_from_sources,grouped_folds,verify_release
from meta_acdc.data.functional_benchmark import validate_rows
from meta_acdc.data.clinical_database import normalize_trial,molecular_label
from meta_acdc.data.release import sha256

ROOT=Path(__file__).resolve().parents[1]

class Measurement(unittest.TestCase):
    def test_exact_microM_to_pkd(self):self.assertAlmostEqual(kd_measurement('10')['pkd'],5)
    def test_reported_uncertainty_is_preserved(self):self.assertEqual(kd_measurement('12.5 +/- 2.1')['error_um'],2.1)
    def test_censored_bounds_are_not_exact(self):
        x=kd_measurement('>200');self.assertIsNone(x['pkd']);self.assertEqual(x['kd_um'],200);self.assertEqual(x['qualifier'],'>')
    def test_not_determined_is_not_a_nonbinder(self):self.assertEqual(kd_measurement('n.d.')['status'],'missing_or_unparsed')
    def test_nan_infinite_zero_do_not_train(self):
        for x in ['nan','inf','0','-1']:self.assertIsNone(kd_measurement(x)['pkd'])
    def test_approximate_is_not_exact(self):self.assertIsNone(kd_measurement('~5')['pkd'])
    def test_fasta_author_chain_is_used(self):self.assertEqual(fasta_records('>1_1|Chains A[auth D], B[auth E]|TCR|Homo sapiens (9606)\nACDEF\n')[0]['chains'],['D','E'])
    def test_missing_paired_chain_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'Alpha'):paired_from_sources(fasta_records('>1_1|Chain D|TCR|Mus musculus (10090)\n'+'A'*90),[])

class Splits(unittest.TestCase):
    def rows(self):
        data=json.loads((ROOT/'docs/database_binding/affinity_dataset.json').read_text());return copy.deepcopy(data)
    def test_connected_blocks_never_cross_folds(self):
        rows=self.rows();folds,info=grouped_folds(rows);seen={}
        for row in rows:
            for key in [('study',row['pmid']),('pair',row['pair_id']),('block',row['split_block'])]:
                if key in seen:self.assertEqual(seen[key],row['fold'])
                seen[key]=row['fold']
        self.assertEqual(set(sum(folds.values(),[])),set(range(len(rows))))
    def test_exact_peptide_never_crosses_folds(self):
        rows=self.rows();grouped_folds(rows);peptides={}
        for row in rows:
            if row['peptide'] in peptides:self.assertEqual(peptides[row['peptide']],row['fold'])
            peptides[row['peptide']]=row['fold']
    def test_too_few_independent_blocks_fails(self):
        rows=self.rows()[:2]
        with self.assertRaisesRegex(ValueError,'three'):grouped_folds(rows)

class Functional(unittest.TestCase):
    def group(self):
        aa='A'*80+'CASSFGAG';peps=['ACDEFGHIK','CDEFGHIKL','DEFGHIKLM','EFGHIKLMN','FGHIKLMNP','GHIKLMNPQ','HIKLMNPQR','IKLMNPQRS','KLMNPQRST','LMNPQRSTV']
        return [dict(tcra_trimmed=aa,tcrb_trimmed=aa,tcra_cdr3='CASSF',tcrb_cdr3='CASSF',peptide=p,hla='A*02:01',label='1' if i==0 else '0') for i,p in enumerate(peps)]
    def test_activation_label_does_not_become_binding_label(self):self.assertTrue(all(x['binding_label'] is None and x['clinical_label'] is None for x in validate_rows(self.group())))
    def test_incomplete_negative_matrix_rejected(self):
        with self.assertRaisesRegex(ValueError,'matrix'):validate_rows(self.group()[:-1])
    def test_duplicate_pair_rejected(self):
        rows=self.group();rows.append(rows[0])
        with self.assertRaisesRegex(ValueError,'duplicate'):validate_rows(rows)
    def test_two_positive_per_tcr_rejected(self):
        rows=self.group();rows[1]['label']='1'
        with self.assertRaisesRegex(ValueError,'matrix'):validate_rows(rows)
    def test_cdr3_mismatch_rejected(self):
        rows=self.group();rows[0]['tcra_cdr3']='WY'
        with self.assertRaisesRegex(ValueError,'CDR3'):validate_rows(rows)

class Clinical(unittest.TestCase):
    def trial(self):return {'protocolSection':{'identificationModule':{'nctId':'NCT00000000','briefTitle':'Trial'}},'resultsSection':{'adverseEventsModule':{'eventGroups':[{'id':'G','title':'Arm'}],'seriousEvents':[{'term':'Cardiac event','stats':[{'groupId':'G','numAffected':1,'numAtRisk':5}]}]}}}
    def test_aggregate_outcome_is_not_binding_or_causal_label(self):
        trial,events=normalize_trial(self.trial(),{});self.assertIsNone(events[0]['binding_label']);self.assertIsNone(events[0]['causal_offtarget_label']);self.assertIsNone(events[0]['peptide'])
        with self.assertRaisesRegex(ValueError,'Aggregate'):molecular_label(events[0])
    def test_missing_count_stays_unknown_not_zero(self):
        trial=self.trial();trial['resultsSection']['adverseEventsModule']['seriousEvents'][0]['stats'][0].pop('numAffected');self.assertIsNone(normalize_trial(trial,{})[1][0]['num_affected'])
    def test_affected_cannot_exceed_denominator(self):
        trial=self.trial();trial['resultsSection']['adverseEventsModule']['seriousEvents'][0]['stats'][0]['numAffected']=6
        with self.assertRaisesRegex(ValueError,'invalid'):normalize_trial(trial,{})
    def test_unknown_group_rejected(self):
        trial=self.trial();trial['resultsSection']['adverseEventsModule']['seriousEvents'][0]['stats'][0]['groupId']='unknown'
        with self.assertRaisesRegex(ValueError,'unknown'):normalize_trial(trial,{})

class Provenance(unittest.TestCase):
    def test_tampered_affinity_artifact_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);(root/'raw').mkdir();(root/'data.json').write_text('original');(root/'registry.json').write_text('{}')
            m={'clinical_registry_sha256':sha256(ROOT/'docs/clinical_primary_20261018/clinical_reservations.json'),'task':'experimental_SPR_pKD_regression','code_files':{'src/meta_acdc/data/atlas_binding.py':sha256(ROOT/'src/meta_acdc/data/atlas_binding.py')},'source_snapshot':str(root/'raw'),'files':{'data.json':sha256(root/'data.json')},'sources':[],'registry_path':str(root/'registry.json'),'external_registry_sha256':sha256(root/'registry.json')};(root/'manifest.json').write_text(json.dumps(m));verify_release(root);(root/'data.json').write_text('changed')
            with self.assertRaisesRegex(ValueError,'artifact'):verify_release(root)

if __name__=='__main__':unittest.main()
