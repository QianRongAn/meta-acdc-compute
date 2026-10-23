"""Mandatory source, cell, sequence and exclusion checks; no optional-data skips."""
import hashlib,json,tempfile,unittest
from pathlib import Path
from scripts.curate_admission_followup import patent_tables,st25
from scripts.audit_followup_identities import audit
from test_six_blockers import cells
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'docs/admission_followup_20261022'
def load(name):return json.loads((D/'derived'/name).read_text())

class AdmissionFollowupTests(unittest.TestCase):
    def test_all_required_primary_sources_and_publisher_versions(self):
        document=load('source_checks.json')
        for p,h in document['source_sha256'].items():
            self.assertEqual(hashlib.sha256((ROOT/p).read_bytes()).hexdigest(),h,p)
        for pair in document['publisher_Zenodo_pairs']:
            self.assertEqual((ROOT/pair['publisher_file']).read_bytes(),(ROOT/pair['zenodo_file']).read_bytes())
        self.assertEqual(len(document['publisher_Zenodo_pairs']),2)

    def test_sequences_match_published_entries_and_st25_declared_lengths(self):
        d=load('sequence_evidence.json');self.assertEqual(len(d['receptors']),6)
        for row in d['receptors']:
            self.assertIn(row['cdr3_alpha'],row['alpha_published'])
            self.assertIn(row['cdr3_beta'],row['beta_published'])
            if row['receptor_name'].startswith('A3-'):
                for chain in ['alpha','beta']:
                    sid=row[chain+'_seq_id']
                    self.assertEqual(st25(d['selected_ST25_text'][str(sid)],sid),row[chain+'_published'])
            else:
                ids=[row[k] for k in ['alpha_seq_id','beta_seq_id','paired_seq_id']]
                table=patent_tables(D/'sources'/(row['source']+'.html'),ids)
                self.assertEqual(table[row['alpha_seq_id']],row['alpha_published'])
                self.assertEqual(table[row['beta_seq_id']],row['beta_published'])
            self.assertFalse(row['training_allowed']);self.assertFalse(row['independent_acceptance_allowed'])

    def test_parser_rejects_conflicting_ids_and_incomplete_st25(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'table.html'
            path.write_text('<table>(SEQ ID NO: 1)ACD</table><table>(SEQ ID NO: 1)ACE</table>')
            with self.assertRaises(ValueError):patent_tables(path,[1])
        with self.assertRaises(ValueError):st25('<210> 9\n<211> 3\n<400> 9\nAla Gly\n',9)

    def test_all_760_raw_cells_controls_and_derived_values(self):
        d=load('kras_raw_cell_ledger.json');source=cells(ROOT/d['source_file'],d['sheet'])
        wt='VVVGAVGVGK'
        for row in d['observations']:
            self.assertEqual(source[row['source_cell']],row['raw_NFAT_EGFP_value'])
            low=source[row['minimum_cell']];high=source[row['maximum_cell']]
            self.assertEqual(low,row['minimum_Jurkat_background']);self.assertEqual(high,row['maximum_PMA_ionomycin'])
            self.assertAlmostEqual(row['derived_PMA_percent_unclipped'],100*(row['raw_NFAT_EGFP_value']-low)/(high-low),places=10)
            if row['source_cell'].startswith('P'):
                label=source['O'+''.join(c for c in row['source_cell'] if c.isdigit())]
                position=int(label[:-1])-7;aa=label[-1]
            else:
                position=ord(row['source_cell'][0])-ord('C')
                aa=source['B'+''.join(c for c in row['source_cell'] if c.isdigit())]
            self.assertEqual(row['peptide'],wt[:position]+aa+wt[position+1:])
            self.assertFalse(row['training_allowed'])
        self.assertEqual(len(d['observations']),760)
        self.assertEqual(len({(r['receptor_name'],r['peptide']) for r in d['observations']}),744)
        self.assertEqual(sum(r['is_parent_sequence'] for r in d['observations']),20)
        self.assertFalse(d['binding_labels_created'])

    def test_identity_overlay_recomputed_and_reserved_rows_excluded(self):
        d=load('identity_overlay.json');self.assertEqual(d,audit(ROOT))
        self.assertEqual(d['summary']['provisional_joint_blocks'],7)
        for x in d['exposure']:
            if x['receptor_name'] in ['A3-05','A3-10']:
                self.assertTrue(x['hits']['docs/database_binding/combined_reservation_registry.json'])
            self.assertFalse(x['independent_acceptance_allowed'])
        self.assertIsNone(d['partition_assignment']);self.assertFalse(d['test_labels_used']);self.assertFalse(d['model_scored'])

if __name__=='__main__':unittest.main()
