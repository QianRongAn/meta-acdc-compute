"""Required checks for prospective plans and source conflicts; no outcome simulations."""
from pathlib import Path
import hashlib,json,math,unittest,xml.etree.ElementTree as ET
from scripts.validate_e5_e6_protocol import validate,stream_seed
from scripts.curate_pteam_primary import translate
from test_six_blockers import cells
R=Path(__file__).resolve().parents[1];D=R/'docs/protocol_followup_20261023'

class ProtocolConstructTests(unittest.TestCase):
    def test_protocol_parameters_and_population_error_truth(self):
        result=validate(R)
        self.assertEqual(result['E6_effect_cells_checked'],78)
        self.assertEqual(result['E5_grid_cells'],72)
        self.assertFalse(result['study_runs_executed'])
        self.assertEqual(stream_seed('E5','1.0','n14_rho0.5',0,'labels'),int(hashlib.sha256(b'E5|1.0|n14_rho0.5|0|labels').hexdigest()[:16],16))

    def test_new_required_sources_are_present_and_unchanged(self):
        d=json.loads((D/'construct_condition_audit_v2.json').read_text())
        for p,h in d['source_sha256'].items():self.assertEqual(hashlib.sha256((R/p).read_bytes()).hexdigest(),h,p)
        receipt=json.loads((D/'pteam_full_retrieval.json').read_text())
        self.assertEqual(hashlib.sha256((R/receipt['path']).read_bytes()).hexdigest(),receipt['sha256'])
        self.assertEqual(receipt['pages'],35)

    def test_cmv_conflict_reproduced_from_all_twenty_primary_DNA_cells(self):
        d=json.loads((D/'construct_condition_audit_v2.json').read_text());refs=d['constant_references']
        raw=cells(R/'docs/six_blockers_20261021/sources/pteam_cmv.xlsx','TCR sequences')
        rows=[x for x in d['pteam_constructs'] if x['source_file']=='pteam_cmv.xlsx']
        self.assertEqual(len(rows),20)
        self.assertEqual(sum(x['prior_DNA_pair_passes'] for x in rows),19)
        for x in rows:
            aa=translate(raw[x['DNA_alpha_cell']])['amino_acids']
            constant=aa[x['constant_anchor_zero_index']:]
            self.assertEqual(constant,x['observed_alpha_constant'])
            self.assertNotEqual(constant,refs['human_P01848'])
            differences=[(i+1,a,b) for i,(a,b) in enumerate(zip(refs['mouse_P01849'],constant)) if a!=b]
            self.assertEqual(len(constant),len(refs['mouse_P01849']))
            self.assertEqual(differences,[(47,'T','C')])
            self.assertTrue(x['methods_sequence_conflict']);self.assertFalse(x['training_allowed'])

    def test_both_kras_control_views_trace_to_all_760_source_cells(self):
        d=json.loads((D/'construct_condition_audit_v2.json').read_text())
        raw=cells(R/'docs/measured_intake_20261021/source_metadata/kras_measurements.xlsx','Figure 2C')
        previous=json.loads((R/'docs/admission_followup_20261022/derived/kras_raw_cell_ledger.json').read_text())
        old={(x['receptor_name'],x['source_cell']):x for x in previous['observations']}
        self.assertEqual(len(d['kras_normalization_views']),760)
        for x in d['kras_normalization_views']:
            v=raw[x['source_cell']];hi=raw[x['cognate_cell']];lo=raw[x['no_peptide_APC_cell']]
            self.assertEqual(hi,x['cognate_value']);self.assertEqual(lo,x['no_peptide_APC_value'])
            y=100*(v-lo)/(hi-lo)
            self.assertAlmostEqual(y,x['methods_interpretation_unclipped'],places=10)
            self.assertEqual(max(0,min(100,y)),x['methods_interpretation_clipped'])
            self.assertEqual(x['figure_PMA_Jurkat_interpretation_unclipped'],old[(x['receptor_name'],x['source_cell'])]['derived_PMA_percent_unclipped'])
            self.assertFalse(x['authors_exact_implemented_normalization_verified']);self.assertFalse(x['training_allowed'])

    def test_primary_methods_and_failures_not_promoted_to_new_truth(self):
        d=json.loads((D/'construct_condition_audit_v2.json').read_text())
        self.assertIn('hTRAC',d['primary_methods_locators']['sec4.4.2'])
        self.assertEqual(d['summary']['CMV_construct_conflicts'],20)
        self.assertIsNone(d['partition_assignment']);self.assertFalse(d['test_labels_used']);self.assertFalse(d['model_scored'])
        self.assertTrue(all(x['peptide_pulse_uM'] is None and not x['training_allowed'] for x in d['pteam_constructs']))
        receipts=json.loads((D/'a3_publisher_retrieval.json').read_text())
        metadata=next(x for x in receipts if x['name']=='a3_publisher_api.xml')
        self.assertFalse(metadata['usable_full_text'])
        t=ET.parse(D/'sources/a3_publisher_api.xml')
        self.assertIsNone(next((e for e in t.iter() if e.tag.split('}')[-1] in ['body','originalText']),None))

    def test_kras_peptide_identity_conflicts_cannot_become_training_labels(self):
        d=json.loads((D/'kras_raw_cell_ledger_v2.json').read_text())
        raw=cells(R/d['primary_peptide_library_file'],'Supplemental File 2')
        library={}
        wt='VVVGAVGVGK';conflicts=[]
        for row in range(4,194):
            sequence=raw['B'+str(row)].removeprefix('H-').removesuffix('-OH')
            p=int(raw['C'+str(row)])-1;aa=sequence[p]
            self.assertNotEqual(aa,wt[p]);self.assertEqual(sequence,wt[:p]+aa+wt[p+1:])
            library[sequence]=row
            if aa!=raw['E'+str(row)]:conflicts.append(row)
        self.assertEqual(len(library),190)
        self.assertEqual(len(conflicts),62)
        self.assertEqual({x['source_row'] for x in d['primary_library_annotation_conflicts']},set(conflicts))
        self.assertEqual(d['verified_identity_cells'],400)
        self.assertEqual(d['quarantined_identity_cells'],360)
        self.assertEqual(d['literal_vs_library_order_disagreements'],244)
        self.assertEqual(d['literal_parent_absent_from_library_cells'],20)
        for x in d['observations']:
            if x['peptide_identity_verified']:
                self.assertIn(x['peptide'],library)
                self.assertEqual(x['literal_library_source_row'],library[x['peptide']])
                self.assertFalse(x['library_order_mapping_is_assumed'])
            else:
                self.assertIsNone(x['peptide']);self.assertTrue(x['library_order_mapping_is_assumed'])
            self.assertIn(x['library_order_candidate_peptide'],library)
            self.assertFalse(x['training_allowed'])

if __name__=='__main__':unittest.main()
