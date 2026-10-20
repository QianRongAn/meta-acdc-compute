import copy,json,tempfile,unittest
from pathlib import Path
from meta_acdc.data.evidence_completion import CURATION
from meta_acdc.data.evidence_sequence import validate_1d3_reference,validate_pending_pap
from meta_acdc.data.spr_followup import paired_receptor_fasta
from meta_acdc.data.experimental_evidence import assert_reserved_absent
ROOT=Path(__file__).resolve().parents[1]
class EvidenceSequence(unittest.TestCase):
    def load(self,n):return json.loads((CURATION/n).read_text())
    def ref(self):return self.load('1D3_published_reference_DNA.json')
    def pap(self):return self.load('PAP_SPR_pending.json')
    def test_published_DNA_reference_is_complete_but_not_clinical_lot(self):
        x=self.ref();self.assertTrue(validate_1d3_reference(x,check_sources=False));self.assertEqual((len(x['vector_DNA']),len(x['full_insert_DNA']),len(x['reference_protein'])),(6338,1797,598));self.assertEqual(x['parts']['alpha']['protein'][-1],'*');self.assertFalse(x['clinical_exact_insert_verified'])
    def test_reference_cannot_promote_to_lot_or_training(self):
        for key in ['clinical_exact_insert_verified','clinical_lot_verified','clinical_vector_sequence_verified','training_allowed','model_scored']:
            x=self.ref();x[key]=True
            with self.assertRaisesRegex(ValueError,'clinical proof'):validate_1d3_reference(x,check_sources=False)
    def test_source_coordinates_or_chain_order_cannot_be_changed(self):
        x=self.ref();x['parts']['beta']['start']+=1
        with self.assertRaisesRegex(ValueError,'coordinates'):validate_1d3_reference(x,check_sources=False)
        x=self.ref();x['chain_order']='alpha-P2A-beta'
        with self.assertRaisesRegex(ValueError,'order'):validate_1d3_reference(x,check_sources=False)
    def test_ambiguous_DNA_or_source_rows_cannot_be_repaired_silently(self):
        x=self.ref();x['vector_DNA']='N'+x['vector_DNA'][1:]
        with self.assertRaisesRegex(ValueError,'transcription'):validate_1d3_reference(x,check_sources=False)
        x=self.ref();x['rows'][0]['text']='A'+x['rows'][0]['text'][1:]
        with self.assertRaisesRegex(ValueError,'transcription'):validate_1d3_reference(x,check_sources=False)
    def test_explicit_engineered_cysteine_positions(self):
        x=self.ref();self.assertEqual(x['parts']['alpha']['protein'][179],'C');self.assertEqual(x['parts']['beta']['protein'][190],'C');x['alpha_cysteine_position']=181
        with self.assertRaisesRegex(ValueError,'cysteines'):validate_1d3_reference(x,check_sources=False)
    def test_numbered_TCR_header_accepted_but_MHC_not_selected(self):
        s='>9NMU_1|Chain A|HLA alpha chain|Homo sapiens\nAAAA\n>9NMU_4|Chain D|TCR156 alpha chain|Homo sapiens\nCAVNNARLMF\n>9NMU_5|Chain E|TCR156 beta chain|Homo sapiens\nCASSVAGSPEAFF\n'
        self.assertEqual(paired_receptor_fasta(s,'TCR156')['alpha'],'CAVNNARLMF')
        with self.assertRaisesRegex(ValueError,'missing'):paired_receptor_fasta(s,'TCR15')
    def test_same_clone_substring_not_selected(self):
        s='>1|TCR1567 alpha chain|Homo sapiens\nACDE\n>2|TCR1567 beta chain|Homo sapiens\nACDE\n'
        with self.assertRaisesRegex(ValueError,'missing'):paired_receptor_fasta(s,'TCR156')
    def test_five_PAP_measurements_are_one_pending_lineage(self):
        rows=self.pap();self.assertEqual(validate_pending_pap(rows,check_sources=False),dict(measurements=5,lineages=1,studies=1,eligible_classical_core=0));self.assertEqual(len({x['lineage'] for x in rows}),1)
    def test_kinetic_KD_must_not_replace_equilibrium_column(self):
        x=self.pap();x[0]['kd_um']=x[0]['kinetic_derived_kd_um']
        with self.assertRaisesRegex(ValueError,'column'):validate_pending_pap(x,check_sources=False)
    def test_format_and_unreported_conditions_remain_explicit(self):
        for key,value in [('classical_core_eligible',True),('buffer','PBS'),('immobilized_ligand_reported','monomer')]:
            rows=self.pap();rows[0][key]=value
            with self.assertRaises(ValueError):validate_pending_pap(rows,check_sources=False)
    def test_partial_candidates_and_old_core_are_unchanged(self):
        old=ROOT/'docs/diversity_followup_20261018'
        for name in ['external_differentV_candidates.json','external_classical_panel.json','review_ledger.json']:
            self.assertEqual(self.load(name),json.loads((old/name).read_text()))
        x=self.load('TW10_workbook_audit.json');self.assertFalse(x['has_full_variable_sequences']);self.assertIn('ambiguous',x['matched_rows'][1]['decision']);self.assertEqual(x['pdb_name_collision']['pdb'],'8V4Z')
    def test_new_reservations_cannot_be_bypassed_by_empty_registry(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'empty.json';p.write_text('{"receptors":[]}')
            rows=self.pap();x=self.ref();rows+=[dict(cdr3_alpha=x['cdr3_alpha'],cdr3_beta=x['cdr3_beta'],alpha_variable=x['parts']['alpha']['protein'][:-1],beta_variable=x['parts']['beta']['protein'])]
            for x in rows:
                with self.assertRaisesRegex(ValueError,'reserved'):assert_reserved_absent([x],p)
    def test_PAP_deposited_position_numbering_is_not_assumed(self):
        rows=self.pap();wt=rows[0]['alpha_variable'];m=rows[2]['alpha_variable'];self.assertEqual([(i+1,a,b) for i,(a,b) in enumerate(zip(wt,m)) if a!=b],[(33,'S','M')]);self.assertTrue(all(rows[0]['beta_variable']==x['beta_variable'] for x in rows))
    def test_local_primary_sources_and_deposited_variants(self):
        if not (ROOT/self.ref()['sources'][0]['path']).exists():self.skipTest('raw thirdparty primary sources localonly')
        validate_1d3_reference(self.ref());validate_pending_pap(self.pap())
if __name__=='__main__':unittest.main()
