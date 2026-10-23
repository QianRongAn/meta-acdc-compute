import unittest,json,tempfile
from pathlib import Path
from meta_acdc.structure.paired_variable_identity import *
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'docs/fingerprint_paired_identity_20261024/sources'
class PairedVariableTests(unittest.TestCase):
    def test_real_numbered_padding(self):
        domains=numbered_domains((SRC/'1ao7.html').read_text())
        self.assertEqual(domains['E']['sequence'],fasta_records((SRC/'1ao7_E.fasta').read_text())[0][1])
        self.assertEqual(domains['D']['positions'][-1],'A128')
    def test_unmatched_padding_fails(self):
        text='<td>VA</td><td>A</td><table><tr><td>A1</td><td></td></tr><tr><td>Q</td><td>E</td></tr></table>'
        with self.assertRaises(ValueError):numbered_domains(text)
    def test_exact_pair_order_and_mutation(self):
        self.assertNotEqual(pair_key('ACD','EFG'),pair_key('EFG','ACD'))
        self.assertNotEqual(pair_key('ACD','EFG'),pair_key('ACE','EFG'))
        with self.assertRaises(ValueError):pair_key('ACX','EFG')
    def test_germline_translation(self):
        refs=germline_references([SRC/'Homo_sapiens_TRAV.fasta',SRC/'Homo_sapiens_TRBJ.fasta'])
        self.assertTrue(refs[('human','TRAV12-2*01')][0]['amino_acids'].startswith('QKEVEQNSGP'))
        self.assertTrue(refs[('human','TRBJ2-7*01')][0]['amino_acids'].endswith('FGPGTRLTVT'))
    def test_missing_terminal_residue_not_imputed(self):
        audit=json.loads((ROOT/'docs/fingerprint_hard_family_20261024/identity_audit/identity_audit.json').read_text())
        row=next(x for x in audit['rows'] if x['pdb']=='1ao7');text=(SRC/'1ao7.html').read_text()
        refs=germline_references(list(SRC.glob('Homo_sapiens_TR*.fasta')))
        with self.assertRaisesRegex(ValueError,'N terminus'):
            validate_domain(row['seqres_alpha'],numbered_domains(text)['D'],(SRC/'1ao7_D.fasta').read_text(),text,refs)
    def test_same_construct_observed_split_repaired(self):
        registry=json.loads((ROOT/'docs/fingerprint_paired_identity_20261024/curation_v3/identity_registry.json').read_text())
        rows=[x for x in registry['rows'] if x['complex_id'] in ['5e9d_0','5e9d_1']]
        self.assertEqual(len(rows),2);self.assertTrue(all(x['status']=='admitted' for x in rows))
        self.assertEqual(len({x['receptor_id'] for x in rows}),1)
        self.assertEqual(len({x['construct_id'] for x in rows}),1)
        self.assertEqual(len({x['old_receptor_id'] for x in rows}),2)
    def test_cif_declaration_rejects_gap(self):
        text='loop_\n_entity_poly.entity_id\n_entity_poly.pdbx_strand_id\n_entity_poly.type\n1 A polypeptide(L)\n#\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.hetero\n_entity_poly_seq.mon_id\n_entity_poly_seq.num\n1 n ALA 1\n1 n GLY 2\n#'
        self.assertEqual(cif_declared_sequences(text),{'A':'AG'})
        with self.assertRaises(ValueError):cif_declared_sequences(text.replace('GLY 2','GLY 3'))
    def test_development_fold_isolation_and_control(self):
        import torch
        from scripts.run_paired_identity_development import training_records,validate
        p=validate(ROOT/'docs/fingerprint_paired_identity_20261024/development/protocol.json')
        data=torch.load(ROOT/p['dataset'],map_location='cpu',weights_only=False);records=data['records']
        self.assertEqual(len(records),25);self.assertEqual(len(p['folds']),2)
        for fold in p['folds']:
            fit,val=fold['fit_indices'],fold['validation_indices']
            for key in ['receptor_id','family_id','alpha_variable','beta_variable']:
                self.assertFalse({records[i][key] for i in fit}&{records[i][key] for i in val})
            for seed in p['seeds']:
                shuffled=training_records(records,fit,seed,'identity_preserved_shuffled')
                self.assertEqual([shuffled[i]['receptor_id'] for i in val],[records[i]['receptor_id'] for i in val])
                self.assertEqual(sorted(shuffled[i]['receptor_id'] for i in fit),sorted(records[i]['receptor_id'] for i in fit))
if __name__=='__main__':unittest.main()
