import unittest
from scripts.audit_np_adjudication import unambiguous_pair,suffix_length
class NaturalPairEvidenceTests(unittest.TestCase):
    def setUp(self):
        common=dict(is_cell='True',high_confidence='True',full_length='True',productive='True',barcode='cell1',raw_clonotype_id='clone1')
        self.a=dict(common,chain='TRA',cdr3='CA',cdr3_nt='TGTGCT',j_gene='TRAJ23',c_gene='TRAC')
        self.b=dict(common,chain='TRB',cdr3='CF',cdr3_nt='TGTTTT',j_gene='TRBJ1-6',c_gene='TRBC1')
    def check(self,rows):return unambiguous_pair(rows,'CA','CF')
    def test_strict_same_cell_pair_supported(self):self.assertIsNotNone(self.check([self.a,self.b]))
    def test_missing_beta_constant_stays_unresolved_in_strict_tier(self):
        b=dict(self.b,c_gene='None');self.assertIsNone(self.check([self.a,b]));self.assertIsNotNone(unambiguous_pair([self.a,b],'CA','CF',require_beta_constant=False))
    def test_extra_productive_alpha_is_ambiguous(self):self.assertIsNone(self.check([self.a,self.b,dict(self.a,cdr3='CC',cdr3_nt='TGTTGT')]))
    def test_different_clonotype_does_not_establish_pair(self):self.assertIsNone(self.check([self.a,dict(self.b,raw_clonotype_id='clone2')]))
    def test_wrong_nucleotide_translation_is_rejected(self):self.assertIsNone(self.check([self.a,dict(self.b,cdr3_nt='TGTTGT')]))
    def test_different_cell_barcodes_are_rejected(self):self.assertIsNone(self.check([self.a,dict(self.b,barcode='cell2')]))
    def test_wrong_alpha_constant_does_not_establish_alpha_role(self):self.assertIsNone(self.check([dict(self.a,c_gene='TRDC'),self.b]))
    def test_exact_suffix_stops_at_first_mismatch(self):self.assertEqual(suffix_length('CAMREFFG','CASSQFFG'),3)
