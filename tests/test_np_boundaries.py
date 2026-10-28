import unittest
from scripts.audit_np_boundaries import checked_numbering
class NumberedBoundaryTests(unittest.TestCase):
    def call(self,full='GGACDTT',numbering=None,start=2,end=4,chain='A'):
        return checked_numbering(full,numbering or [([1,' '],'A'),([2,' '],'C'),([128,' '],'D')],start,end,chain,'A')
    def test_exact_interval_preserves_extensions(self):
        r=self.call();self.assertEqual((r['start0'],r['end0_exclusive'],r['sequence']),(2,5,'ACD'))
    def test_inclusive_end_is_checked(self):
        with self.assertRaises(ValueError):self.call(end=5)
    def test_duplicated_substring_is_ambiguous(self):
        with self.assertRaises(ValueError):self.call(full='GGACDACDTT')
    def test_role_conflict_cannot_be_relabelled(self):
        with self.assertRaises(ValueError):self.call(chain='D')
    def test_incomplete_terminal_number_is_rejected(self):
        with self.assertRaises(ValueError):self.call(numbering=[([1,' '],'A'),([2,' '],'C'),([120,' '],'D')])
    def test_duplicate_positions_are_rejected(self):
        with self.assertRaises(ValueError):self.call(numbering=[([1,' '],'A'),([1,' '],'C'),([128,' '],'D')])

from scripts.extract_np_final_sources import paired_repertoire_support
class PairedRepertoireTests(unittest.TestCase):
    def setUp(self):self.target=dict(CDR3_alpha='AAA',CDR3_beta='BBB',V_alpha='VA',J_alpha='JA',V_beta='VB',J_beta='JB')
    def test_beta_only_does_not_establish_pair(self):self.assertFalse(paired_repertoire_support(dict(self.target,CDR3_alpha='CCC'),self.target))
    def test_same_loops_with_different_V_does_not_establish_context(self):self.assertFalse(paired_repertoire_support(dict(self.target,V_alpha='VD'),self.target))
