import unittest
from scripts.validate_expanded_controls import complete_pair_ready
class ExpandedPairReadinessTests(unittest.TestCase):
    def test_partner_alone_is_not_a_pair(self):self.assertFalse(complete_pair_ready(False,True,True))
    def test_anchor_alone_is_not_a_pair(self):self.assertFalse(complete_pair_ready(True,False,True))
    def test_fragments_are_not_ready(self):self.assertFalse(complete_pair_ready(True,True,False))
    def test_both_complete_fixed_passes(self):self.assertTrue(complete_pair_ready(True,True,True))
