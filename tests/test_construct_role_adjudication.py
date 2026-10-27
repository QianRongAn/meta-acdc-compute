import unittest
from scripts.adjudicate_construct_sources import linked_roles
from meta_acdc.structure.primary_role_evidence import classify_role


class ConstructRoleEvidenceTests(unittest.TestCase):
    def test_orthogonal_primary_motifs_identify_distinct_chains(self):
        self.assertEqual(linked_roles({'D': 'XXALPHAMOTIFYY', 'E': 'ZZBETAMOTIF'},
                                     ['ALPHA', 'MOTIFYY'], ['BETA']), {'alpha': 'D', 'beta': 'E'})

    def test_ambiguous_published_motif_is_rejected(self):
        with self.assertRaises(ValueError):
            linked_roles({'D': 'SHARED', 'E': 'SHARED'}, ['SHARED'], ['SHARED'])

    def test_one_chain_cannot_supply_both_roles(self):
        with self.assertRaises(ValueError):
            linked_roles({'D': 'ALPHABET'}, ['ALPHA'], ['BET'])

    def test_missing_sequence_anchor_rejected(self):
        with self.assertRaises(ValueError):
            linked_roles({'D': 'A', 'E': 'B'}, [], ['B'])

    def test_general_conflict_gate_remains_strict(self):
        self.assertEqual(classify_role('T-cell receptor alpha', 'T-cell receptor beta'),
                         'conflicting_alpha_beta')


if __name__ == '__main__':
    unittest.main()
