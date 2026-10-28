import unittest
from scripts.join_native_search import parse_hits,pair_candidates
from scripts.validate_native_search import exact_entity_pairs


class NativeSearchTests(unittest.TestCase):
    def test_failure_cannot_be_reported_as_no_hits(self):
        with self.assertRaises(ValueError):
            parse_hits(500,b'{"status":500}')
        self.assertIsNone(pair_candidates(None,['7N5P_5']))

    def test_true_no_hits_is_distinct(self):
        self.assertEqual(parse_hits(204,b''),[])
        self.assertEqual(pair_candidates([],['7N5P_5']),[])

    def test_local_score_is_not_whole_domain_identity(self):
        alpha = parse_hits(200,b'{"result_type":"polymer_entity","total_count":1,"result_set":[{"identifier":"3PQY_4","score":0.0}]}')
        self.assertEqual(pair_candidates(alpha,['3PQY_5']),['3pqy'])

    def test_incomplete_results_rejected(self):
        with self.assertRaises(ValueError):
            parse_hits(200,b'{"result_type":"polymer_entity","total_count":2,"result_set":[]}')

    def test_local_fragment_not_exact_full_query(self):
        entry={'polymer_entities':[{'rcsb_id':'X_1','entity_poly':{'pdbx_seq_one_letter_code_can':'KTT'}},
                                   {'rcsb_id':'X_2','entity_poly':{'pdbx_seq_one_letter_code_can':'BETA'}}]}
        self.assertEqual(exact_entity_pairs(entry,'DAKTT','BETA'),[])

    def test_single_entity_fusion_does_not_establish_paired_chains(self):
        entry={'polymer_entities':[{'rcsb_id':'X_1','entity_poly':{'pdbx_seq_one_letter_code_can':'ALPHABET'}}]}
        self.assertEqual(exact_entity_pairs(entry,'ALPHA','BET'),[])

    def test_two_complete_distinct_entities_are_discovery_lead(self):
        entry={'polymer_entities':[{'rcsb_id':'X_1','entity_poly':{'pdbx_seq_one_letter_code_can':'TAGALPHA'}},
                                   {'rcsb_id':'X_2','entity_poly':{'pdbx_seq_one_letter_code_can':'BETATAIL'}}]}
        self.assertEqual(len(exact_entity_pairs(entry,'ALPHA','BETA')),1)


if __name__=='__main__':
    unittest.main()
