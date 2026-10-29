import unittest
from scripts.audit_np105_context import domain_context, outcome
class ContextAuditTests(unittest.TestCase):
    def test_original_interval_is_preserved(self):
        self.assertEqual(domain_context(dict(ORF="NABCQ",boundary=dict(start0=1,end0_exclusive=4,sequence="ABC"))),("N","ABC","Q"))
    def test_nonunique_domain_is_rejected(self):
        with self.assertRaises(ValueError):
            domain_context(dict(ORF="ABCABC",boundary=dict(start0=0,end0_exclusive=3,sequence="ABC")))
    def test_wrong_offset_is_rejected(self):
        with self.assertRaises(ValueError):
            domain_context(dict(ORF="NABCQ",boundary=dict(start0=0,end0_exclusive=3,sequence="ABC")))
    def test_failure_and_wrong_chain_are_distinct(self):
        self.assertNotEqual(outcome(dict(neural=dict(numbering=None,error="cutoff",chain_type="F"))),outcome(dict(neural=dict(numbering=None,error=None,chain_type="A"))))
    def test_score_or_query_offset_alone_not_new_numbering_outcome(self):
        a=dict(numbering=[[[1," "],"A"],[[2," "],"-"]],chain_type="A",error=None,score=80,query_start=3)
        b=dict(a,score=90,query_start=5)
        self.assertEqual(outcome(dict(neural=a)),outcome(dict(neural=b)))
if __name__=="__main__": unittest.main()
