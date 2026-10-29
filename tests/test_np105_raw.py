from pathlib import Path
import tempfile,unittest,json
from scripts.select_np105_raw_reads import records,trim_ends
from scripts.validate_np105_raw import aligned_bases,parse_cigar
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'docs/fingerprint_np105_raw_20261029'
class RawAuditTests(unittest.TestCase):
 def test_partial_terminal_record_is_never_rescued(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'reads.fastq';p.write_text('@good\nACGT\n+\nIIII\n@partial\nACGT\n+\nII')
   self.assertEqual(list(records(p)),[('@good','ACGT','IIII')])
 def test_malformed_complete_record_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'reads.fastq';p.write_text('@bad\nACGT\n+\nIII\n')
   with self.assertRaises(ValueError):list(records(p))
 def test_trim_never_substitutes_internal_bases(self):
  self.assertEqual(trim_ends('AACNGT','!II!I!'),('ACNG','II!I'))
 def test_cigar_indels_and_softclip_preserve_coordinates(self):
  self.assertEqual(list(aligned_bases(11,'1S2M1I1M1D1M','AACGTT','IIIIII')),[(10,'A',40),(11,'C',40),(12,'T',40),(14,'T',40)])
 def test_invalid_cigar_rejected(self):
  with self.assertRaises(ValueError):parse_cigar('2Mgarbage')
 def test_six_independent_domains_and_no_filled_bases(self):
  r=json.loads((B/'raw_sequence_validation.json').read_text())
  self.assertFalse(r['biological_reference_used_for_assembly']);self.assertEqual(r['bases_substituted'],0)
  for x in r['domains'].values():
   self.assertEqual(len(x['independent_barcodes']),3);self.assertTrue(x['same_nt_in_three_cells']);self.assertTrue(x['all_bases_supported_in_top3_cells']);self.assertTrue(x['non_gap_numbering_identical']);self.assertGreater(x['minimum_Q30_cells'],70);self.assertTrue(x['primary_labeled_variable_identical_nt_in_three_cells']);self.assertTrue(all(y['top3_cells_supported'] and y['Q30_unambiguous_molecules']>0 for y in x['post_numbered_source_base_evidence']))
 def test_ORF_J_cannot_be_promoted_by_productive_source_flag(self):
  x=json.loads((B/'raw_sequence_validation.json').read_text());self.assertEqual(x['domains']['TRB']['longest_exact_J_suffix_matches'][0]['function'],'ORF');self.assertFalse(x['identity_admitted']);self.assertFalse(x['protocol_frozen'])
 def test_every_published_base_has_observed_support(self):
  r=json.loads((B/'raw_sequence_validation.json').read_text())
  for x in r['domains'].values():
   self.assertEqual(len(x['base_evidence']),len(x['nucleotides']))
   self.assertTrue(all(y['Q30_unambiguous_molecules']>0 for y in x['base_evidence']))
if __name__=='__main__':unittest.main()
