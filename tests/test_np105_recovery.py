from pathlib import Path
import unittest,tempfile,hashlib,io,json
from unittest.mock import patch
from types import SimpleNamespace
from scripts.resume_np105_native import assess,fingerprint
from scripts.audit_np105_whole_variants import haplotype_umis
from scripts.number_np105_whole_neural import occupied_numbering
from scripts.validate_np105_whole_pairs import unique_canonical_pair
from scripts.select_np105_full_sra import read_fastq,selected_pair,check_destination,SIZE
R=Path(__file__).resolve().parents[1];B=R/'docs/fingerprint_np105_recovery_20261029'
class RecoveryTests(unittest.TestCase):
 def test_haplotype_rejects_conflicting_or_missing_bases(self):
  calls={'good':{1:{'A'},2:{'C'}},'conflict':{1:{'A','G'},2:{'C'}},'missing':{1:{'A'}}};self.assertEqual(haplotype_umis(calls,[1,2],'AC'),{'good'})
 def test_separate_molecules_cannot_fill_one_haplotype(self):
  self.assertEqual(haplotype_umis({'one':{1:{'A'}},'two':{2:{'C'}}},[1,2],'AC'),set())
 def test_failed_neural_numbering_is_retained_as_no_map(self):
  self.assertEqual(occupied_numbering(dict(numbering=None,chain_type='F',error='Score less than cut off.')),{})
 def test_single_exact_pair_is_distinct_from_ambiguous_assemblies(self):
  cell={role:[dict(complete=True,exact_previous_nt=True)] for role in ['TRA','TRB']};cell.update(TRA_complete_unique_nt=1,TRB_complete_unique_nt=1);self.assertTrue(unique_canonical_pair(cell));cell['TRA_complete_unique_nt']=2;cell['TRA'].append(dict(complete=True,exact_previous_nt=False));self.assertFalse(unique_canonical_pair(cell))
 def test_missing_full_chain_cannot_pass_canonical_pair_check(self):
  cell=dict(TRA=[],TRB=[dict(complete=True,exact_previous_nt=True)],TRA_complete_unique_nt=1,TRB_complete_unique_nt=1);self.assertFalse(unique_canonical_pair(cell))
 def test_guest_filesystem_cannot_be_used_for_full_export(self):
  with self.assertRaisesRegex(ValueError,'D volume'):check_destination(Path('/home/administrator/unsafe-output'))
 def test_D_free_space_must_pass_before_export(self):
  with patch('scripts.select_np105_full_sra.shutil.disk_usage',return_value=SimpleNamespace(free=SIZE)):
   with self.assertRaisesRegex(ValueError,'free space'):check_destination(Path('/mnt/d/ChatGPT/fresh-output'))
 def test_D_with_sufficient_space_passes(self):
  with patch('scripts.select_np105_full_sra.shutil.disk_usage',return_value=SimpleNamespace(free=9*SIZE)):self.assertEqual(check_destination(Path('/mnt/d/ChatGPT/fresh-output')),9*SIZE)
 def test_complete_file_needs_whole_md5_and_preserved_prefix(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'data';p.write_bytes(b'ACGTACGT');x=assess(p,4,hashlib.sha256(b'ACGT').hexdigest(),8,hashlib.md5(b'ACGTACGT').hexdigest());self.assertTrue(x['complete'])
 def test_changed_preserved_prefix_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'data';p.write_bytes(b'TCGTACGT');x=assess(p,4,hashlib.sha256(b'ACGT').hexdigest(),8,hashlib.md5(b'TCGTACGT').hexdigest());self.assertFalse(x['complete']);self.assertFalse(x['prefix_unchanged'])
 def test_prefix_is_not_a_whole_file(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'data';p.write_bytes(b'ACGT');self.assertFalse(assess(p,4,hashlib.sha256(b'ACGT').hexdigest(),8,hashlib.md5(b'ACGTACGT').hexdigest())['complete'])
 def test_wrong_whole_file_checksum_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'data';p.write_bytes(b'ACGTACGT');self.assertFalse(assess(p,4,hashlib.sha256(b'ACGT').hexdigest(),8,'0'*32)['complete'])
 def test_partial_source_record_is_never_repaired(self):
  with self.assertRaises(ValueError):read_fastq(io.StringIO('@SRR34543498.1\nACGT\n+\nII'))
 def test_pair_identifier_mismatch_rejected(self):
  with self.assertRaises(ValueError):selected_pair(('@one','A'*150,'I'*150),('@two','A'*150,'I'*150),{'A'*16})
 def test_selection_preserves_observed_bases(self):
  one=('@SRR34543498.1','A'*16+'C'*10+'GT'*62,'I'*150);two=('@SRR34543498.1','TG'*75,'I'*150);x,reason=selected_pair(one,two,{'A'*16});self.assertIsNone(reason);self.assertEqual(x[1],'C'*10);self.assertEqual(x[2],[('GT'*62,'I'*124),('TG'*75,'I'*150)])
 def test_paper_generation_source_not_structure_study_or_native_pair_proof(self):
  r=json.loads((B/'tcr1414_source_audit.json').read_text());self.assertTrue(r['all_four_chain_queries_successful']);self.assertTrue(r['receptor_generation_DOI_not_imputed_as_structure_DOI']);self.assertEqual(r['reported_chain_candidates']['tested_combinations'],15)
  for x in r['receptors']:self.assertIsNone(x['structure_publication_DOI']);self.assertFalse(x['native_single_cell_pairing_established_by_reported_method']);self.assertFalse(x['full_variable_sequence_identity_across_papers_established'])
if __name__=='__main__':unittest.main()
