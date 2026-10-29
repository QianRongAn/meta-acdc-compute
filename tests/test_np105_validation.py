from pathlib import Path
import unittest,json,gzip,tempfile,sys,io,hashlib,contextlib
from unittest.mock import patch
from scripts.number_np105_all_cells import strict_boundary
from scripts import fetch_np105_complete
R=Path(__file__).resolve().parents[1];B=R/'docs/fingerprint_np105_validation_20261029';P=R/'docs/fingerprint_np105_raw_20261029'
class Response(io.BytesIO):
 def __init__(self,data,content_range):super().__init__(data);self.status=206;self.headers={'Content-Length':str(len(data)),'Content-Range':content_range}
class ValidationTests(unittest.TestCase):
 def alpha(self):
  r=json.loads((P/'raw_classic_numbering.json').read_text());i=next(i for i,x in enumerate(r['inputs']) if x['role']=='TRA');return r['inputs'][i],r['numbered'][i][0],r['details'][i][0]
 def test_strict_boundary_accepts_observed_alpha(self):
  x,(n,start,end),d=self.alpha();self.assertEqual(strict_boundary(x['ORF'],n,start,end,d,'TRA')['first_IMGT'],'1')
 def test_truncated_FR1_is_never_extended(self):
  x,(n,start,end),d=self.alpha();n=[(p,a) for p,a in n if p[0]>2]
  with self.assertRaisesRegex(ValueError,'endpoints'):strict_boundary(x['ORF'],n,start+2,end,d,'TRA')
 def test_species_fallback_cannot_be_admitted(self):
  x,(n,start,end),d=self.alpha();d=dict(d,species='human')
  with self.assertRaisesRegex(ValueError,'species'):strict_boundary(x['ORF'],n,start,end,d,'TRA')
 def test_low_quality_variant_is_retained_and_quarantined(self):
  r=json.loads((B/'variant_raw_support.json').read_text());self.assertEqual(len(r['targets']),1);x=r['targets'][0];self.assertEqual(x['variant_codon'],'AAC');self.assertEqual(x['variant_codon_target_read_qualities'][0]['qualities'][-1],12);self.assertNotIn('AAC',x['target_cell_Q30_read_codon_counts']);self.assertIn('quarantine',x['decision'])
 def test_all92_are_not_counted_as_independent_receptors(self):
  r=json.loads((B/'cell_summary.json').read_text());self.assertEqual(sum(r['cell_counts'].values()),92);self.assertEqual(r['n_complete_paired_cells_exact_prior_nt'],37);self.assertEqual(r['n_independent_receptors_added'],0);self.assertFalse(r['identity_admitted']);q=json.loads((B/'paired_cells_base_support.json').read_text());self.assertEqual(q['n_independently_assembled_pairs'],37);self.assertEqual(q['n_pairs_with_Q30_support_at_every_base'],33);self.assertEqual(sum(not c['every_base_in_both_roles_has_Q30_support'] for c in q['cells']),4)
 def download(self,data,md5,content_range):
  with tempfile.TemporaryDirectory() as d:
   out=Path(d)/'fresh';argv=['fetch','--out',str(out)]
   with patch.object(sys,'argv',argv),patch.object(fetch_np105_complete,'FILES',[(1,len(data),md5)]),patch.object(fetch_np105_complete.urllib.request,'urlopen',return_value=Response(data,content_range)):
    error=None
    try:
     with contextlib.redirect_stdout(io.StringIO()):fetch_np105_complete.main()
    except Exception as e:error=e
   return error,json.loads((out/'download_receipt.json').read_text())
 def test_whole_file_checksum_is_required(self):
  error,r=self.download(b'ACGT','0'*32,'bytes 0-3/4');self.assertIsInstance(error,AssertionError);self.assertFalse(r['complete']);self.assertFalse(r['files'][0]['complete'])
 def test_wrong_HTTP_range_cannot_pass(self):
  error,r=self.download(b'ACGT',hashlib.md5(b'ACGT').hexdigest(),'bytes 0-3/400');self.assertIsInstance(error,ValueError);self.assertFalse(r['complete'])
 def test_complete_mock_file_passes_size_and_md5(self):
  error,r=self.download(b'ACGT',hashlib.md5(b'ACGT').hexdigest(),'bytes 0-3/4');self.assertIsNone(error);self.assertTrue(r['complete']);self.assertEqual(r['files'][0]['sha256'],hashlib.sha256(b'ACGT').hexdigest())
if __name__=='__main__':unittest.main()
