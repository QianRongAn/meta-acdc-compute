import unittest,json,hashlib
from pathlib import Path
from scripts.summarize_cross_source_extension import qualified_rows,control_groups
ROOT=Path(__file__).resolve().parents[1]
class ExtensionChecks(unittest.TestCase):
 def row(self,pdb='a',rid='one',study='doi:x',pmhc='p'):
  return dict(pdb=pdb,alpha_chain='A',beta_chain='B',status='identity_admitted_structure_QC_pending',receptor_id=rid,source_study_key=study,V_context='human|TRAV1|TRBV1',unit_metadata=dict(pmhc_construct_key=pmhc,peptide='PEP'))
 def qc(self,pdb='a',rid='one',status='native_graph_QC_pass'):
  return dict(pdb=pdb,alpha_chain='A',beta_chain='B',receptor_id=rid,status=status,graph_sha256='hash',unit_metadata=dict(pmhc_construct_key='p',peptide='PEP'))
 def test_qc_identity_mismatch_rejected(self):
  with self.assertRaises(ValueError):qualified_rows([self.row()],[self.qc(rid='other')])
 def test_duplicate_key_rejected(self):
  with self.assertRaises(ValueError):qualified_rows([self.row(),self.row()],[self.qc()])
 def test_flagged_unit_excluded(self):self.assertEqual(qualified_rows([self.row()],[self.qc(status='quarantined')]),[])
 def test_different_source_is_not_same_source_control(self):self.assertEqual(control_groups([self.row(),self.row('b','two','doi:y')]),[])
 def test_same_identity_is_not_negative(self):self.assertEqual(control_groups([self.row(),self.row('b')]),[])
 def test_pmhc_match_and_lineage_not_assumed(self):
  g=control_groups([self.row(),self.row('b','two')])[0];self.assertTrue(g['pairs'][0]['same_complete_pmhc_construct']);self.assertFalse(g['independent_clonotypes_proven'])
 def test_actual_jdi_control_and_blocked_release(self):
  d=json.loads((ROOT/'docs/fingerprint_cross_source_extension_20261024/joint_readiness.json').read_text());self.assertEqual(d['summary']['n_QC_cross_study_pairs'],0);self.assertEqual(d['summary']['n_QC_source_balanced_blocks'],0);self.assertFalse(d['training_release_created']);g=d['same_study_same_V_controls'][0];self.assertEqual(set(g['pdbs']),{'7ow5','7ow6','7pb2'});self.assertTrue(any(p['same_complete_pmhc_construct'] and set([p['pdb_a'],p['pdb_b']])=={'7ow6','7pb2'} for p in g['pairs']))
 def test_actual_index_mapping(self):
  d=json.loads((ROOT/'docs/fingerprint_cross_source_extension_20261024/legacy_mapping_review_final/records_index_exact.json').read_text());self.assertEqual([r['graph_index'] for r in d['records']],list(range(331)));self.assertEqual(d['summary']['n_legacy_observed_keys_crossing_folds'],0);self.assertEqual(d['summary']['n_admitted_canonical_ids_crossing_legacy_folds'],1)
 def test_attachment_multipart_integrity(self):
  base=ROOT/'docs/fingerprint_cross_source_extension_20261024/sources/attachment_probe';m=json.loads((base/'multipart_manifest.json').read_text());raw=b''.join((base/x['file']).read_bytes() for x in m['parts']);self.assertEqual(len(raw),m['original_bytes']);self.assertEqual(hashlib.sha256(raw).hexdigest(),m['original_sha256'])
 def test_actual_confidence_mask_is_both_domains(self):
  d=json.loads((ROOT/'docs/fingerprint_cross_source_extension_20261024/legacy_mapping_review_final/records_index_exact.json').read_text())
  for v in d['summary']['confidence_column_by_domain'].values():self.assertEqual(v['raw_unique'],[0.5]);self.assertEqual(v['normalized_min'],0.0);self.assertEqual(v['normalized_max'],0.0)
if __name__=='__main__':unittest.main()
