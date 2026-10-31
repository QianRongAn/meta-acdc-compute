"""Regression checks against downloaded public candidate sources; no model inputs."""
import unittest
from scripts.audit_control_closure import audit
class ControlClosureTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.report=audit()
 def test_exact_substrings_do_not_follow_local_identity_score(self):
  r=self.report;self.assertEqual(r['n_fresh_metadata_PDBs'],24);self.assertEqual(r['n_exact_fragment_pair_leads'],20);self.assertEqual(r['n_local_alignment_false_pair_leads'],4)
  self.assertEqual({x['pdb'] for x in r['records'] if x['local_alignment_only']},{'1QRN','1QSE','1QSF','4PRI'})
 def test_sequence_rescue_exists_but_has_no_class_I_components(self):
  self.assertEqual(self.report['new_external_sequence_gate_pass_PDBs'],['7AMP','7AMR']);keys=set()
  for x in self.report['records']:
   if x['pdb'] in ['7AMP','7AMR']:
    self.assertTrue(x['component_inventory']['Fab']);self.assertFalse(x['component_inventory']['explicit_class_I']);self.assertFalse(x['component_inventory']['beta2m'])
    for c in x['external_sequence_checks']:
     self.assertTrue(c['complete_external_sequence_gate_pass']);keys.add(c['paired_variable_sequence_key']);self.assertFalse(c['identity_admitted'])
  self.assertEqual(len(keys),1)
 def test_TK3_candidate_is_not_promoted_without_numbering(self):
  x=next(x for x in self.report['records'] if x['pdb']=='4PRH');self.assertTrue(x['exact_distinct_entity_fragment_pair']);self.assertTrue(x['component_inventory']['explicit_class_I']);self.assertTrue(x['component_inventory']['beta2m']);self.assertTrue(x['external_annotation']['explicit_missing_entry']);self.assertEqual(x['external_annotation']['domains'],0);self.assertFalse(x['identity_admitted'])
 def test_NP_missing_entry_pages_are_not_valid_numbering(self):
  for x in self.report['np_external_numbering'].values():self.assertTrue(x['explicit_missing_entry']);self.assertEqual(x['domains'],0);self.assertEqual(x['explicit_pairs'],[])
 def test_class_II_and_unbound_not_used_as_class_I(self):
  records={x['pdb']:x for x in self.report['tcr1414_scope']};self.assertEqual(records['9XC3']['n_explicit_class_II_entities'],2);self.assertEqual(records['9X9I']['n_explicit_class_II_entities'],2);self.assertTrue(all(not x['qualifies_current_class_I'] for x in records.values()))
 def test_no_release_or_acceptance_access(self):
  for key in ['new_identity_admissions','new_QC_cross_study_pairs','new_native_same_source_same_V_groups','protocol_frozen','training_started','independent_acceptance_accessed']:self.assertFalse(self.report[key])
 def test_unchanged_revisions_do_not_trigger_coordinate_refetch(self):
  self.assertEqual(len(self.report['coordinate_revision_recheck']),2)
  for x in self.report['coordinate_revision_recheck']:self.assertFalse(x['changed_since_previous_coordinate_recheck']);self.assertFalse(x['coordinate_download_repeated'])
 def test_full_input_TK3_numbering_replays_without_sequence_edits(self):
  diagnosis=self.report['independent_4PRH_numbering'];self.assertFalse(diagnosis['inputs_trimmed']);self.assertFalse(diagnosis['imputation']);self.assertEqual(len(diagnosis['rows']),2)
  for row in diagnosis['rows']:self.assertTrue(row['two_method_positions_exact']);self.assertTrue(row['functional_human_V_matches']);self.assertFalse(row['identity_admitted'])
 def test_TK3_alpha_J_failure_is_not_hidden_by_N_terminal_rescue(self):
  rows={r['role']:r for r in self.report['independent_4PRH_numbering']['rows']};self.assertEqual(rows['alpha']['last10'],'GEGTQLTVNP');self.assertEqual(rows['alpha']['functional_human_J_matches'],[]);self.assertEqual(rows['alpha']['all_function_status_human_or_mouse_J_exact_matches'],[]);self.assertEqual(rows['beta']['functional_human_J_matches'],['TRBJ2-2*01'])
if __name__=='__main__':unittest.main()
