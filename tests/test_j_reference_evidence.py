"""Unfiltered audit evidence must retain stops without changing admission rules."""
from pathlib import Path
import json,unittest
from meta_acdc.structure.j_reference_evidence import read_j_evidence,translate_reference,suffix_evidence
from meta_acdc.structure.paired_variable_identity import germline_references
from scripts.audit_np_boundaries import raw_j_references
from scripts.audit_tk3_j_evidence import audit
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_tk3_j_evidence_20261101'
class JReferenceEvidenceTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.paths=sorted((ROOT/'docs/fingerprint_paired_identity_20261024/sources').glob('*_TRAJ.fasta'))+sorted((ROOT/'docs/fingerprint_paired_identity_20261024/sources').glob('*_TRBJ.fasta'));cls.records=read_j_evidence(cls.paths);cls.report=audit()
 def test_germline_stop_is_retained(self):
  r=next(x for x in self.records if x['gene']=='TRAJ58*01' and x['species']=='Homo sapiens');self.assertEqual(r['functionality'],'ORF');self.assertEqual(r['protein'],'*ETSGSRLTFGEGTQLTVNP');self.assertEqual(r['stop_positions0'],[0]);self.assertFalse(any(x['gene']=='TRAJ58*01' and '|Homo sapiens|' in x['header'] for x in raw_j_references(self.paths)))
 def test_TK3_tail_matches_downstream_of_stop(self):
  r=next(x for x in self.report['TK3_cases'] if x['role']=='alpha');self.assertEqual(r['old_reported_all_status_matches'],[]);self.assertEqual(len(r['corrected_all_status_suffix_evidence']),1);m=r['corrected_all_status_suffix_evidence'][0];self.assertEqual(m['gene'],'TRAJ58*01');self.assertEqual(m['longest_exact_suffix_length'],18);self.assertTrue(m['all_stops_before_matching_suffix']);self.assertFalse(m['identity_admitted'])
 def test_existing_F_gate_is_not_relaxed(self):
  refs=germline_references(self.paths);self.assertNotIn(('human','TRAJ58*01'),refs);self.assertFalse(self.report['admission_gates_changed']);self.assertEqual(self.report['new_identity_admissions'],0)
 def test_ambiguous_codon_retained_only_as_evidence(self):
  h='test|TRAJ999*01|Homo sapiens|ORF|J-REGION|1..33|33 nt|1|';r=translate_reference(h,'NNN'+'GCT'*10);self.assertEqual(r['protein'],'X'+'A'*10);self.assertEqual(r['ambiguous_positions0'],[0]);m=suffix_evidence('A'*10,[r],'Homo sapiens','TRAJ')[0];self.assertTrue(m['old_stopfree_filter_would_drop']);self.assertFalse(m['identity_admitted'])
 def test_stop_in_observed_variable_is_rejected(self):
  with self.assertRaises(ValueError):suffix_evidence('A'*10+'*',self.records,'Homo sapiens','TRAJ')
 def test_malformed_header_and_nucleotides_fail(self):
  with self.assertRaises(ValueError):translate_reference('invalid','ATG')
  with self.assertRaises(ValueError):translate_reference('test|TRAJ1*01|Homo sapiens|F|J-REGION|1..3|3 nt|4|','ATG')
  with self.assertRaises(ValueError):translate_reference('test|TRAJ1*01|Homo sapiens|F|J-REGION|1..3|3 nt|1|','<html>')
 def test_short_fragment_is_not_tail_evidence(self):self.assertEqual(suffix_evidence('G'*9,self.records,'Homo sapiens','TRAJ'),[])
 def test_bulk_impact_does_not_mean_admission(self):
  self.assertEqual(self.report['n_observed_external_domains'],778);self.assertEqual(self.report['n_domains_with_newly_visible_references'],7);self.assertEqual(set(self.report['affected_PDBs']),{'3mv7','3mv8','3mv9','4pri','4prp','7n2n','7n2o'});self.assertEqual(self.report['affected_genes'],['TRAJ58*01']);self.assertEqual(self.report['new_QC_cross_study_pairs'],0)
 def test_primary_database_conflict_stays_explicit(self):
  p=json.loads((BASE/'primary_source_audit_final.json').read_text());self.assertEqual(p['primary_paper']['reported_V_J']['alpha'],['TRAV20','TRAJ58']);self.assertEqual(p['IEDB_conflict']['calculated_alpha_J'],'TRAJ57*01');self.assertEqual(p['IEDB_conflict']['annotation_status'],'source_conflict_quarantined');self.assertFalse(p['IEDB_conflict']['external_record_edited'])
 def test_clone_name_homonyms_do_not_count_as_TCR_controls(self):
  p=json.loads((BASE/'primary_source_audit_final.json').read_text());self.assertEqual({x['pdb'] for x in p['natural_clone_search']['excluded_homonyms']},{'4JEG','9E7K'});self.assertFalse(p['new_natural_control_groups']);self.assertFalse(p['protocol_frozen']);self.assertFalse(p['independent_acceptance_accessed'])
if __name__=='__main__':unittest.main()
