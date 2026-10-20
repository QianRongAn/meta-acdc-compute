"""Mandatory checks of retained evidence, not scientific acceptance tests."""
import csv,hashlib,importlib.util,json,statistics,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'docs/forensic_followup_20261021'
spec=importlib.util.spec_from_file_location('panel_audit',D/'audit_pending_panel.py')
panel=importlib.util.module_from_spec(spec);spec.loader.exec_module(panel)
class ForensicFollowupTests(unittest.TestCase):
 def test_panel_all_original_bytes_and_counts(self):
  d=json.loads((D/'pending_panel.json').read_text());self.assertEqual(len(d['rows']),347)
  for r in d['rows']:
   data=(ROOT/r['path']).read_bytes()
   self.assertEqual(hashlib.sha256(data).hexdigest(),r['sha256'])
   self.assertEqual(panel.parse(data.decode()),r['sequences'])
   self.assertTrue(r['matches_first_commit'])
  self.assertEqual(d['counts']['exact_receptor_pairs'],98)
  self.assertEqual(d['counts']['pdb_prefixes'],111)
  self.assertEqual(d['counts']['beta_starts_NA'],53)
  self.assertEqual(d['counts']['mhc_starts_GS'],279)
  self.assertTrue(d['all_bytes_match_first_commit'])
  self.assertTrue(d['manifest']['matches_fasta_job_set'])
  self.assertTrue(all(v['exact_inventory_match'] for v in d['copy_comparisons'].values()))
 def test_parser_fails_instead_of_merging_or_guessing_roles(self):
  for bad in ['ACD','>alpha\nACD\n>alpha\nDEF','>alpha\nACD']:
   with self.assertRaises(ValueError):panel.parse(bad)
 def test_family_sensitivity_is_explicit_not_independence(self):
  d=json.loads((D/'pending_panel.json').read_text())
  self.assertFalse(d['biological_independence_established'])
  pairs=[(r['alpha'],r['beta']) for r in d['pair_index']]
  for method,count in [('equal_length_hamming',89),('matching_blocks_max_length',57)]:
   computed=panel.components(pairs,method)
   self.assertEqual(computed,d['clustering_sensitivity'][method]);self.assertEqual(computed['n_components'],count)
 def test_probe_three_disputed_jobs_match_all_seeds(self):
  d=json.loads((D/'legacy_probe.json').read_text())
  self.assertEqual((d['matches'],d['comparisons']),(15,15))
  self.assertFalse(d['clinical_interpretation_allowed']);self.assertFalse(d['historical_causal_attestation'])
 def test_full_replay_counts_rounding_and_coverage(self):
  d=json.loads((D/'legacy_full_replay.json').read_text())
  self.assertEqual((d['matches'],d['comparisons']),(938,940))
  self.assertEqual(len(d['input_sha256']),884)
  self.assertEqual(len({r['job_id'] for r in d['results']}),188)
  self.assertEqual(len(d['checkpoints']),5)
  for r in d['results']:
   self.assertEqual(r['replay_n_models'],r['historical_n_models']);self.assertTrue(r['model_coverage_matches'])
   mean=statistics.mean(float(f'{v:.4f}') for v in r['per_af3_model_scores'])
   self.assertEqual(f'{mean:.4f}',r['replay_mean'])
   self.assertLessEqual(abs(float(r['replay_mean'])-float(r['historical_mean'])),.000100001)
  self.assertEqual(d['record_status'],'FORENSIC_WITHDRAWN')
  self.assertFalse(d['clinical_interpretation_allowed']);self.assertFalse(d['historical_causal_attestation'])
 def test_original_per_structure_tables_and_replay(self):
  d=json.loads((D/'per_structure_reconciliation.json').read_text())
  full=json.loads((D/'legacy_full_replay.json').read_text())
  self.assertEqual((d['matches'],d['comparisons']),(4417,4420));self.assertEqual(d['uncovered_historical_rows'],0)
  n=0;bad=0
  for name,meta in d['tables'].items():
   path=D/'recovered_per_structure'/name
   self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),meta['sha256'])
   with path.open() as f:old={r['job_id']:r for r in csv.DictReader(f,delimiter='\t')}
   seed=int(name.removeprefix('da_scores_seed').removesuffix('.tsv'))
   for r in full['results']:
    if r['seed']!=seed:continue
    paths=sorted(p for p in full['input_sha256'] if Path(p).stem.startswith(r['job_id']+'_model_'))
    for p,value in zip(paths,r['per_af3_model_scores']):
     n+=1;prior=old[Path(p).stem]['score'];bad+=f'{value:.4f}'!=prior
     self.assertLessEqual(abs(float(f'{value:.4f}')-float(prior)),.000100001)
  self.assertEqual((n,bad),(4420,3));self.assertFalse(d['clinical_interpretation_allowed'])
 def test_failed_21_candidate_retained_without_promotion(self):
  d=json.loads((D/'legacy_21_probe.json').read_text())
  self.assertEqual((d['matches'],d['comparisons']),(0,3))
  self.assertEqual(d['record_status'],'FORENSIC_WITHDRAWN_FAILED_CANDIDATE')
  self.assertFalse(d['clinical_interpretation_allowed'])
if __name__=='__main__':unittest.main()
