"""Fresh-reserve tests: exclusion, range role, partial-cohort denial and identity."""
import copy,json,unittest
from meta_acdc.data.fresh_spr_reserve import ROOT,read,verify,assert_development_safe,require_acceptance_ready
from meta_acdc.data.independent_spr import validate_partition
class FreshSPRReserveTests(unittest.TestCase):
 def test_partial_reserve_cannot_be_scored(self):
  result=verify();self.assertFalse(result['evaluation_allowed']);self.assertEqual(result['counts']['lineages'],1)
  with self.assertRaisesRegex(ValueError,'new acceptance blocked'):require_acceptance_ready()
 def test_shared_author_validation_is_excluded(self):
  parent=json.loads((ROOT/'docs/affinity_development_20261018/development_input.json').read_text());x=next(x for x in parent if x['evidence_id']=='ATLAS:485')
  with self.assertRaisesRegex(ValueError,'reserved study source'):assert_development_safe([x])
 def test_reserved_family_cannot_reenter(self):
  dev=json.loads((ROOT/'docs/affinity_development_20261018/development_source_safe.json').read_text());x=dict(dev[0]);x.update({k:v for k,v in read('reserved_measurements.json')[0].items() if k in ('alpha_variable','beta_variable','cdr3_alpha','cdr3_beta')})
  with self.assertRaisesRegex(ValueError,'reserved receptor family'):assert_development_safe([x])
 def test_reserved_peptide_cannot_reenter(self):
  dev=json.loads((ROOT/'docs/affinity_development_20261018/development_source_safe.json').read_text());x=dict(dev[0],peptide='ILKEPVHGVYY')
  with self.assertRaisesRegex(ValueError,'reserved peptide family'):assert_development_safe([x])
 def test_outside_range_mutant_cannot_be_primary(self):
  rows=copy.deepcopy(read('reserved_measurements.json'));rows[1]['metric_role']='primary'
  with self.assertRaisesRegex(ValueError,'range extrapolation'):validate_partition(rows,[],read('study_sources.json'),[],[])
 def test_two_measurements_do_not_make_two_receptors(self):
  rows=read('reserved_measurements.json');self.assertEqual(len(rows),2);self.assertEqual(len({(x['alpha_variable'],x['beta_variable']) for x in rows}),1);self.assertEqual(len({x['lineage'] for x in rows}),1)
 def test_development_rows_keep_parent_fold_and_labels(self):
  parent={x['evidence_id']:x for x in json.loads((ROOT/'docs/affinity_development_20261018/development_input.json').read_text())};dev=json.loads((ROOT/'docs/affinity_development_20261018/development_source_safe.json').read_text())
  self.assertEqual(set(parent)-{x['evidence_id'] for x in dev},{'ATLAS:485'})
  for x in dev:self.assertEqual(x,parent[x['evidence_id']])
if __name__=='__main__':unittest.main()
