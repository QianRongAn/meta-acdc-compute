"""Incident reproductions using committed real structure and recovered score tables.

Fixtures are mandatory. Missing files fail rather than skip these checks.
"""
import csv,copy,hashlib,json,statistics,tempfile,unittest
from pathlib import Path
from dataclasses import replace
from meta_acdc.structure.artifact_status import ROOT,assert_not_withdrawn,validate_score_receipt
D=ROOT/'docs/delivery_audit_20261018'
FIX=ROOT/'tests/fixtures/delivery_audit'
def table(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
class DeliveryAuditTests(unittest.TestCase):
 def test_real_five_chain_identity_matches_official_pdb(self):
  from meta_acdc.structure.cif import parse_cif
  from meta_acdc.structure.identity import chain_sequence
  seal=json.loads((FIX/'manifest.json').read_text())
  for name,h in seal['files'].items():self.assertEqual(hashlib.sha256((FIX/name).read_bytes()).hexdigest(),h)
  refs={}
  for rec in (FIX/'5BRZ.fasta').read_text().split('>')[1:]:
   h,*lines=rec.splitlines();refs[''.join(lines)]=h
  by={}
  for x in parse_cif(FIX/'5brz_EVDPIGHLY_model_0.cif'):by.setdefault(x.chain,[]).append(x)
  self.assertEqual(sorted(len(x) for x in by.values()),[9,100,197,241,275])
  for rs in by.values():self.assertIn(chain_sequence(rs),refs)
  self.assertFalse(seal['clinical_a3a_equivalent'])
 def test_recovered_all_five_seeds_reproduce_all_seven_clinical_rows(self):
  seeds=[{x['job_id']:float(x['score']) for x in table(D/'recovered'/f'da_scores_seed{i}_ensemble.tsv')} for i in range(5)]
  for x in table(D/'historical/data/processed/kn11_clinical_scan.tsv'):
   job=x['pdb']+'_'+x['peptide'];self.assertEqual(f'{statistics.mean(s[job] for s in seeds):.4f}',x['score'])
 def test_dispersion_counts_and_four_qc_flags_reproduce(self):
  rows=table(D/'historical/data/processed/prediction_scores_ensemble.relabeled.tsv');sd=[float(x['std']) for x in rows]
  self.assertEqual(len(sd),21);self.assertEqual(statistics.median(sd),.2512);self.assertEqual(sum(x>.2 for x in sd),13)
  audit=json.loads((D/'clinical_table_reconciliation.json').read_text());self.assertEqual(sum(bool(x['qc_flag']) for x in audit),4)
 def test_current_clinical_table_has_no_score_or_hazard(self):
  for x in table(ROOT/'data/processed/kn11_clinical_scan.tsv'):
   self.assertEqual(x['record_status'],'WITHDRAWN');self.assertEqual(x['clinical_interpretation_allowed'],'no')
   for k in ['score','target_score','delta_vs_target','hazard_flag']:self.assertEqual(x[k],'')
 def test_reference_generator_separates_constructs_and_never_labels_risk(self):
  from meta_acdc.data.clinical_cases import load_references
  rows=load_references();by={x['record_id']:x for x in rows}
  self.assertEqual(by['structural-MAG-IC3']['structural_receptor'],'MAG-IC3')
  self.assertEqual(by['clinical-a3a']['pdb'],'')
  self.assertTrue(all(not x['clinical_prediction_label'] for x in rows))
 def test_bad_pmids_removed_from_active_references(self):
  text=(ROOT/'data/processed/clinical_gold_standard.tsv').read_text()
  for bad in ['23999400','23999405','23470321','34525232']:self.assertNotIn(bad,text)
  for good in ['23926201','26758806','23377668','34551229']:self.assertIn(good,text)
 def test_legacy_scores_blocked_even_under_a_new_filename(self):
  original=D/'recovered/prediction_scores_ensemble.tsv'
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'fresh.tsv';p.write_bytes(original.read_bytes())
   with self.assertRaisesRegex(ValueError,'withdrawn historical'):assert_not_withdrawn(p)
 def test_source_receipt_missing_or_stale_is_rejected(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'scores.tsv';p.write_text('job_id\tscore\tqc\nx\t.6\tqc-pass\n')
   with self.assertRaisesRegex(ValueError,'missing score provenance'):validate_score_receipt(p,'abc')
   Path(str(p)+'.provenance.json').write_text(json.dumps(dict(output_sha256='wrong')))
   with self.assertRaisesRegex(ValueError,'hash mismatch'):validate_score_receipt(p,'abc')
 def test_constant_seed_retained_and_std_axis_explicit(self):
  from meta_acdc.structure.ensemble_mean import combine
  rows,diag=combine([dict(a=.1,b=.8),dict(a=.5,b=.5)],['a','b'])
  self.assertEqual(diag['degenerate_seeds_retained'],[1]);self.assertTrue(all(x['n_seeds']==2 and x['quality_eligible']=='no' and x['std_axis']=='between_seed_means' for x in rows))
 def test_qc_failure_cannot_enter_seed_aggregation(self):
  from meta_acdc.structure.ensemble_mean import load_seed
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'bad.tsv';p.write_text('job_id\tscore\tqc\nx\t.9\tlow_iptm\n')
   with self.assertRaisesRegex(ValueError,'QC-pass'):load_seed(p)
 def test_grouping_eight_complex_aliases_is_one_group(self):
  from meta_acdc.models.train_domain_adapt import pdb_group
  aliases=['5brz','5brz_0','5brz_1','5brz-disp','5brz<-3uts','af3-5brz_EVDPIGHLY','af3-5brz_EVDPIGHLY-disp','tcr-5brz_EVDPIGHLY']
  self.assertEqual({pdb_group(x) for x in aliases},{'5brz'})
 def test_real_graft_length_and_missing_backbone_rejected(self):
  from meta_acdc.models.dataset import make_decoy,_transplanted_atoms
  from meta_acdc.structure.cif import parse_cif
  res=parse_cif(FIX/'5brz_EVDPIGHLY_model_0.cif');roles={'A':'tcr','B':'tcr','C':'peptide','D':'mhc','E':'b2m'}
  for x in res:x.role=roles[x.chain]
  pep=[x for x in res if x.chain=='C']
  with self.assertRaisesRegex(ValueError,'lengths must match'):make_decoy(res,pep+pep[:1])
  bad=replace(pep[0],atoms=[a for a in pep[0].atoms if a[0]!='N'])
  with self.assertRaisesRegex(ValueError,'complete N/CA/C'):_transplanted_atoms(pep[0],bad)
 def test_figure_reads_registered_report_and_rejects_unbound_input(self):
  from meta_acdc.analysis.verified_figures import load_validation_comparison
  report=ROOT/'docs/affinity_development_20261018/bounded_run/comparison.json';manifest=ROOT/'docs/affinity_development_20261018/delivery_manifest.json'
  values=load_validation_comparison(report,manifest);raw=json.loads(report.read_text());self.assertEqual(dict(values)['anchor_raw'],raw['anchor_validation']['rmse_pkd'])
  with self.assertRaisesRegex(ValueError,'not bound'):load_validation_comparison(D/'dispersion_audit.json',manifest)
 def test_clinical_scan_refuses_old_mode_and_audit_does_not_score(self):
  from meta_acdc.structure.clinical_scan import main
  with self.assertRaisesRegex(ValueError,'Clinical scoring blocked'):main([])
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'receipt.json';main(['--audit-withdrawal','--out',str(p)]);receipt=json.loads(p.read_text());self.assertEqual(receipt['scores_emitted'],0);self.assertEqual(receipt['hazard_flags_emitted'],0)
   with self.assertRaisesRegex(ValueError,'overwrite'):main(['--audit-withdrawal','--out',str(p)])
if __name__=='__main__':unittest.main()
