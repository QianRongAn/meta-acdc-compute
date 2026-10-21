"""Publication evidence checks; no test-label consumption or training."""
import importlib.util,json,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/'docs/preprint_20261021'
spec=importlib.util.spec_from_file_location('preprint_evidence',DIR/'build_evidence.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class PreprintEvidenceTests(unittest.TestCase):
    def test_source_hashes_and_no_reserved_inputs(self):
        source=m.load_sources();self.assertEqual(len(source),15)
        self.assertFalse(any('fresh_spr_reserve' in p or 'test_predictions' in p or 'frozen_test.json' in p for p in source))
    def test_all_pairs_and_fixed_qc_calculation(self):
        s,f,p=m.derive(m.load_sources());self.assertEqual(len(p),588)
        a={r['kind']:r for r in s['atom']};self.assertEqual(a['graft']['n_severe_overlap'],146)
        self.assertEqual(a['native']['n_severe_overlap'],0);self.assertEqual(a['displaced']['n_severe_overlap'],300)
        self.assertAlmostEqual(a['graft']['equal_family_overlap_fraction'],0.5798573080199586)
        self.assertAlmostEqual(s['qc_only_discriminator'][0]['balanced_accuracy'],(1+146/288)/2)
    def test_every_candidate_and_mean_comparator_retained(self):
        s,_,_=m.derive(m.load_sources());self.assertEqual(len(s['validation']),5)
        anchor=s['validation'][0]
        self.assertTrue(all(r['rmse_pkd']>anchor['rmse_pkd'] for r in s['validation'][1:4]))
        self.assertAlmostEqual(s['rmse_relative_improvement'],0.0073902603815091,places=5)
    def test_distance_grid_boundary_duplicates_and_negative_cells(self):
        spec=importlib.util.spec_from_file_location('atomic_grid',DIR/'verify_atomic_inputs.py')
        v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
        points=[(0.,0.,0.),(1.2,0.,0.),(1.2001,0.,0.),(-1.2,0.,0.),(0.,0.,0.)]
        grid=v.spatial_grid(points)
        self.assertEqual(v.count_close_pairs([(0.,0.,0.)],grid),4)
        self.assertEqual(v.count_close_pairs([(20.,0.,0.)],grid),0)
        with self.assertRaises(ValueError):v.spatial_grid([(float('nan'),0,0)])
        with self.assertRaises(ValueError):v.spatial_grid(points,0)
    def test_later_context_has_separate_endpoint_and_no_new_acceptance(self):
        import hashlib
        manifest=json.loads((DIR/'context_manifest.json').read_text())
        for name,digest in manifest['files'].items():self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(),digest)
        r=json.loads((ROOT/'docs/six_blockers_20261021/ergo_comparison/report.json').read_text())
        self.assertFalse(r['new_independent_acceptance']);self.assertFalse(r['model_fit_or_tuning_performed'])
        self.assertFalse(r['physical_binding_validation_passed']);self.assertEqual(r['n_observations'],10000)
    def test_export_requires_historical_helpers_not_new_main_generator(self):
        spec=importlib.util.spec_from_file_location('historical_export',DIR/'export_atomic_inputs.py')
        v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
        self.assertTrue(v.verify_snapshot_helpers())
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for name in json.loads((DIR/'atomic_bundle_v1/export_receipt.json').read_text())['code_files']:
                p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'changed helper')
            with self.assertRaises(ValueError):v.verify_snapshot_helpers(root)
    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):m.build(Path(d),figures=False)
    def test_references_and_excluded_claims(self):
        refs=json.loads((DIR/'references.json').read_text());self.assertEqual(len(refs),13)
        self.assertEqual(len({r['doi'].lower() for r in refs}),13)
        claims=json.loads((DIR/'claim_evidence.json').read_text())
        self.assertTrue(all(not r['supported'] for r in claims if r['id'].startswith('X')))
        readiness=json.loads((DIR/'submission_readiness.json').read_text());self.assertFalse(readiness['submission_ready'])
        self.assertFalse(readiness['validated_binding_or_clinical_model'])
if __name__=='__main__':unittest.main()
