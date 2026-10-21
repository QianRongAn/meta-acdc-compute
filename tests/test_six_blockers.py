"""Mandatory source-backed data and benchmark checks; no optional skips."""
import csv,hashlib,importlib.util,io,json,math,pickle,re,statistics,unittest,zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from sklearn.metrics import roc_auc_score,average_precision_score
from meta_acdc.data.clinical_primary import CODONS
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'docs/six_blockers_20261021';S=D/'sources'
def load(name):return json.loads((D/name).read_text())
def cells(path,sheet):
 ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
 with zipfile.ZipFile(path) as z:
  book=ET.fromstring(z.read('xl/workbook.xml'))
  rid=next(x for x in book.findall('.//s:sheet',ns) if x.get('name')==sheet).get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
  target=next(x for x in ET.fromstring(z.read('xl/_rels/workbook.xml.rels')) if x.get('Id')==rid).get('Target')
  target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
  strings=[]
  if 'xl/sharedStrings.xml' in z.namelist():strings=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
  result={}
  for c in ET.fromstring(z.read(target)).findall('.//s:sheetData/s:row/s:c',ns):
   v=c.findtext('s:v',default='',namespaces=ns)
   if c.get('t')=='s':value=strings[int(v)]
   elif c.get('t')=='inlineStr':value=''.join(c.find('s:is',ns).itertext())
   else:
    try:value=float(v)
    except ValueError:value=v
   result[c.get('r')]=value
  return result
class SixBlockerTests(unittest.TestCase):
 def test_mandatory_original_source_bytes(self):
  for name,meta in load('required_sources.json').items():
   p=ROOT/meta['path'];self.assertTrue(p.is_file(),name)
   self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),meta['sha256'],name)
 def test_finite_checkpoint_search_cannot_restore_21_row_claim(self):
  plan=load('legacy_21_search_plan.json');result=load('legacy_21_search_results.json')
  self.assertEqual(len(plan['eligible_candidates']),8);self.assertEqual(len(result['input_sha256']),105)
  self.assertEqual(len({x['sha256'] for x in plan['eligible_candidates']}),8)
  self.assertEqual(result['complete_matches'],[])
  self.assertTrue(all(x['comparisons']==105 and x['graph_counts_match']==105 and x['matches']<=1 for x in result['results']))
  self.assertFalse(result['scientific_model_promotion_allowed'])
 def test_primary_response_and_temperature_match_real_source(self):
  d=load('tcr4_primary_challenge.json');c=cells(S/'pi3kca_source_data.xlsx','Fig 1c')
  self.assertEqual((c['A3'],c['B3']),(208.24,-9.68758))
  self.assertEqual(d['measured_source_curves']['WT'][0],dict(concentration_uM=c['A3'],response_RU=c['B3']))
  self.assertEqual(d['measured_source_curves']['neo'][0]['concentration_uM'],c['F3'])
  article=ET.fromstring((S/'PMC11756396.xml').read_bytes())
  captions=[''.join(f.itertext()) for f in article.iter('fig') if f.get('id')=='Fig1']
  self.assertTrue(any('4' in x and '°C' in x and 'six' in x for x in captions))
  self.assertTrue(all(r['temperature_K']==277.15 for r in d['observations']))
 def test_non_detection_has_no_invented_affinity_or_safety_label(self):
  d=load('tcr4_primary_challenge.json');negative=next(x for x in d['observations'] if not x['assay_binding_detected'])
  for k in ['reported_KD_uM','valid_KD_lower_bound_uM','quantified_detection_limit_RU','universal_binding_label','clinical_risk_label']:self.assertIsNone(negative[k])
  self.assertFalse(negative['point_regression_target_available']);self.assertFalse(negative['training_allowed'])
  self.assertEqual(negative['assayed_analyte_max_uM'],208.24)
 def test_name_collision_and_consumed_family_not_relabelled_independent(self):
  d=load('tcr4_exposure_audit.json');self.assertIn('TCR4',d['registry_hits']['consumed_frozen_test_identity_only'])
  self.assertEqual(d['new_receptor_lineages'],0);self.assertEqual(d['new_fresh_independent_acceptance_records'],0)
  candidates=load('batcave_candidate_index.json')['constructs']
  other=next(x for x in candidates if x['pmid']=='35589842' and x['tcr_name']=='TCR4')
  self.assertNotEqual(other['cdr3_alpha'],d['identity_fields']['cdr3_alpha'])
 def test_batcave_assays_species_and_censoring_stay_separate(self):
  d=load('batcave_candidate_index.json');kd=d['affinity_candidate_summary']
  self.assertEqual((kd['rows'],kd['censored_rows'],kd['human_censored_rows']),(560,67,0))
  self.assertTrue(kd['all_censored_rows_nonhuman']);self.assertFalse(kd['training_allowed'])
  self.assertEqual(d['counts']['batcave_mhci.xlsx']['human_canonical_paired_sequence_rows'],5798)
  for r in d['constructs']:
   self.assertIsNone(r['physical_binding_label']);self.assertIsNone(r['clinical_risk_label']);self.assertFalse(r['training_allowed'])
 def test_primary_research_DNA_exact_cells_and_author_exclusion(self):
  d=load('pteam_primary_sequence_index.json');by={x['name']:x for x in d['receptors']}
  self.assertEqual(len(by),27);self.assertEqual(sum(x['DNA_pair_passes'] for x in by.values()),25)
  cmv=cells(S/'pteam_cmv.xlsx','TCR sequences');neo=cells(S/'pteam_rnf43.xlsx','Sequences')
  for r in by.values():
   table=cmv if r['source_file']=='pteam_cmv.xlsx' else neo;i=r['source_row']
   for chain,col in [('alpha','K'),('beta','L' if table is cmv else 'M')]:
    dna=re.sub(r'\s+','',str(table[col+str(i)])).upper();self.assertEqual(dna,r[chain]['dna'])
    aa=''.join(CODONS[dna[j:j+3]] for j in range(0,len(dna)-2,3)).removesuffix('*')
    self.assertEqual(aa,r[chain]['amino_acids'])
   self.assertIsNone(r['clinical_administered_product_link']);self.assertFalse(r['training_allowed']);self.assertFalse(r['model_scored'])
  self.assertFalse(by['TCR1-4']['beta']['published_translation_matches']);self.assertFalse(by['TCR1-4']['DNA_pair_passes'])
  self.assertEqual(by['R27']['alpha']['trailing_nt'],'A');self.assertTrue(by['R27']['author_excluded'])
  self.assertEqual(sum(x['eligible_for_further_primary_review'] for x in d['summary'].values()),4072)
 def test_ergo_keeps_all_pairs_control_and_fixed_baseline(self):
  report=load('ergo_comparison/report.json');rows=load('ergo_comparison/predictions.json')
  with (D/'ergo_comparison/queries.csv').open() as f:q=list(csv.reader(f))
  with (D/'ergo_comparison/official_stdout.tsv').open() as f:raw=list(csv.reader(f,delimiter='\t'))
  self.assertEqual((len(rows),len(q),len(raw)),(10000,10001,10001));self.assertEqual(len({x['evidence_id'] for x in rows}),10000)
  self.assertEqual([x[:2] for x in raw],q);self.assertEqual(q[0],q[-1]);self.assertAlmostEqual(float(raw[0][2]),float(raw[-1][2]),places=5)
  for row,values in zip(rows,raw):self.assertEqual(row['ergo_score'],float(values[2]));self.assertEqual(row['peptide'],values[1])
  for col in ['ergo_score','spr_reference_pkd_score']:
   auc=[];ap=[]
   for pep in sorted({x['peptide'] for x in rows}):
    group=[x for x in rows if x['peptide']==pep];y=[x['observed_activation'] for x in group];v=[x[col] for x in group]
    auc.append(roc_auc_score(y,v));ap.append(average_precision_score(y,v))
   self.assertAlmostEqual(statistics.mean(auc),report['metrics'][col]['macro_auroc'],places=12)
   self.assertAlmostEqual(statistics.mean(ap),report['metrics'][col]['macro_average_precision'],places=12)
  self.assertFalse(report['new_independent_acceptance']);self.assertFalse(report['clinical_validation_passed'])
  self.assertFalse(report['model_fit_or_tuning_performed']);self.assertTrue(report['control_excluded_from_metrics'])
 def test_official_training_exposure_audit_and_blob_identity(self):
  class DataOnly(pickle.Unpickler):
   def find_class(self,module,name):raise ValueError('Executable pickle')
  p=S/'ergo_lstm_mcpas_train1.pickle';data=p.read_bytes()
  self.assertEqual(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest(),'d6864bd6330cb462584d9b809af1fb58d6d7b1f1')
  training=DataOnly(io.BytesIO(data)).load();pairs={(x[0],x[1][0] if isinstance(x[1],(list,tuple)) else x[1]) for x in training};betas={a for a,b in pairs}
  with (D/'ergo_comparison/queries.csv').open() as f:queries=list(csv.reader(f))[:10000]
  self.assertEqual(sum(x[0] in betas for x in queries),120);self.assertEqual(sum(tuple(x) in pairs for x in queries),0)
  self.assertFalse(load('tcrmatch_feasibility.json')['labels_inferred_from_no_hit'])
 def test_clinical_and_historical_archive_limits_remain_explicit(self):
  self.assertFalse(load('clinical_gap_decisions.json')['clinical_model_released'])
  self.assertFalse(load('archive_status.json')['original_archive_byte_identity_verified'])
if __name__=='__main__':unittest.main()
