import sys,tempfile,unittest
from pathlib import Path
import numpy as np
from meta_acdc.models.development_affinity import train,predict,chemical_features,feature_map,data_firewall,ROOT,inner_cv
class Tests(unittest.TestCase):
 def row(self,p='ACADACA'):
  return dict(evidence_id='x',split_block='a',cdr3_alpha='CAVNNAGNMLTF',cdr3_beta='CASSETRGAPYGYTF',alpha_variable='ACDEFGHIKLMNPQRSTVWY'*5,beta_variable='ACDEFGHIKLMNPQRSTVWY'*5,peptide=p,primary_assay={'instrument':'Biacore T200'},mhc_allele='HLA-A*02:01',temperature_K=298.15,pkd=5.)
 def test_word_collision_resolved(self):
  a=self.row();b=self.row('ACACADA');self.assertEqual(feature_map(a,'anchor_raw'),feature_map(b,'anchor_raw'));self.assertNotEqual(chemical_features(a),chemical_features(b))
 def test_metadata_omitted(self):
  for name in ('sequence_balanced','cdr3_chemical','cdr3_words'):
   a=self.row();b=dict(a,mhc_allele='HLA-X*99:99',primary_assay={'instrument':'Different'})
   self.assertEqual(feature_map(a,name),feature_map(b,name))
 def test_firewall_blocks_consumed_and_new_holdout(self):
  with tempfile.TemporaryDirectory() as tmp:
   with data_firewall([],tmp):
    with self.assertRaises(PermissionError):(ROOT/'docs/independent_spr_20261018/frozen_test.json').read_text()
    with self.assertRaises(PermissionError):(ROOT/'docs/affinity_baseline_20261018/acceptance/test_predictions.json').read_text()
 def test_weighted_group_unit_trace(self):
  a=self.row();b=self.row('GGDDDDD');b.update(evidence_id='y',split_block='b',pkd=4.)
  m=train([a,b],'cdr3_chemical',.1);self.assertEqual(m['n_fit'],2);self.assertTrue(all(np.isfinite(m['coefficients'])));self.assertTrue(all(np.isfinite(predict(m,[a,b]))))
 def test_inner_cv_block_separation(self):
  rows=[]
  for i,p in enumerate(['ACADACA','GGDDDDD','YYWWWWW']):
   x=self.row(p);x.update(evidence_id=str(i),split_block=str(i),pkd=4+i);rows.append(x)
  alpha,cv=inner_cv(rows,'cdr3_chemical',[.1,1.]);self.assertIn(alpha,[.1,1.])
  for c in cv:
   for f in c['folds']:self.assertFalse(set(f['fit_ids'])&set(f['held_ids']))
if __name__=='__main__':unittest.main()
