"""Regression checks for literal primary roles; no identity/atom repair."""
import unittest
from pathlib import Path
from meta_acdc.structure.primary_role_evidence import classify_role, primary_roles
ROOT=Path(__file__).resolve().parents[1]
class PrimaryRoleEvidenceTests(unittest.TestCase):
 def test_literal_roles(self):
  for text,expected in [('T-CELL RECEPTOR ALFA CHAIN','alpha'),('ALFALFA','unspecified'),('Valpha7.2','alpha'),('Valpha1','alpha'),('alpha24','alpha'),('Vbeta16','beta'),('beta17','beta'),('T-CELL RECEPTOR BETA-2 CHAIN C REGION','beta'),('BETA-2-MICROGLOBULIN','b2m'),('HLA CLASS I ALPHA CHAIN','mhc'),('AS8.4A','unspecified'),('AS8.4B','unspecified')]:
   with self.subTest(text=text): self.assertEqual(classify_role(text),expected)
 def test_conflicting_primary_fields_are_not_resolved_by_preference(self):
  self.assertEqual(classify_role('T-CELL RECEPTOR ALPHA CHAIN','T-CELL RECEPTOR BETA CHAIN'),'conflicting_alpha_beta')
 def header(self,pdb):
  for folder in ['fingerprint_cross_source_extension_20261024','fingerprint_cross_source_20261024']:
   p=ROOT/'docs'/folder/'sources/headers'/(pdb+'.pdb')
   if p.exists():return p
  self.fail('missing deposited header: '+pdb)
 def test_real_alfa_headers(self):
  for pdb in ['6vrm','6vrn']:
   roles=primary_roles(self.header(pdb));self.assertIn('alpha',[x['role'] for x in roles.values()]);self.assertIn('beta',[x['role'] for x in roles.values()])
 def test_actual_conflicts_remain_unresolved(self):
  for pdb in ['5e6i','6at6','7n5c','7n5p']:
   self.assertIn('conflicting_alpha_beta',[x['role'] for x in primary_roles(self.header(pdb)).values()],pdb)
 def test_real_beta2_is_not_b2m(self):
  roles=primary_roles(self.header('5e6i'));matched=[x for x in roles.values() if 'BETA-2 CHAIN C REGION' in x['molecule'].upper()];self.assertTrue(matched);self.assertTrue(all(x['role']=='beta' for x in matched))
 def test_ambiguous_AS_names_remain_unresolved(self):
  roles=primary_roles(self.header('8cx4'));self.assertTrue(any(x['role']=='unspecified' and 'AS8.4' in x['molecule'].upper() for x in roles.values()))
if __name__=='__main__': unittest.main()
