"""Regression tests for real-format ingestion and data-release failure modes."""
import csv,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import torch
from meta_acdc.data.ingest_iedb import ingest,ic50_score,assay_family
from meta_acdc.data.release import load_release,sha256,complex_id
from meta_acdc.data.build_release import assign_folds
from meta_acdc.models.benchmark_baseline import load
from meta_acdc.models.dataset import make_decoy,_transplanted_atoms
from meta_acdc.structure.graph import Residue
from meta_acdc.structure.alanine_scan import to_mutant
from meta_acdc.structure.quality import confidence,qc_file

H1=['Assay ID','Assay','Epitope','MHC Restriction','Assay','Assay','Assay','Host','Reference']
H2=['IEDB IRI','Method','Name','Name','Qualitative Measurement','Quantitative measurement','Units','Name','PMID']

def source(path,sep,rows):
    with path.open('w',newline='') as f:
        w=csv.writer(f,delimiter=sep);w.writerow(H1);w.writerow(H2);w.writerows(rows)

def row(qual='Positive',value='5',method='binding IC50',unit='nM'):
    return ['assay:1',method,'PEPTIDE','HLA-A*01:01',qual,value,unit,'Homo sapiens','123']

class Ingestion(unittest.TestCase):
    def test_csv_and_tsv_preserve_assays_and_conflicts(self):
        for sep in [',','\t']:
            with self.subTest(sep=sep),tempfile.TemporaryDirectory() as td:
                p=Path(td);source(p/'input',sep,[row(),row('Negative')]);ingest(p/'input',p/'out')
                with (p/'out').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
                self.assertEqual({r['label'] for r in rows},{'0','1'});self.assertEqual(len(rows),2)
                self.assertTrue(all(r['assay_id']=='assay:1' and r['label_scope']=='mhc_assay' for r in rows))
    def test_bad_input_does_not_truncate_output(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);(p/'in').write_text('bad\nheader\n');(p/'out').write_text('preserved')
            with self.assertRaises(ValueError):ingest(p/'in',p/'out')
            self.assertEqual((p/'out').read_text(),'preserved')
    def test_max_rows_is_enforced(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);source(p/'in',',',[row()]*5);self.assertEqual(ingest(p/'in',p/'out',max_rows=2).n_records,2)
    def test_endpoint_and_unit_boundaries(self):
        self.assertEqual(assay_family('binding KD'),'binding-kd');self.assertEqual(assay_family('binding koff'),'kinetics_or_stability')
        self.assertEqual(assay_family('purified MHC binding'),'binding-unspecified')
        self.assertAlmostEqual(ic50_score('5','uM')[0],-np.log10(5000))
        for value,unit in [('0','nM'),('-1','nM'),('nan','nM'),('inf','nM'),('<5','nM'),('5','nmol'),('5','')]:self.assertIsNone(ic50_score(value,unit)[0])
    def test_old_vdj_labels_are_refused(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'old.tsv';p.write_text('source\tcdr3_beta\tpeptide\tlabel\nvdjdb\tCASS\tPEPTIDE\t0\n')
            with self.assertRaisesRegex(ValueError,'Unverified'):load(p)

class StructuralSafety(unittest.TestCase):
    def test_prefix_graft_is_refused(self):
        rs=[Residue('P','ALA',i,0,0,0,atoms=[('N',0,0,0),('CA',1,0,0),('C',2,0,0)]) for i in range(9)]
        for r in rs:r.role='peptide'
        with self.assertRaisesRegex(ValueError,'lengths'):make_decoy(rs,rs[:8])
    def test_missing_backbone_graft_refused(self):
        a=Residue('P','ALA',1,0,0,0,atoms=[('CA',0,0,0)])
        with self.assertRaises(ValueError):_transplanted_atoms(a,a)
    def test_gly_to_ala_cannot_fake_a_sidechain(self):
        r=Residue('P','GLY',1,0,0,0,atoms=[('CA',0,0,0)])
        with self.assertRaises(ValueError):to_mutant(r,'ALA')
    def test_split_identity_aliases(self):
        self.assertEqual({complex_id(x) for x in ['5brz','5brz_0','af3-5brz_EVDPIGHLY','5brz_1<-1ao7']},{'5brz'})
    def test_family_components_and_fold_stability(self):
        rng=np.random.default_rng(14);aa=np.array(list('ACDEFGHIKLMNPQRSTVWY'))
        records=[{'pdb':f'p{i:03d}','complex_id':f'p{i:03d}','alpha':''.join(rng.choice(aa,100)),'beta':''.join(rng.choice(aa,100))} for i in range(20)]
        records.append(dict(records[0],complex_id='p000_copy'))
        copy=json.loads(json.dumps(records));assign_folds(records);assign_folds(copy)
        self.assertEqual([r['fold'] for r in records],[r['fold'] for r in copy]);self.assertEqual(records[0]['fold'],records[-1]['fold'])
        self.assertEqual({r['fold'] for r in records},{'train','validation','calibration','test'})
    def test_per_model_confidence_not_borrowed(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);(p/'job_conf.json').write_text('{"iptm":0.99}')
            self.assertEqual(confidence(p/'job_model_0.cif')[0],.99)
            self.assertIsNone(confidence(p/'job_model_1.cif')[0])
    def test_high_confidence_parse_failure_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);f=p/'1abc_PEPTIDE_model_0.cif';f.write_text('bad');(p/'1abc_PEPTIDE_conf.json').write_text('{"iptm":0.99}')
            q=qc_file(f,{})
            self.assertIn('structure_or_identity_failed',q['flag'])

class ManifestSafety(unittest.TestCase):
    def test_manifest_required(self):
        with self.assertRaises(ValueError):load_release(None)
    def test_artifact_tamper_refused(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);(p/'graphs.pt').write_bytes(b'data');code=Path(__file__).resolve().parents[1];name='src/meta_acdc/data/release.py'
            m={'schema_version':1,'status':'validated','task':'native_vs_synthetic_decoy','training_allowed':True,'code_files':{name:sha256(code/name)},'files':{'graphs.pt':{'sha256':sha256(p/'graphs.pt')}},'dataset':'graphs.pt'}
            (p/'release.json').write_text(json.dumps(m));load_release(p/'release.json')
            (p/'graphs.pt').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):load_release(p/'release.json')
    def test_binding_scope_is_not_authorized(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'manifest';p.write_text(json.dumps({'schema_version':1,'status':'validated','task':'native_vs_synthetic_decoy','training_allowed':True}))
            with self.assertRaises(ValueError):load_release(p,task='tcr_binding')

if __name__=='__main__':unittest.main()


class AtomIntegrity(unittest.TestCase):
    def test_missing_atom_is_not_graft_eligible(self):
        from meta_acdc.structure.graph import Residue
        from meta_acdc.structure.identity import validate_peptide_atoms
        r=Residue('P','ALA',1,0,0,0,role='peptide',atoms=[(a,0,0,0) for a in ['N','CA','C','O']])
        with self.assertRaisesRegex(ValueError,'incomplete'):validate_peptide_atoms([r])
        r.atoms.append(('CB',1,0,0));validate_peptide_atoms([r])


class ModelConsistency(unittest.TestCase):
    def test_single_and_mixed_batch_match_scoring(self):
        import numpy as np,torch
        from meta_acdc.models.egnn import EGNN
        from meta_acdc.models.train_egnn import collate
        from meta_acdc.structure.graph import InterfaceGraph
        from meta_acdc.structure.score_predictions import score_graph
        torch.manual_seed(3)
        rng=np.random.default_rng(7)
        def graph(n):
            h=rng.normal(size=(n,26));h[:,24]=0;h[0,24]=1
            x=rng.normal(size=(n,3));ei=[(i,j) for i in range(n) for j in range(n) if i!=j]
            ef=rng.normal(size=(len(ei),16));return h.tolist(),x.tolist(),ei,ef.tolist()
        a,b=graph(8),graph(13);net=EGNN(26,16,depth=2,hidden=16).eval()
        with torch.no_grad():
            solo=net(**collate([a]))[0];mixed=net(**collate([a,b]))[0]
        self.assertTrue(torch.allclose(solo,mixed,atol=1e-6,rtol=1e-5))
        g=InterfaceGraph(n_nodes=8,node_features=a[0],node_coords=a[1],edge_index=a[2],edge_features=a[3],peptide_mask=[True]+[False]*7)
        self.assertAlmostEqual(score_graph(net,g,'cpu'),torch.sigmoid(solo).item(),places=6)
        moved=(a[0],(np.array(a[1])+np.array([20,-7,4])).tolist(),a[2],a[3])
        with torch.no_grad():translated=net(**collate([moved]))[0]
        self.assertTrue(torch.allclose(solo,translated,atol=1e-6,rtol=1e-5))
