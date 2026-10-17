"""New audit counterexamples; synthetic fixtures never constitute biological validation."""
import copy, csv, json, tempfile, unittest
from pathlib import Path
import numpy as np
import torch
from meta_acdc.structure.graph import (Residue, InterfaceGraph, TCR, PEPTIDE, MHC, OTHER,
    classify_chains, is_tcr_chain, build_interface_graph_from_residues, _node_features)
from meta_acdc.structure.feature_schema import AA_ORDER, PEPTIDE_FLAG_INDEX, validate_node_features
from meta_acdc.structure.identity import missing_heavy_atoms
from meta_acdc.structure.artifact_qc import construction_screen
from meta_acdc.structure.da_stability import analyse, load_jobs
from meta_acdc.models.train_hybrid import collate
from meta_acdc.models.hybrid import AA_IDX, HybridModel
from meta_acdc.active_learning.simulate_al import fixed_budget_batches
from meta_acdc.structure.vdjdb_structure_map import association_hits


def residue(chain, n, role, offset=0):
    atoms=[('N',offset,0,0),('CA',offset+1,0,0),('C',offset+2,0,0),('O',offset+2,1,0),('CB',offset+1,1.5,0)]
    return Residue(chain,'ALA',n,offset+1,0,0,atoms=atoms,role=role)


def graph_residues():
    return [residue(c,i,role) for c,role in [('T',TCR),('H',MHC),('P',PEPTIDE)] for i in range(1,6)]


class IdentityAndAtoms(unittest.TestCase):
    def test_motif_length_and_letters_do_not_assign_roles(self):
        residues=[Residue('A','TYR',i,0,0,0) for i in range(100)]
        residues += [Residue('B',aa,i,0,0,0) for i,aa in enumerate(['TYR','PHE','CYS'])]
        self.assertEqual(set(classify_chains(residues).values()),{OTHER})
        with self.assertRaises(ValueError):is_tcr_chain('A')
    def test_summary_uses_roles_not_letters(self):
        g=InterfaceGraph(n_nodes=3,node_roles=[MHC,TCR,PEPTIDE],residue_info=[('A','ALA',1),('X','ALA',1),('B','ALA',1)])
        self.assertIn('TCR=1, pMHC=2',g.summary())
    def test_missing_receptor_atoms_are_rejected(self):
        rs=graph_residues();rs[0].atoms=[a for a in rs[0].atoms if a[0]!='CB']
        with self.assertRaisesRegex(ValueError,'incomplete interface'):build_interface_graph_from_residues(rs)
        diagnostic=build_interface_graph_from_residues(rs,require_complete_atoms=False)
        self.assertTrue(diagnostic.missing_atom_nodes)
    def test_complete_gly_is_not_missing_sidechain(self):
        r=residue('P',1,PEPTIDE);r.resname='GLY';r.atoms=[a for a in r.atoms if a[0]!='CB']
        self.assertEqual(missing_heavy_atoms(r),[])
    def test_wrong_feature_width_rejected(self):
        with self.assertRaises(ValueError):validate_node_features([[0]*27])
    def test_file_order_does_not_change_peptide_order(self):
        a=build_interface_graph_from_residues(graph_residues())
        b=build_interface_graph_from_residues(list(reversed(graph_residues())))
        self.assertEqual(a.node_keys,b.node_keys)
    def test_duplicate_residue_identity_rejected(self):
        rs=graph_residues()
        with self.assertRaisesRegex(ValueError,'duplicate'):build_interface_graph_from_residues(rs+[copy.deepcopy(rs[0])])
    def test_overlap_screen_distinguishes_controls_without_nonbinding_claim(self):
        a=residue('P',1,PEPTIDE);b=residue('T',1,TCR,offset=20)
        self.assertTrue(construction_screen([a,b])['passed'])
        b=residue('T',1,TCR)
        report=construction_screen([a,b]);self.assertFalse(report['passed'])
        self.assertFalse(report['validates_nonbinding']);self.assertFalse(report['validates_rotamers'])
    def test_bad_graft_ca_cb_bond_is_rejected(self):
        a=residue('P',1,PEPTIDE);a.atoms[-1]=('CB',10,10,10)
        report=construction_screen([a,residue('T',1,TCR,offset=30)])
        self.assertTrue(report['ca_cb_outliers'])


class HybridAndInvariance(unittest.TestCase):
    def test_hybrid_sequence_reads_peptide_and_correct_aa_order(self):
        rs=graph_residues();g=build_interface_graph_from_residues(rs)
        raw=(g.node_features,g.node_coords,g.edge_index,g.edge_features)
        b=collate([raw]);self.assertEqual(int(b['peptide_mask'].sum()),5)
        self.assertEqual(b['seq'][0,:5,AA_IDX['A']].tolist(),[1]*5)
        self.assertEqual(float(b['seq'].sum()),5)
        torch.manual_seed(2);net=HybridModel(hidden=16,egnn_depth=2).eval()
        with torch.no_grad():a=net(**b)[0];mixed=net(**collate([raw,raw]))[0]
        self.assertTrue(torch.allclose(a,mixed,atol=1e-6))
    def test_fixed_coordinate_baseline_preserves_coordinates_and_scalar_invariance(self):
        from meta_acdc.models.egnn import EGNN
        from meta_acdc.models.train_egnn import collate as graph_collate
        g=build_interface_graph_from_residues(graph_residues())
        raw=(g.node_features,g.node_coords,g.edge_index,g.edge_features)
        b=graph_collate([raw]);net=EGNN(26,16,depth=2,hidden=16,update_coordinates=False).eval()
        with torch.no_grad():
            a=net(**b);rotation=torch.tensor([[0.,-1,0],[1,0,0],[0,0,1]])
            b['x']=b['x']@rotation+torch.tensor([5.,6,7]);other=net(**b)
        self.assertTrue(torch.allclose(a,other,atol=1e-6))


class AllSeedAccounting(unittest.TestCase):
    def seeds(self):return [{str(i):float(i) for i in range(6)},{str(i):float(i+1) for i in range(6)}]
    def families(self):return {str(i):'F'+str(i//2) for i in range(6)}
    def test_degenerate_seed_blocks_even_if_survivors_agree(self):
        seeds=self.seeds()+[{str(i):0. for i in range(6)}]
        report=analyse(seeds,self.families(),expected_seeds=3,resamples=10)
        self.assertEqual(report['verdict'],'DEGENERATE_SEED');self.assertIsNone(report['mean_r'])
        self.assertAlmostEqual(report['conditional_nonconstant_mean_r'],1)
        self.assertEqual(len(report['pairs']),3)
    def test_missing_preplanned_seed_rejected(self):
        with self.assertRaises(ValueError):analyse(self.seeds(),self.families(),expected_seeds=3)
    def test_unequal_job_coverage_cannot_pass(self):
        seeds=self.seeds();seeds[1]['extra']=9
        report=analyse(seeds,self.families(),expected_seeds=2,resamples=10)
        self.assertEqual(report['verdict'],'INCOMPLETE_COVERAGE')
    def test_undefined_bootstraps_are_counted(self):
        seeds=self.seeds();families={str(i):'F'+str(i//3) for i in range(6)}
        for s in seeds:
            for i in range(6):s[str(i)]=float(i//3)
        report=analyse(seeds,families,expected_seeds=2,resamples=100)
        self.assertGreater(report['bootstrap']['undefined'],0)
        self.assertEqual(report['bootstrap']['defined']+report['bootstrap']['undefined'],100)
    def test_duplicate_and_nonfinite_scores_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'x.tsv'
            for rows in ['a\t0.1\na\t0.2\n','a\tnan\n']:
                p.write_text('job_id\tscore\n'+rows)
                with self.assertRaises(ValueError):load_jobs(p)


class BudgetAndEvidence(unittest.TestCase):
    def test_fixed_budget_uses_short_final_batch(self):
        for b in [1,3,7,30]:
            sizes=fixed_budget_batches(100,b,20)
            self.assertEqual(sum(sizes),20);self.assertTrue(all(0<n<=b for n in sizes))
        self.assertEqual(fixed_budget_batches(100,7,20),[7,7,6])
    def test_invalid_budget_fails(self):
        for pool,batch,budget in [(10,0,5),(10,3,0),(10,3,11)]:
            with self.assertRaises(ValueError):fixed_budget_batches(pool,batch,budget)
    def test_zero_confidence_is_association_not_nonbinder_and_chain_specific(self):
        records=[dict(pdb='xxxx',complex_id='xxxx',alpha='AACASSABCDEFG',beta='DIFFERENT',alpha_chain='X',beta_chain='Y',peptide='PEPTIDE')]
        rows=[dict(cdr3_alpha='CASSABCD',cdr3_beta='',peptide='FOREIGN',confidence='0')]
        hits=list(association_hits(records,rows));self.assertEqual(len(hits),1)
        self.assertEqual(hits[0]['confidence'],'0');self.assertEqual(hits[0]['gene'],'TRA')
        rows[0]['cdr3_beta']=rows[0].pop('cdr3_alpha')
        self.assertEqual(list(association_hits(records,rows)),[])
    def test_iedb_unit_white_list(self):
        from meta_acdc.data.ingest_iedb import ic50_score
        for unit in ['nM','uM','µM','mM','M']:
            # Test only supported spellings; unsupported concentration units remain raw evidence.
            score,_=ic50_score('5',unit)
            if unit in ['nM','uM']:self.assertIsNotNone(score)
        for unit in ['mg/mL','nM-ish','nmol','', 'unknown']:
            self.assertIsNone(ic50_score('5',unit)[0])
    def test_filtering_inherits_frozen_folds_without_reshuffling(self):
        from meta_acdc.data.build_release import inherit_folds
        before=[dict(complex_id=str(i),pdb=str(i),alpha='A',beta='B',peptide='P',domain='crystal',family_id='F'+str(i),fold='test' if i==1 else 'train') for i in range(3)]
        after=[{k:v for k,v in before[1].items() if k not in {'fold','family_id'}}]
        kept,excluded=inherit_folds(after,before)
        self.assertEqual(after[0]['fold'],'test');self.assertEqual(after[0]['family_id'],'F1')
        self.assertEqual(kept,{'1'});self.assertEqual(excluded,[])
        after[0]['alpha']='CHANGED'
        with self.assertRaises(ValueError):inherit_folds(after,before)
    def test_similar_cdr3_transitive_clusters_and_unseen_both(self):
        from meta_acdc.data.splits import cluster_sequences, group_split
        seqs=['CASSLGQETQYF','CASSLGQETQYA','CASSLGQETQYG','CAWRGGANQF','CAVRDSNYQLI','CSARDTGGYF']
        clusters=cluster_sequences(seqs);self.assertEqual(len(set(clusters[:3])),1)
        peptides=['AAAAAAAAA','CCCCCCCCC','DDDDDDDDD','EEEEEEEEE','FFFFFFFFF','GGGGGGGGG']
        split=group_split(seqs,peptides,group_by='unseen-both',test_frac=.4,seed=4)
        train,test=split.train_idx,split.test_idx
        self.assertTrue(train and test)
        self.assertFalse({clusters[i] for i in train}&{clusters[i] for i in test})
        self.assertFalse({peptides[i] for i in train}&{peptides[i] for i in test})


class ConfidenceProvenance(unittest.TestCase):
    def cif(self,path,value):
        path.write_text('data_fixture\nloop_\n'+'\n'.join([
            '_atom_site.group_PDB','_atom_site.label_atom_id','_atom_site.auth_asym_id',
            '_atom_site.auth_comp_id','_atom_site.auth_seq_id','_atom_site.Cartn_x',
            '_atom_site.Cartn_y','_atom_site.Cartn_z','_atom_site.B_iso_or_equiv'])+
            '\nATOM CA A ALA 1 0 0 0 '+str(value)+'\n#\n')
    def test_unspecified_cif_does_not_convert_b_factor(self):
        from meta_acdc.structure.cif import parse_cif
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'x.cif';self.cif(p,85)
            self.assertIsNone(parse_cif(p)[0].plddt)
            self.assertEqual(parse_cif(p,confidence_kind='af3_plddt')[0].plddt,85)
    def test_zero_plddt_is_not_missing(self):
        from meta_acdc.structure.cif import parse_cif
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'x.cif';self.cif(p,0)
            self.assertEqual(parse_cif(p,confidence_kind='af3_plddt')[0].plddt,0)
    def test_out_of_range_plddt_is_not_silently_clipped(self):
        from meta_acdc.structure.cif import parse_cif
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'x.cif';self.cif(p,110)
            with self.assertRaises(ValueError):parse_cif(p,confidence_kind='af3_plddt')


if __name__=='__main__':unittest.main()
