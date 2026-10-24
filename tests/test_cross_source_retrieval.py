import unittest,json,copy
from pathlib import Path
from meta_acdc.structure.source_balanced_retrieval import source_balanced_blocks,representatives,joint_components,author_keys
ROOT=Path(__file__).resolve().parents[1]
def row(pdb,rid,study,author='A',context='V1'):
    return dict(pdb=pdb,receptor_id=rid,source_study_key=study,alpha_chain='D',beta_chain='E',V_context=context,primary_citation={'rcsb_authors':[author]},alpha={'variable':rid+'A'},beta={'variable':rid+'B'})
class CrossSourceChecks(unittest.TestCase):
    def test_different_pdb_same_paper_not_independent(self):
        self.assertFalse(source_balanced_blocks([row('1aaa','R1','S1'),row('2aaa','R1','S1'),row('3aaa','R2','S1')]))
    def test_positive_negative_different_gallery_studies_rejected(self):
        self.assertFalse(source_balanced_blocks([row('1aaa','R1','S1'),row('2aaa','R1','S2'),row('3aaa','R2','S3')]))
    def test_source_balanced_and_author_disjoint(self):
        rows=[row('1aaa','R1','S1','A'),row('2aaa','R1','S2','B'),row('3aaa','R2','S2','B')]
        blocks=source_balanced_blocks(rows,True);self.assertEqual(len(blocks),1);self.assertEqual(blocks[0]['positive_pdbs'],['2aaa']);self.assertEqual(blocks[0]['negative_pdbs'],['3aaa'])
    def test_shared_author_not_independent_team(self):
        rows=[row('1aaa','R1','S1','Same, A.'),row('2aaa','R1','S2','same A'),row('3aaa','R2','S2','same A')]
        self.assertEqual(len(source_balanced_blocks(rows)),1);self.assertFalse(source_balanced_blocks(rows,True))
    def test_author_diacritics_not_false_independence(self):
        self.assertEqual(author_keys(row('1aaa','R','S','García, K.')),author_keys(row('2aaa','R','S','Garcia K')))
    def test_unknown_author_not_independent(self):
        rows=[row('1aaa','R1','S1','?'),row('2aaa','R1','S2','A'),row('3aaa','R2','S2','A')]
        self.assertFalse(source_balanced_blocks(rows,True))
    def test_gene_context_mismatch_excluded(self):
        rows=[row('1aaa','R1','S1'),row('2aaa','R1','S2'),row('3aaa','R2','S2',context='V2')]
        self.assertFalse(source_balanced_blocks(rows))
    def test_duplicate_copies_not_new_source_units(self):
        a=row('1aaa','R1','S1');b=copy.deepcopy(a);b['alpha_chain']='I';b['beta_chain']='J'
        self.assertEqual(len(representatives([b,a])),1);self.assertEqual(representatives([b,a])[0]['alpha_chain'],'D')
    def test_shared_chain_joint_exclusion(self):
        a=row('1aaa','R1','S1',context='V1');b=row('2aaa','R2','S2',context='V2');b['alpha']['variable']=a['alpha']['variable']
        self.assertEqual(len(joint_components([a,b])),1)
    def test_actual_search_gate_and_missing_atoms(self):
        base=ROOT/'docs/fingerprint_cross_source_20261024';d=json.loads((base/'source_balanced_design.json').read_text());self.assertEqual(d['summary']['n_source_balanced_query_blocks'],0);self.assertFalse(d['training_release_created'])
        qc=json.loads((base/'native_atomic_qc.json').read_text());r=next(x for x in qc['rows'] if x['pdb']=='7rtr' and x.get('reasons'))
        self.assertIn('incomplete interface heavy atoms',r['reasons'][0]);self.assertIn('OE1',r['reasons'][0])
    def test_exact_entity_annotation_inheritance(self):
        data=json.loads((ROOT/'docs/fingerprint_cross_source_20261024/sequence_curation_v2.json').read_text());rows=[x for x in data['rows'] if x['status']=='sequence_supported_pending_primary_header']
        inherited=[x for x in rows if any(x[r]['annotation_chain']!=x[r+'_chain'] for r in ['alpha','beta'])]
        self.assertTrue(inherited)
        for x in inherited:
            for r in ['alpha','beta']:
                if x[r]['annotation_chain']!=x[r+'_chain']:self.assertIn('same RCSB polymer entity',x[r]['annotation_inheritance'])
if __name__=='__main__':unittest.main()
