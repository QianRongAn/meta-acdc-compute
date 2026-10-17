"""Atomic reconstruction invariants, with synthetic fixtures only for mechanics."""
from __future__ import annotations
import copy,json,unittest
from dataclasses import replace
import numpy as np
from meta_acdc.structure.graph import Residue,build_interface_graph_from_residues
from meta_acdc.structure.rigid_controls import (transfer_peptide,atomic_screen,possible_hbond,
    donor_acceptor,peptide_geometry,policy,policy_sha256)
from meta_acdc.data.repair_release import without_confidence


def peptide(n=4):
    rows=[]
    for i in range(n):
        x=i*3.8
        atoms=[('N',x,0.,0.),('CA',x+1.45,.6,.2),('C',x+2.5,.2,0.),
               ('O',x+2.8,1.2,.1),('CB',x+1.45,2.1,.7)]
        rows.append(Residue('P','ALA',i+1,x,0,0,atoms=atoms,role='peptide'))
    return rows


def transform(rows,q,t):
    out=[]
    for r in rows:
        atoms=[(a[0],*map(float,np.asarray(a[1:])@q.T+t)) for a in r.atoms]
        out.append(replace(r,atoms=atoms))
    return out


def xyz(rows):return np.asarray([a[1:] for r in rows if r.role=='peptide' for a in r.atoms])


class RigidTransfer(unittest.TestCase):
    def test_all_internal_distances_and_chirality_preserved(self):
        receiver=peptide();donor=peptide()
        q=np.array([[0.,-1,0],[1,0,0],[0,0,1]])
        donor=transform(donor,q,np.array([10.,2,-5]))
        unchanged=copy.deepcopy(donor)
        control,report=transfer_peptide(receiver,donor)
        a,b=xyz(donor),xyz(control)
        np.testing.assert_allclose(np.linalg.norm(a[:,None]-a[None,:],axis=-1),np.linalg.norm(b[:,None]-b[None,:],axis=-1),atol=1e-10)
        self.assertAlmostEqual(np.linalg.det(report['rotation']),1)
        self.assertEqual(donor,unchanged)
        self.assertFalse(report['validates_nonbinding'])
    def test_receiver_environment_and_residue_keys_unchanged(self):
        receiver=peptide()+[Residue('T','ALA',30,0,0,0,atoms=[('CA',30,30,30)],role='tcr')]
        old=copy.deepcopy(receiver);control,_=transfer_peptide(receiver,peptide())
        self.assertEqual(control[-1],old[-1]);self.assertEqual(receiver,old)
        self.assertEqual([(r.chain,r.resid,r.icode) for r in control],[(r.chain,r.resid,r.icode) for r in receiver])
    def test_transfer_commutes_with_global_rotation_translation(self):
        receiver=peptide();donor=peptide()
        original,_=transfer_peptide(receiver,donor)
        q=np.array([[0.,0,1],[1,0,0],[0,1,0]]);t=np.array([30.,-8,5])
        shifted,_=transfer_peptide(transform(receiver,q,t),transform(donor,q,t))
        np.testing.assert_allclose(xyz(shifted),xyz(original)@q.T+t,atol=1e-10)
    def test_backbone_break_refused(self):
        receiver=peptide();receiver[1].atoms=[(a[0],a[1]+20,a[2],a[3]) for a in receiver[1].atoms]
        with self.assertRaisesRegex(ValueError,'noncontinuous'):transfer_peptide(receiver,peptide())
    def test_mismatched_length_and_missing_donor_atoms_refused(self):
        with self.assertRaisesRegex(ValueError,'lengths'):transfer_peptide(peptide(),peptide(3))
        donor=peptide();donor[0].atoms=[a for a in donor[0].atoms if a[0]!='CB']
        with self.assertRaisesRegex(ValueError,'incomplete'):transfer_peptide(peptide(),donor)
    def test_geometry_mismatch_is_not_hidden_by_local_stitching(self):
        donor=peptide();donor[0].atoms=[(a[0],a[1],a[2]+10,a[3]) if a[0] in {'CA','O'} else a for a in donor[0].atoms]
        with self.assertRaisesRegex(ValueError,'alignment'):transfer_peptide(peptide(),donor)


class ContactChemistry(unittest.TestCase):
    def test_hbond_eligibility_uses_atom_chemistry(self):
        n=Residue('P','ALA',1,0,0,0);o=Residue('M','ASP',1,0,0,0)
        self.assertTrue(possible_hbond((n,('N',0,0,0)),(o,('OD1',2.5,0,0))))
        self.assertFalse(possible_hbond((n,('CB',0,0,0)),(o,('OD1',2.5,0,0))))
        self.assertFalse(donor_acceptor(replace(n,resname='PRO'),'N')[0])
    def screen(self,environment_atom,distance):
        rows=peptide(1)
        # Only a single environment atom is needed to verify this contact classification.
        rows.append(Residue('M','ASP',1,0,0,0,atoms=[(environment_atom,distance,0,0)],role='mhc'))
        return atomic_screen(rows)
    def test_extreme_overlap_always_rejected_even_for_hbond(self):
        result=self.screen('OD1',.5)
        self.assertFalse(result['passed']);self.assertGreater(result['gross_screen']['cross_chain_pairs_below_cutoff'],0)
    def test_carbon_contacts_do_not_get_hbond_exemption(self):
        result=self.screen('CB',2.5)
        self.assertFalse(result['passed']);self.assertGreater(result['vdw_overlap_pairs'],0)
    def test_policy_does_not_claim_full_molprobity(self):
        rows=peptide(1)+[Residue('M','ALA',1,30,30,30,atoms=[('CB',30,30,30)],role='mhc')]
        result=atomic_screen(rows);self.assertTrue(result['passed'])
        self.assertTrue(result['heavy_atom_only']);self.assertFalse(result['is_molprobity'])
    def test_confidence_removed_from_both_classes_without_mutating_sources(self):
        rows=peptide(1);rows[0].plddt=95
        masked=without_confidence(rows)
        self.assertEqual(rows[0].plddt,95);self.assertIsNone(masked[0].plddt)
    def test_legacy_builder_cannot_silently_generate_retired_controls(self):
        from meta_acdc.data.build_release import build
        with self.assertRaisesRegex(ValueError,'retired'):build(None,None)
        from meta_acdc.models.dataset import build_dataset
        with self.assertRaisesRegex(ValueError,'retired'):build_dataset(None,None)


if __name__=='__main__':unittest.main()
