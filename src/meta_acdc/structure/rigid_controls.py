"""Whole-observed-peptide controls; never physical nonbinding ground truth.

One proper rigid transformation preserves all donor peptide internal distances,
bond angles, torsions and chirality; this avoids per-residue backbone stitching.
Source conformations can still have errors. The conservative heavy-atom screen
is not MolProbity and does not certify thermodynamic or biological validity.
"""
from __future__ import annotations
from dataclasses import replace
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .identity import missing_heavy_atoms

POLICY_PATH=Path(__file__).resolve().parents[3]/'docs/structure_repair_20261017/policy_v2.json'


def policy():
    return json.loads(POLICY_PATH.read_text())


def policy_sha256():
    return hashlib.sha256(POLICY_PATH.read_bytes()).hexdigest()


def peptide_rows(residues):
    rows=sorted([r for r in residues if r.role=='peptide'],key=lambda r:(r.resid,r.icode))
    if not rows or len({r.chain for r in rows})!=1:
        raise ValueError('exactly one explicit peptide required')
    if any(missing_heavy_atoms(r) for r in rows):raise ValueError('incomplete donor/receiver peptide')
    return rows


def peptide_geometry(residues):
    rows=peptide_rows(residues);bad=[]
    lo,hi=policy()['peptide_CN_bond_range_A']
    for a,b in zip(rows,rows[1:]):
        aa={x[0]:np.asarray(x[1:]) for x in a.atoms};bb={x[0]:np.asarray(x[1:]) for x in b.atoms}
        distance=float(np.linalg.norm(aa['C']-bb['N']))
        if not lo<=distance<=hi:bad.append(dict(left=a.resid,right=b.resid,distance_A=distance))
    return bad


def transfer_peptide(receiver,donor):
    native=peptide_rows(receiver);foreign=peptide_rows(donor)
    if len(native)!=len(foreign):raise ValueError('peptide lengths differ')
    if peptide_geometry(receiver) or peptide_geometry(donor):raise ValueError('noncontinuous peptide backbone')
    from ..models.dataset import _kabsch
    def backbone(rows):
        return np.asarray([dict((a[0],a[1:]) for a in r.atoms)[n]
                           for r in rows for n in ['N','CA','C','O']])
    source=backbone(foreign);target=backbone(native)
    rotation,translation=_kabsch(source,target)
    rmsd=float(np.sqrt(np.mean(np.sum((source@rotation.T+translation-target)**2,axis=1))))
    if rmsd>policy()['max_backbone_fit_rmsd_A']:raise ValueError('backbone alignment outside preset tolerance')
    if not np.isclose(np.linalg.det(rotation),1,atol=1e-10):raise ValueError('improper rotation')
    transplanted=[]
    for a,b in zip(native,foreign):
        xyz=np.asarray([x[1:] for x in b.atoms])@rotation.T+translation
        atoms=[(old[0],*map(float,coord)) for old,coord in zip(b.atoms,xyz)]
        center=xyz.mean(0)
        transplanted.append(replace(a,resname=b.resname,atoms=atoms,x=float(center[0]),y=float(center[1]),z=float(center[2]),plddt=None))
    # All pair distances (including between peptide residues) must be preserved.
    before=np.asarray([a[1:] for r in foreign for a in r.atoms]);after=np.asarray([a[1:] for r in transplanted for a in r.atoms])
    d0=np.linalg.norm(before[:,None]-before[None,:],axis=-1)
    d1=np.linalg.norm(after[:,None]-after[None,:],axis=-1)
    error=float(np.max(np.abs(d0-d1)))
    if error>policy()['internal_distance_preservation_max_error_A']:raise ValueError('donor geometry not preserved')
    swaps={id(a):b for a,b in zip(native,transplanted)}
    result=[swaps.get(id(r),r) for r in receiver]
    return result,dict(policy_id=policy()['policy_id'],rotation=rotation.tolist(),translation=translation.tolist(),
                       backbone_fit_rmsd_A=rmsd,internal_distance_max_error_A=error,
                       donor_geometry_preserved=True,relaxed=False,validates_nonbinding=False)


# Potential donors/acceptors only; tautomer/protonation and H angles remain unknown.
_SIDE_DONORS={'ARG':{'NE','NH1','NH2'},'LYS':{'NZ'},'ASN':{'ND2'},'GLN':{'NE2'},
              'HIS':{'ND1','NE2'},'SER':{'OG'},'THR':{'OG1'},'TYR':{'OH'},
              'TRP':{'NE1'},'CYS':{'SG'}}
_SIDE_ACCEPTORS={'ASP':{'OD1','OD2'},'GLU':{'OE1','OE2'},'ASN':{'OD1'},'GLN':{'OE1'},
                 'HIS':{'ND1','NE2'},'SER':{'OG'},'THR':{'OG1'},'TYR':{'OH'},
                 'CYS':{'SG'},'MET':{'SD'}}


def donor_acceptor(residue,atom_name):
    name=residue.resname.upper();atom_name=atom_name.upper()
    donor=(atom_name=='N' and name!='PRO') or atom_name in _SIDE_DONORS.get(name,set())
    acceptor=atom_name in {'O','OXT'} or atom_name in _SIDE_ACCEPTORS.get(name,set())
    return donor,acceptor


def possible_hbond(a,b):
    da,aa=donor_acceptor(a[0],a[1][0]);db,ab=donor_acceptor(b[0],b[1][0])
    return (da and ab) or (db and aa)


def atomic_screen(residues):
    from .artifact_qc import construction_screen
    initial=construction_screen(residues);settings=policy()
    pep=[(r,a) for r in residues if r.role=='peptide' for a in r.atoms if not a[0].startswith('H')]
    env=[(r,a) for r in residues if r.role!='peptide' for a in r.atoms if not a[0].startswith('H')]
    radii=settings['bondi_vdw_radii_A']
    def element(atom):
        e=atom[0][0].upper()
        if e not in radii:raise ValueError('unsupported atom element: '+atom[0])
        return e
    pr=np.asarray([radii[element(a)] for r,a in pep]);er=np.asarray([radii[element(a)] for r,a in env])
    px=np.asarray([a[1:] for r,a in pep]);ex=np.asarray([a[1:] for r,a in env])
    if not np.isfinite(px).all() or not np.isfinite(ex).all():raise ValueError('nonfinite atom coordinates')
    hits=cKDTree(px).query_ball_tree(cKDTree(ex),2*max(radii.values()))
    bad=[];max_overlap=0.;potential_hbonds=0;short_hbonds=0
    for i,js in enumerate(hits):
        for j in js:
            distance=float(np.linalg.norm(px[i]-ex[j]));overlap=float(pr[i]+er[j]-distance)
            max_overlap=max(max_overlap,overlap)
            candidate_hbond=possible_hbond(pep[i],env[j])
            if candidate_hbond and overlap>settings['max_cross_chain_vdw_overlap_A']+1e-10:
                if distance>=settings['potential_hbond_min_distance_A']:
                    potential_hbonds+=1;continue
                short_hbonds+=1
            if overlap>settings['max_cross_chain_vdw_overlap_A']+1e-10:
                a,b=pep[i],env[j]
                bad.append(dict(peptide=[a[0].chain,a[0].resid,a[0].icode,a[1][0]],
                                environment=[b[0].chain,b[0].resid,b[0].icode,b[1][0]],
                                distance_A=distance,overlap_A=overlap))
    bad_bonds=peptide_geometry(residues)
    return dict(policy_id=settings['policy_id'],passed=initial['passed'] and not bad and not bad_bonds,
                gross_screen=initial,vdw_overlap_pairs=len(bad),max_vdw_overlap_A=max_overlap,
                worst_pairs=sorted(bad,key=lambda x:-x['overlap_A'])[:8],peptide_CN_outliers=bad_bonds,
                potential_hbond_contacts=potential_hbonds,too_short_donor_acceptor_contacts=short_hbonds,
                heavy_atom_only=True,is_molprobity=False,validates_nonbinding=False)
