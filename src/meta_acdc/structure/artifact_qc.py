"""Conservative construction screens; passing is NOT conformational or binding validation."""
from __future__ import annotations
import numpy as np
from scipy.spatial import cKDTree


def construction_screen(residues, severe_distance=1.2):
    peptide = [r for r in residues if r.role == 'peptide']
    environment = [a[1:] for r in residues if r.role != 'peptide' for a in r.atoms
                   if not a[0].upper().startswith('H')]
    patoms = [a[1:] for r in peptide for a in r.atoms if not a[0].upper().startswith('H')]
    if not patoms or not environment:
        raise ValueError('missing peptide or environment heavy atoms')
    # Cross-chain peptide/environment pairs cannot be ordinary intra-residue bonds.
    hits = cKDTree(np.asarray(patoms)).query_ball_tree(cKDTree(np.asarray(environment)), severe_distance)
    clashes = sum(len(x) for x in hits)
    bad_bonds = []
    for r in peptide:
        atoms = {a[0]:np.asarray(a[1:]) for a in r.atoms}
        if r.resname == 'GLY': continue
        if 'CA' not in atoms or 'CB' not in atoms:
            bad_bonds.append([r.chain,r.resid,'missing_CA_CB']); continue
        distance = float(np.linalg.norm(atoms['CA']-atoms['CB']))
        if not 1.1 <= distance <= 1.9:
            bad_bonds.append([r.chain,r.resid,distance])
    return dict(policy='construction-screen-v1', severe_distance_A=severe_distance,
                cross_chain_pairs_below_cutoff=clashes, ca_cb_outliers=bad_bonds,
                passed=clashes==0 and not bad_bonds,
                validates_rotamers=False, validates_nonbinding=False)


def require_construction_screen(residues):
    report = construction_screen(residues)
    if not report['passed']:
        raise ValueError('construction artifacts: ' + str(report))
    return report


def load_record_residues(record, *, strict_identity=True):
    """Reconstruct the five declared chains and verify their observed sequences."""
    from pathlib import Path
    from .graph import parse_pdb
    from .cif import parse_cif
    from .identity import chain_sequence
    path = Path(record['source_path'])
    residues = parse_pdb(path) if record['domain']=='crystal' else parse_cif(path, confidence_kind="af3_plddt")
    roles = dict(alpha_chain=('tcr','tra','alpha'), beta_chain=('tcr','trb','beta'),
                 peptide_chain=('peptide','peptide','peptide'),
                 mhc_chain=('mhc','mhc_i','mhc_sequence'), b2m_chain=('b2m','b2m','b2m_sequence'))
    chains = [record[k] for k in roles]
    if len(set(chains)) != 5: raise ValueError('chain roles must be distinct')
    selected = []
    for field,(role,kind,sequence) in roles.items():
        rs = [r for r in residues if r.chain==record[field]]
        if not rs or (strict_identity and chain_sequence(rs)!=record[sequence]):
            raise ValueError('source sequence / declared chain identity mismatch: '+field)
        for r in rs: r.role=role; r.chain_kind=kind
        selected.extend(rs)
    return selected
