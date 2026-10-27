"""Apply unchanged complete class-I graph and atomic gates to adjudicated entries."""
from dataclasses import replace
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import torch

from meta_acdc.structure.graph import parse_pdb, build_interface_graph_from_residues
from meta_acdc.structure.identity import annotate, class_i_units
from meta_acdc.structure.receptor_identity_audit import seqres_sequences, is_subsequence
from meta_acdc.structure.rigid_controls import atomic_screen, policy_sha256

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'docs/fingerprint_construct_followup_20261027'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--identity', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    identity = json.loads(args.identity.read_text())
    rows, graphs, records = [], [], []
    for item in identity['rows']:
        path = ROOT / item['primary_header']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item['primary_header_sha256']
        residues = parse_pdb(path)
        if not residues:
            raise ValueError('source has no protein coordinates; header-only evidence is not atomic QC')
        annotate(residues, path)
        seqres = seqres_sequences(path.read_text())
        for role, kind in [('alpha', 'tra'), ('beta', 'trb')]:
            assert seqres[item[role+'_chain']] == item[role]['full_construct']
            for residue in residues:
                if residue.chain == item[role+'_chain']:
                    residue.role, residue.chain_kind = 'tcr', kind
        units, excluded = class_i_units(residues)
        rows.append(dict(pdb=item['pdb'], stage='unit_extraction', n_units=len(units), excluded=excluded))
        for unit, meta in units:
            result = dict(pdb=item['pdb'], alpha_chain=meta['alpha_chain'], beta_chain=meta['beta_chain'],
                          status='quarantined', reasons=[], source_sha256=item['primary_header_sha256'])
            try:
                for role in ['alpha', 'beta']:
                    if meta[role+'_chain'] != item[role+'_chain']:
                        raise ValueError('contact-selected unit differs from adjudicated author chains')
                    if not is_subsequence(meta[role], item[role]['full_construct']):
                        raise ValueError('observed sequence conflicts with deposited construct')
                graph = build_interface_graph_from_residues([replace(r, plddt=None) for r in unit])
                screen = atomic_screen(unit)
                result['atomic_qc'] = screen
                if not screen['passed']:
                    raise ValueError('existing fixed atomic screen rejected unit')
                result.update(status='native_graph_QC_pass', receptor_id=item['receptor_id'],
                              V_context=item['V_context'], source_study_key=item['source_study_key'], n_nodes=graph.n_nodes,
                              unit_metadata=dict(meta, pmhc_construct_key=hashlib.sha256(
                                  '|'.join(seqres[meta[k]] for k in ['mhc_chain','b2m_chain','peptide_chain']).encode()).hexdigest()))
                arrays = (graph.node_features, graph.node_coords, graph.edge_index, graph.edge_features)
                digest = hashlib.sha256()
                for array in arrays:
                    digest.update(str(np.asarray(array).shape).encode()); digest.update(np.asarray(array).tobytes())
                result['graph_sha256'] = digest.hexdigest()
                graphs.append(arrays)
                records.append(dict(result, peptide_mask=graph.peptide_mask, node_keys=graph.node_keys))
            except ValueError as error:
                result['reasons'].append(str(error))
            rows.append(result)
    torch.save(dict(task='sequence-adjudicated native candidates; no training release', trainable=False,
                    binding_labels_available=False, graphs=graphs, records=records,
                    atomic_policy_sha256=policy_sha256()), args.out/'native_graph_candidates.pt')
    report = dict(rows=rows, n_graphs=len(graphs), n_pdbs_with_pass=len({r['pdb'] for r in records}),
                  atomic_policy_sha256=policy_sha256(), model_training_started=False)
    (args.out/'native_atomic_qc.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(n_graphs=len(graphs), decisions=[{k: r.get(k) for k in ['pdb','stage','status','reasons','n_units']} for r in rows]), indent=2))


if __name__ == '__main__':
    main()
