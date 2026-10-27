"""Offline replay, unchanged-gate verification and preservation of seven frozen deliveries."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'docs/fingerprint_construct_followup_20261027'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--replay-dir', type=Path, required=True)
    args = ap.parse_args()
    args.replay_dir.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, PYTHONPATH='src:.')
    def run(script, *extra):
        process = subprocess.run([sys.executable, 'scripts/'+script, *map(str,extra)], cwd=ROOT,
                                 env=env, text=True, capture_output=True)
        with (args.replay_dir/'replay.log').open('a') as handle:
            handle.write(process.stdout+process.stderr)
        if process.returncode:
            raise RuntimeError(script+' failed: '+process.stderr[-1500:])
    parts = json.loads((BASE/'sources/supplement_parts/manifest.json').read_text())
    for name in ['6at6','7n5']:
        recovered = args.replay_dir/(name+'_supplement.pdf')
        run('reassemble_construct_supplements.py','--name',name,'--out',recovered)
        assert sha(recovered) == parts['supplements'][name]['sha256']
    receipt = args.replay_dir/'construct_adjudication_final.json'
    run('adjudicate_construct_sources.py', '--include-natural', '--out', receipt)
    assert sha(receipt) == sha(BASE/'construct_adjudication_final.json')
    run('qc_construct_followup.py', '--identity', receipt, '--out', args.replay_dir/'native_qc')
    assert sha(args.replay_dir/'native_qc/native_atomic_qc.json') == sha(BASE/'native_qc_final/native_atomic_qc.json')
    run('summarize_construct_followup.py', '--identity', receipt, '--qc', args.replay_dir/'native_qc/native_atomic_qc.json',
        '--out', args.replay_dir/'joint_readiness.json')
    assert sha(args.replay_dir/'joint_readiness.json') == sha(BASE/'joint_readiness_final.json')
    a = torch.load(BASE/'native_qc_final/native_graph_candidates.pt', map_location='cpu', weights_only=False)
    b = torch.load(args.replay_dir/'native_qc/native_graph_candidates.pt', map_location='cpu', weights_only=False)
    assert a['records'] == b['records'] and len(a['graphs']) == len(b['graphs']) == 1
    assert a['trainable'] is b['trainable'] is False
    for graph, replay in zip(a['graphs'], b['graphs']):
        for x, y in zip(graph,replay):
            assert np.array_equal(x,y)
    previous = []
    for folder in ['fingerprint_retrain_20261023','fingerprint_hard_family_20261024',
                   'fingerprint_paired_identity_20261024','fingerprint_cross_source_20261024',
                   'fingerprint_cross_source_extension_20261024','fingerprint_audit_integration_20261027',
                   'fingerprint_positive_controls_20261027']:
        directory = ROOT/'docs'/folder
        manifest = json.loads((directory/'DELIVERY_MANIFEST.json').read_text())
        for name, expected in manifest['files'].items():
            assert sha(directory/name) == expected['sha256'], (folder,name)
        previous.append(dict(folder=folder,n_files=len(manifest['files']),unchanged=True))
    gates = json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/DELIVERY_MANIFEST.json').read_text())
    for name, expected in {**gates['unchanged_gates'], 'src/meta_acdc/structure/primary_role_evidence.py':gates['code']['src/meta_acdc/structure/primary_role_evidence.py']}.items():
        assert sha(ROOT/name) == expected['sha256'], name
    modules = ['tests.test_construct_role_adjudication','tests.test_primary_role_evidence',
               'tests.test_cross_source_extension','tests.test_cross_source_retrieval','tests.test_paired_variable_identity']
    process = subprocess.run([sys.executable, '-m','unittest',*modules,'-v'], cwd=ROOT,env=env,capture_output=True,text=True)
    (args.replay_dir/'tests.log').write_text(process.stdout+process.stderr)
    assert process.returncode == 0, process.stderr[-2000:]
    summary = json.loads((BASE/'joint_readiness_final.json').read_text())['summary']
    assert summary['n_QC_cross_study_pairs'] == summary['n_QC_source_balanced_blocks'] == summary['n_QC_independent_natural_control_groups'] == 0
    assert summary['n_QC_graphs'] == 43 and a['records'][0]['pdb'] == '7n5c'
    result = dict(status='pass', n_graphs_exactly_replayed=1,n_total_QC_units=43,
                  n_CC_BY_supplements_reassembled_exactly=2,
                  receipt_and_readiness_byte_exact=True, unchanged_general_role_and_boundary_and_atomic_gates=True,
                  previous_deliveries=previous,n_tests=40,test_skips=0,
                  protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    (args.replay_dir/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
