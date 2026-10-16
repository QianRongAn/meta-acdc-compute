"""Read-only pretraining checks. A pass is necessary, not a clinical validation.

Use the existing project Python environment:
  python scripts/pretrain_audit.py --project-root /path/to/meta-acdc --report audit.json
Never repairs or overwrites input data. Exit 2 means training must not proceed.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def audit(root: Path) -> dict:
    import numpy as np
    import torch
    sys.path[:0] = [str(root / 'src'), str(root / 'src/meta_acdc/data')]
    from meta_acdc.structure.graph import parse_pdb, classify_chains
    from meta_acdc.structure.family_map import receptor_sequences
    from meta_acdc.models.train_domain_adapt import pdb_group
    import ingest_iedb
    checks = []
    def check(name, passed, **detail):
        checks.append(dict(name=name, passed=bool(passed), **detail))
    sha = subprocess.run(['git','-C',str(root),'rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip()
    processed = root / 'data/processed'
    vdj = processed / 'vdjdb.clean.tsv'
    if vdj.exists():
        counts = collections.Counter()
        with vdj.open() as f:
            reader = csv.DictReader(f, delimiter='\t')
            fields = reader.fieldnames or []
            for row in reader:
                counts[row.get('label','')] += 1
        check('vdjdb_not_confidence_derived_labels', 'confidence' in fields and sum(counts.values()) > 0 and not any(v for k,v in counts.items() if k), label_counts=dict(counts), confidence_column='confidence' in fields)
    else:
        check('vdjdb_not_confidence_derived_labels', False, reason='missing input')
    raw = root / 'data/raw/iedb/mhc_ligand_full.csv'
    if raw.exists():
        with raw.open() as f:
            header = f.readline()
        delimiter = '\t' if header.count('\t') > header.count(',') else ','
        with raw.open() as f:
            reader = csv.reader(f, delimiter=delimiter)
            rows = [next(reader) for _ in range(4)]
        with tempfile.TemporaryDirectory(prefix='meta-iedb-audit-') as td:
            source, dest = Path(td)/'excerpt.csv', Path(td)/'normalized.tsv'
            with source.open('w', newline='') as f:
                csv.writer(f, delimiter=delimiter).writerows(rows)
            try:
                ingest_iedb.ingest(source, dest)
                check('iedb_real_format_readable', True, delimiter=repr(delimiter))
            except Exception as exc:
                check('iedb_real_format_readable', False, delimiter=repr(delimiter), error=f'{type(exc).__name__}: {exc}')
    else:
        check('iedb_real_format_readable', False, reason='missing input')
    path = processed / 'egnn_dataset.pt'
    if path.exists():
        data = torch.load(path, map_location='cpu', weights_only=False)
        check('graph_source_provenance_present', 'provenance' in data, keys=list(data))
        check('graph_arrays_aligned', len(data['graphs']) > 0 and len(data['graphs']) == len(data['labels']) == len(data['pdbs']))
        signatures = collections.defaultdict(list)
        no_peptide, invalid = [], []
        for name, label, graph in zip(data['pdbs'], data['labels'], data['graphs']):
            nf, co, ei, ef = (np.asarray(a) for a in graph)
            if not all(np.isfinite(a).all() for a in (nf,co,ei,ef)):
                invalid.append(name)
            if nf.ndim != 2 or nf.shape[1] <= 24 or not np.any(nf[:,24] == 1):
                no_peptide.append(name)
            digest = hashlib.sha256()
            for a in (nf,co,ei,ef):
                digest.update(str(a.shape).encode())
                digest.update(a.tobytes())
            signatures[digest.hexdigest()].append((name,int(label)))
        conflicts = [v for v in signatures.values() if len({x[1] for x in v}) > 1]
        check('graph_values_finite', not invalid, examples=invalid[:20])
        check('every_graph_contains_peptide_nodes', not no_peptide, n=len(no_peptide), examples=no_peptide[:20])
        check('no_identical_graph_with_opposite_labels', not conflicts, groups=len(conflicts), examples=conflicts[:20])
    else:
        check('graph_input_exists', False)
    mapping = processed / 'pdb_family_map.tsv'
    if mapping.exists():
        mismatches, missing = [], []
        for row in csv.DictReader(mapping.open(), delimiter='\t'):
            source = root/'data/raw/structures'/f"{row['pdb']}.pdb"
            if not source.exists():
                missing.append(row['pdb'])
                continue
            residues = parse_pdb(source)
            sequences = receptor_sequences(residues, classify_chains(residues))
            key = hashlib.sha1('|'.join(sequences).encode()).hexdigest()[:12]
            if row['receptor_key'] != key:
                mismatches.append(row['pdb'])
        check('family_map_matches_current_source', not mismatches and not missing, mismatched=len(mismatches), missing=missing, examples=mismatches[:20])
    else:
        check('family_map_matches_current_source', False, reason='missing map')
    groups = [pdb_group(x) for x in ['5brz','5brz_0','af3-5brz_EVDPIGHLY']]
    check('domain_adaptation_complex_identity_grouping', len(set(groups)) == 1, groups=groups)
    clinical = processed/'clinical_index_cases.tsv'
    if clinical.exists():
        wrong = [x for x in csv.DictReader(clinical.open(),delimiter='\t') if 'A118T' in x.get('tcr','') and x.get('pdb')]
        check('clinical_A118T_not_linked_to_5BRZ', not wrong, examples=wrong)
    else:
        check('clinical_A118T_not_linked_to_5BRZ', False, reason='missing cases')
    check('frozen_release_manifest_present', (processed/'release_manifest.json').exists(), reason='Manifest contents and source/split hashes still need independent validation; file existence alone does not certify a release.')
    passed = all(c['passed'] for c in checks)
    return dict(source_commit=sha, project_root=str(root), machine_checks_passed=passed, training_allowed=False, decision='NO_GO' if not passed else 'REQUIRES_INDEPENDENT_ACCEPTANCE', checks=checks, note='Read-only necessary checks; no training authorization. Full structure identity, QC coverage, split isolation, task-label semantics and any release manifest need independent acceptance.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    p.add_argument("--release",type=Path)
    args = p.parse_args()
    try:
        if args.release:
            from meta_acdc.data.release import load_release
            from meta_acdc.data.validate_release import validate
            m,_=load_release(args.release,require_training=False)
            v=validate(args.release.parent)
            result=dict(machine_checks_passed=v['passed'],training_allowed=v['passed'] and m.get('training_allowed',False),binding_training_allowed=False,decision=('GO_STRUCTURAL_ONLY' if m.get('training_allowed',False) else 'WAIT_TRAINING_RELEASE_DECISION') if v['passed'] else 'NO_GO',task=m['task'],validation=v)
        else:result = audit(args.project_root.resolve())
    except Exception as exc:
        result = dict(machine_checks_passed=False,training_allowed=False,decision='NO_GO',error=f'{type(exc).__name__}: {exc}')
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},ensure_ascii=False))
    for item in result.get('checks',[]):
        print(('PASS ' if item['passed'] else 'FAIL ')+item['name'])
    return 0 if result.get('training_allowed') else 2

if __name__ == '__main__':
    raise SystemExit(main())
