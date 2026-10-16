"""Frozen data release contracts. Never infer task semantics from filenames."""
from __future__ import annotations
import hashlib,json
from pathlib import Path

SCHEMA_VERSION=1
STRUCTURAL_TASK='native_vs_synthetic_decoy'

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def complex_id(name):
    n=name.split('<-')[0].split('-disp')[0]
    for pre in ('af3-','tcr-','inst-','nativetcr-'):
        if n.startswith(pre):n=n[len(pre):]
    return n.split('_')[0].lower()


def load_release(path, *, data_path=None, task=STRUCTURAL_TASK, require_training=True):
    if path is None:raise ValueError('A verified --release manifest is required; legacy caches are refused')
    path=Path(path).resolve();m=json.loads(path.read_text())
    if m.get('schema_version')!=SCHEMA_VERSION or m.get('status')!='validated':raise ValueError('release is not validated')
    if m.get('task')!=task or (require_training and not m.get('training_allowed',False)):raise ValueError('release task not eligible for this training endpoint')
    if not m.get('code_files') or not m.get('files'):raise ValueError('release provenance missing')
    code=Path(__file__).resolve().parents[3]
    for name,digest in m['code_files'].items():
        if not (code/name).resolve().is_relative_to(code) or sha256(code/name)!=digest:raise ValueError('builder/consumer code changed since data acceptance: '+name)
    for name,record in m['files'].items():
        p=(path.parent/name).resolve()
        if not p.is_relative_to(path.parent) or sha256(p)!=record['sha256']:raise ValueError('release artifact hash mismatch: '+name)
    for name,digest in m.get('source_structures',{}).items():
        if sha256(name)!=digest:raise ValueError('source structure changed: '+name)
    if 'input_sources.json' in m['files']:
        for name,record in json.loads((path.parent/'input_sources.json').read_text()).items():
            if sha256(name)!=record['sha256']:raise ValueError('source evidence changed: '+name)
    actual=(path.parent/m['dataset']).resolve()
    if not actual.is_relative_to(path.parent) or m['dataset'] not in m['files']:raise ValueError('dataset is outside approved artifacts')
    if 'splits.json' in m['files'] and json.loads((path.parent/'splits.json').read_text())!=m.get('splits'):raise ValueError('manifest/frozen split mismatch')
    if data_path is not None and Path(data_path).resolve()!=actual:raise ValueError('data argument is not the approved dataset')
    return m,actual


def frozen_indices(manifest, data):
    splits=manifest['splits'];seen=set()
    for fold in ('train','validation','calibration','test'):
        idx=splits[fold]
        if not idx or len(set(idx))!=len(idx) or seen&set(idx):raise ValueError('invalid/non-disjoint frozen folds')
        if set(data['labels'][i] for i in idx)!={0,1}:raise ValueError('each fold needs both structural classes')
        if any(data['records'][i]['fold']!=fold for i in idx):raise ValueError('fold/record mismatch')
        seen.update(idx)
    if seen!=set(range(len(data['graphs']))):raise ValueError('frozen split does not cover dataset')
    return splits
