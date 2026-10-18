"""Withdrawn structural/clinical artifacts cannot re-enter current reports."""
from __future__ import annotations
import json
from pathlib import Path
from meta_acdc.data.release import sha256
ROOT=Path(__file__).resolve().parents[3]
CATALOG=ROOT/'data/processed/artifact_status.json'
def catalog():return json.loads(CATALOG.read_text())
def assert_not_withdrawn(path):
    path=Path(path).resolve();m=catalog()
    forbidden={x['sha256'] for x in m['withdrawn']}|set(m['recovered_historical_hashes'].values())
    if any(path==(ROOT/x['path']).resolve() for x in m['withdrawn']) or sha256(path) in forbidden:
        raise ValueError('withdrawn historical artifact; cannot produce current scores or clinical conclusions: '+str(path))
    return path

def validate_score_receipt(path,release_hash):
    path=assert_not_withdrawn(path);receipt_path=Path(str(path)+'.provenance.json')
    if not receipt_path.is_file():raise ValueError('missing score provenance receipt')
    receipt=json.loads(receipt_path.read_text())
    if receipt.get('output_sha256')!=sha256(path):raise ValueError('score receipt hash mismatch')
    if receipt.get('release_sha256')!=release_hash:raise ValueError('score receipt release mismatch')
    if receipt.get('std_axis')!='AF3_models_of_checkpoint_means':raise ValueError('expected within-job AF3-model summary')
    if len(receipt.get('seed_ids',[]))!=1 or receipt['seed_ids'][0] is None:raise ValueError('one declared seed per input is required')
    if not receipt.get('checkpoint_sha256') or not receipt.get('qc_sha256'):raise ValueError('missing checkpoint/QC provenance')
    for n,h in receipt['checkpoint_sha256'].items():
        f=Path(n);f=f if f.is_absolute() else ROOT/f
        if not f.is_file() or sha256(f)!=h:raise ValueError('checkpoint artifact unavailable or changed')
    q=Path(receipt.get('qc_path',''))
    if not q.is_file() or sha256(q)!=receipt['qc_sha256']:raise ValueError('QC artifact unavailable or changed')
    return receipt
