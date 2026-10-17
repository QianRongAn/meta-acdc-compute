"""Primary-paper SPR ledger. Raw database candidates never equal approved training data."""
from __future__ import annotations
import argparse,copy,json,math
from collections import Counter
from pathlib import Path
from .release import sha256
from .experimental_evidence import reserved_hits
ROOT=Path(__file__).resolve().parents[3]
from .evidence_completion import CURATION
REVIEW=CURATION/'review_ledger.json'
ADDITIONS=CURATION/'reservation_additions.json'

def validate_review(rows,parent,*,check_sources=True):
    old={x['evidence_id']:x for x in parent}
    if len(old)!=len(parent) or len(rows)!=len(parent) or {x['evidence_id'] for x in rows}!=set(old):raise ValueError('review coverage/duplicates')
    for x in rows:
        original=old[x['evidence_id']]
        if x['original_record']!=original:raise ValueError('immutable original changed')
        if (x['parent_fold'],x['parent_block'])!=(original['fold'],original['split_block']):raise ValueError('inherited split changed')
        if x['status'] not in {'primary_verified_soluble_SPR','quarantine_primary_unretrieved','quarantine_temperature_unverified','quarantine_wrong_ligand_value','separate_full_length_TCR_CD3','quarantine_secondary_measurement_origin','quarantine_SPR_Octet_method_conflict','quarantine_citation_and_method_conflict'}:raise ValueError('unknown review status')
        allowed=x['status']=='primary_verified_soluble_SPR'
        if x.get('training_allowed') is not allowed:raise ValueError('quarantine bypass')
        if x.get('clinical_label') is not None or x.get('binary_binding_label') is not None:raise ValueError('endpoint promotion')
        if allowed:
            for key in ['primary_locator','primary_source_sha256','primary_source_file','orientation','instrument','fit']:
                if not x.get(key):raise ValueError('primary evidence missing '+key)
            if x['fit'] not in ['equilibrium','kinetic'] or 'Octet' in x['instrument']:raise ValueError('not verified soluble SPR')
            if x.get('temperature_K') is None or not 273<=x['temperature_K']<=320:raise ValueError('temperature not verified')
            factor={'uM':1.,'nM':.001,'pM':.000001}.get(x.get('primary_unit'))
            if factor is None or not math.isfinite(x['primary_value']) or x['primary_value']<=0:raise ValueError('invalid primary unit/value')
            value=x['primary_value']*factor
            if not math.isclose(value,original['measurement']['kd_um'],rel_tol=1e-9) or not math.isclose(value,x['primary_kd_um'],rel_tol=1e-9):raise ValueError('mismatched primary value')
        if check_sources:
            for key in ['primary','numeric','method']:
                filename=x.get(key+'_source_file');digest=x.get(key+'_source_sha256')
                if not filename:continue
                f=(ROOT/filename).resolve()
                if not f.is_relative_to(ROOT) or sha256(f)!=digest:raise ValueError('primary source changed')
    return True

def build(snapshot,output):
    snapshot=snapshot.resolve()
    if output.exists():raise ValueError('refusing to overwrite release')
    curation=CURATION
    parent_path=ROOT/'docs/clinical_primary_20261018/affinity_dataset.json'
    parent=json.loads(parent_path.read_text());review=json.loads(REVIEW.read_text());validate_review(review,parent)
    registry_path=ROOT/'docs/clinical_primary_20261018/combined_reservation_registry.json'
    clinical_path=ROOT/'docs/clinical_primary_20261018/clinical_reservations.json'
    refs=json.loads(registry_path.read_text())['receptors']+json.loads(ADDITIONS.read_text())['receptors']
    records=[]
    for x in review:
        if not x['training_allowed']:continue
        y=copy.deepcopy(x['original_record'])
        if reserved_hits(y,refs):raise ValueError('clinical/external reserved family in primary dataset')
        y['primary_review']='primary paper value/units/pair/conditions checked; no complete expression-vector claim'
        y['primary_locator']=x['primary_locator'];y['primary_source_sha256']=x['primary_source_sha256']
        y['primary_assay']=dict(instrument=x['instrument'],orientation=x['orientation'],buffer=x['buffer'],chip=x['chip'],fit=x['fit'])
        y['temperature_K']=x['temperature_K'];y['method']=x['instrument']
        for key in ['structure_assay_mhc_match','structure_based_affinity_training_allowed']:
            if key in x:y[key]=x[key]
        y['measurement']['error_um']=None if x['primary_error'] is None else x['primary_error']*{'uM':1,'nM':.001,'pM':.000001}[x['primary_unit']]
        records.append(y)
    records.sort(key=lambda x:x['evidence_id']);splits={f:[i for i,x in enumerate(records) if x['fold']==f] for f in ['train','validation','test']}
    counts={f:len(ids) for f,ids in splits.items()};blocks={f:len({records[i]['split_block'] for i in ids}) for f,ids in splits.items()}
    # No model fitting while a retained independent split lacks any measurements.
    train_ok=all(counts.values()) and all(n>=3 for n in blocks.values())
    summary=dict(n_input=len(parent),n_primary_verified=len(records),n_quarantined=len(parent)-len(records),dispositions=dict(Counter(x['status'] for x in review)),
        fold_counts=counts,fold_blocks=blocks,training_allowed=bool(train_ok),binary_training_allowed=False,clinical_training_allowed=False,
        status='frozen_primary_verified_affinity' if train_ok else 'inspection_only_insufficient_independent_splits',
        scope='record-level source verification; original folds retained, no test replenishment by reshuffle, clinical heldout cohort never consumed',
        missing_condition_policy='unreported buffer details remain explicit; unverified assay temperature fails this release policy')
    output.mkdir(parents=True)
    def save(n,x):(output/n).write_text(json.dumps(x,indent=2)+'\n')
    save('affinity_dataset.json',records);save('splits.json',splits);save('review_ledger.json',review);save('summary.json',summary)
    sourcefiles=json.loads((curation/'source_files.json').read_text());sources=[]
    for name,info in sourcefiles.items():
        f=ROOT/name
        if sha256(f)!=info['sha256']:raise ValueError('snapshot source changed')
        sources.append(dict(file=str(f.relative_to(snapshot)),sha256=info['sha256']))
    code=['src/meta_acdc/data/spr_primary.py','src/meta_acdc/data/evidence_completion.py','src/meta_acdc/data/atlas_binding.py','src/meta_acdc/data/experimental_evidence.py','src/meta_acdc/data/build_release.py','src/meta_acdc/data/release.py','src/meta_acdc/models/affinity_baseline.py']
    manifest=dict(schema_version=2,task='experimental_SPR_pKD_regression',status=summary['status'],training_allowed=bool(train_ok),binary_training_allowed=False,clinical_training_allowed=False,
        primary_review_sha256=sha256(REVIEW),clinical_additions_sha256=sha256(ADDITIONS),clinical_registry_sha256=sha256(clinical_path),
        parent_dataset_path=str(parent_path),parent_dataset_sha256=sha256(parent_path),registry_path=str(registry_path),external_registry_sha256=sha256(registry_path),
        code_files={n:sha256(ROOT/n) for n in code},files={p.name:sha256(p) for p in output.iterdir() if p.is_file()},source_snapshot=str(snapshot.resolve()),sources=sources,summary=summary)
    save('manifest.json',manifest);return summary

def main():
    p=argparse.ArgumentParser();p.add_argument('--snapshot',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args();print(json.dumps(build(a.snapshot,a.out),indent=2))
if __name__=='__main__':main()
