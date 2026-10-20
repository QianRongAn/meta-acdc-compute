"""Offline primary-sequence reconstruction with explicit clinician/paper-level curated facts.
Run with PYTHONPATH=src; input snapshots are preserved locally and source_manifest.json tracks provenance.
Never evaluates a model or changes a frozen panel in place.
"""
import argparse,hashlib,json
from pathlib import Path
from meta_acdc.data.clinical_primary import foy_receptors,e7_construct,validate_bundle
p=argparse.ArgumentParser();p.add_argument('--snapshot',required=True,type=Path);p.add_argument('--curation',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
a=p.parse_args()
if a.out.exists():raise ValueError('refusing to overwrite a release')
manifest=json.loads((a.curation/'source_manifest.json').read_text())
for x in manifest:
    if not x.get('file') or not x.get('sha256') or x.get('status')=='failed':continue
    f=(a.snapshot/x['file']).resolve()
    if not f.is_relative_to(a.snapshot.resolve()) or hashlib.sha256(f.read_bytes()).hexdigest()!=x['sha256']:raise ValueError('source changed '+x['file'])
receptors=foy_receptors(a.snapshot/'Foy_Table2.xlsx')
facts=json.loads((a.curation/'functional_facts.json').read_text())
for x in receptors:x['functional_evidence']=facts[x['receptor_name']]
receptors.append(e7_construct((a.snapshot/'E7_optimized_supplement.txt').read_text()))
bundle=dict(receptors=receptors,patients=json.loads((a.curation/'patient_facts.json').read_text()))
validate_bundle(bundle);a.out.mkdir(parents=True)
(a.out/'reconstructed_bundle.json').write_text(json.dumps(bundle,indent=2)+'\n')
print(json.dumps(dict(receptors=len(receptors),patients=len(bundle['patients']),clinical_training_allowed=False)))
