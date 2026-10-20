"""Verify the inspection-only repair bundle. This command never authorizes training."""
from __future__ import annotations
import argparse,json
from pathlib import Path
from meta_acdc.data.release import sha256,load_release
from meta_acdc.data.repair_release import validate_native_inventory
from meta_acdc.data.validate_release import validate


def verify(folder):
    root=Path(__file__).resolve().parents[1]
    manifest=json.loads((folder/'asset_manifest.json').read_text())
    if manifest.get('status')!='inspection_only' or any(manifest.get(k) is not False for k in
        ['training_allowed','binding_training_allowed','clinical_training_allowed']):
        raise ValueError('inspection bundle cannot become a training release')
    for name,record in manifest['files'].items():
        path=(folder/name).resolve()
        if not path.is_relative_to(folder.resolve()) or sha256(path)!=record['sha256']:
            raise ValueError('asset artifact changed: '+name)
    for name,digest in manifest['code_files'].items():
        path=(root/name).resolve()
        if not path.is_relative_to(root) or sha256(path)!=digest:raise ValueError('repair/consumer code changed: '+name)
    for name,digest in manifest['sources'].items():
        if sha256(name)!=digest:raise ValueError('raw or frozen source changed: '+name)
    native=validate_native_inventory(folder);pilot=validate(folder)
    if not native['passed'] or not pilot['atomic_engineering_passed']:
        raise ValueError('atomic reconstruction checks failed')
    try:load_release(folder/'asset_manifest.json',require_training=False)
    except ValueError:rejected=True
    else:raise ValueError('normal training loader incorrectly accepts unlabelled inspection assets')
    return dict(native_inventory=native,paired_control_pilot=pilot,
                normal_release_loader_refuses_bundle=rejected,training_allowed=False)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    if args.out.exists():raise ValueError('refusing to overwrite verification report')
    report=verify(args.bundle.resolve());args.out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='paired_control_pilot'},indent=2))


if __name__=='__main__':main()
