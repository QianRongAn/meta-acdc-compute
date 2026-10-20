"""Draw a validation-only figure from a hash-bound machine-readable run report."""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from meta_acdc.data.release import sha256
ROOT=Path(__file__).resolve().parents[3]

def load_validation_comparison(report,manifest):
    report=Path(report).resolve();manifest=Path(manifest).resolve()
    if not report.is_relative_to(ROOT):raise ValueError('report must be a repository artifact')
    seal=json.loads(manifest.read_text());key=str(report.relative_to(ROOT))
    if seal.get('artifacts',{}).get(key)!=sha256(report):raise ValueError('report is not bound by supplied run manifest')
    data=json.loads(report.read_text())
    if data.get('test_used_for_selection') is not False or data.get('consumed_test_evaluated') is not False or data.get('new_test_evaluated') is not False:raise ValueError('expected registered development-only comparison')
    values=[(x['name'],float(x['validation']['rmse_pkd'])) for x in data['candidates']]
    values.append(('training_mean',float(data['train_mean_validation']['rmse_pkd'])))
    if len({x[0] for x in values})!=len(values) or any(not math.isfinite(v) or v<0 for _,v in values):raise ValueError('invalid validation comparison')
    return values

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report',type=Path,required=True)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--out-dir',type=Path,required=True)
    a=p.parse_args(argv);values=load_validation_comparison(a.report,a.manifest)
    if a.out_dir.exists():raise ValueError('refusing to overwrite figures')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    a.out_dir.mkdir(parents=True,exist_ok=False)
    fig,ax=plt.subplots(figsize=(8,3.8));names,ys=zip(*values)
    bars=ax.barh(names,ys,color='#397ab8');ax.invert_yaxis();ax.set_xlim(0,max(ys)*1.20)
    ax.set_xlabel('Block-weighted validation RMSE (pKD; lower is better)')
    ax.set_title('Registered development comparison — no independent acceptance',fontsize=11)
    for bar,v in zip(bars,ys):ax.text(v+.01,bar.get_y()+bar.get_height()/2,f'{v:.4f}',va='center')
    fig.tight_layout();fig.savefig(a.out_dir/'validation_comparison.png',dpi=160);fig.savefig(a.out_dir/'validation_comparison.pdf');plt.close(fig)
    receipt=dict(scope='development validation only; no clinical/independent claim',input_report_sha256=sha256(a.report),input_manifest_sha256=sha256(a.manifest),values=[dict(candidate=k,validation_rmse_pkd=v) for k,v in values],outputs={f.name:sha256(f) for f in a.out_dir.iterdir() if f.is_file()})
    (a.out_dir/'provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return 0
if __name__=='__main__':raise SystemExit(main())
