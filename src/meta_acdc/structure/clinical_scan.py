"""Historical clinical score interpretation is withdrawn; no clinical model exists.

The optional audit output contains status/provenance, never replacement scores.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from meta_acdc.structure.artifact_status import catalog

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores',type=Path,help='legacy argument; clinical scoring remains blocked')
    p.add_argument('--clinical',type=Path,help='legacy argument; clinical scoring remains blocked')
    p.add_argument('--out',type=Path)
    p.add_argument('--audit-withdrawal',action='store_true',help='write a fresh withdrawal receipt, without scores')
    a=p.parse_args(argv)
    if not a.audit_withdrawal:raise ValueError('Clinical scoring blocked: MAG-IC3 is not clinical a3a; historical scores/QC do not validate clinical risk. See docs/delivery_audit_20261018/README.md')
    if a.out is None:raise ValueError('audit requires a new --out path')
    if a.out.exists():raise ValueError('refusing to overwrite clinical artifact')
    m=catalog()
    result=dict(status='WITHDRAWN',clinical_model_released=False,scores_emitted=0,hazard_flags_emitted=0,reference_facts=m['active_references'],historical_artifacts=m['withdrawn'])
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('withdrawal receipt only; zero clinical scores or hazard judgments')
    return 0
if __name__=='__main__':raise SystemExit(main())
