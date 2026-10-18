"""Export corrected clinical/structural reference facts, never prediction labels.

MAG-IC3, clinical a3a and the HLA-A02 case remain separate constructs/cases.
Historical CASES with inferred clinical labels are retired in the audit archive.
"""
from __future__ import annotations
import argparse,csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
REFERENCE=ROOT/'docs/delivery_audit_20261018/corrected_clinical_references.tsv'

def load_references():
    with REFERENCE.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    if not rows or any(x['record_status']!='REFERENCE_ONLY_NOT_GOLD_STANDARD' or x['clinical_prediction_label'] for x in rows):raise ValueError('clinical references cannot supply prediction labels')
    return rows

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args(argv)
    if a.out.exists():raise ValueError('refusing to overwrite reference facts')
    rows=load_references();a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    return 0
if __name__=='__main__':raise SystemExit(main())
