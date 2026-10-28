"""Bind current raw coordinates to unchanged previously admitted receptor identities."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_adjudication_20261029'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    old=json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/source_balanced_design.json').read_text());rows=[]
    for r in old['rows']:
        if r['pdb'] not in ['7n1f','7rtr','4jfe','4jfh','6eqb'] or r['status']!='identity_admitted_structure_QC_pending':continue
        item=json.loads(json.dumps(r));path=BASE/'sources/followup'/(r['pdb']+'.pdb')
        import gzip
        raw=path.read_bytes() if path.exists() else gzip.decompress(path.with_suffix('.pdb.gz').read_bytes())
        item['primary_header']=path.relative_to(ROOT).as_posix();item['primary_header_sha256']=hashlib.sha256(raw).hexdigest();rows.append(item)
    assert len(rows)==5
    a.out.write_text(json.dumps(dict(rows=rows,identity_basis='unchanged previously admitted complete identities; current coordinate bytes only',independent_acceptance_accessed=False),indent=2)+'\n')
if __name__=='__main__':main()
