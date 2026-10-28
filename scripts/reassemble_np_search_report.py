"""Recover the exact complete broad-search report from the public gzip payload."""
from pathlib import Path
import argparse,gzip,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    base=ROOT/'docs/fingerprint_np_boundaries_20261029/broad_search';m=json.loads((base/'report_packaging.json').read_text());packed=(base/m['gzip_file']).read_bytes();assert hashlib.sha256(packed).hexdigest()==m['gzip_sha256'];raw=gzip.decompress(packed);assert len(raw)==m['raw_bytes'] and hashlib.sha256(raw).hexdigest()==m['raw_sha256'];a.out.write_bytes(raw);print(m['raw_sha256'])
if __name__=='__main__':main()
