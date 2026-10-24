"""Reassemble the original public supplement from byte-identical archive segments."""
from pathlib import Path
import argparse,json,hashlib
ROOT=Path(__file__).resolve().parents[1]
def reassemble(out):
 base=ROOT/'docs/fingerprint_cross_source_extension_20261024/sources/attachment_probe';m=json.loads((base/'multipart_manifest.json').read_text());pieces=[]
 for row in m['parts']:
  raw=(base/row['file']).read_bytes();assert len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256'];pieces.append(raw)
 raw=b''.join(pieces);assert len(raw)==m['original_bytes'] and hashlib.sha256(raw).hexdigest()==m['original_sha256']
 with Path(out).open('xb') as stream:stream.write(raw)
 return m['original_sha256']
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();print(reassemble(a.out))
