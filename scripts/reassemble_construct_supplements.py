"""Recover unchanged CC BY source PDFs from transport-only segments."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def reassemble(name, out):
    base = ROOT/'docs/fingerprint_construct_followup_20261027/sources/supplement_parts'
    item = json.loads((base/'manifest.json').read_text())['supplements'][name]
    pieces = []
    for part in item['parts']:
        data = (base/part['file']).read_bytes()
        if len(data) != part['bytes'] or hashlib.sha256(data).hexdigest() != part['sha256']:
            raise ValueError('supplement segment integrity failure')
        pieces.append(data)
    data = b''.join(pieces)
    if len(data) != item['bytes'] or hashlib.sha256(data).hexdigest() != item['sha256']:
        raise ValueError('original supplement integrity failure')
    with Path(out).open('xb') as handle:
        handle.write(data)
    return item['sha256']


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', choices=['6at6','7n5'], required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(reassemble(args.name,args.out))
