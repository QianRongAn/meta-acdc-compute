"""Archive public primary-source requests with status, hashes and failures; no training."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import requests

ROOT = Path(__file__).resolve().parents[1]


def fetch(out, item):
    name, url = item
    path = out / name
    if path.exists():
        raise FileExistsError(path)
    row = dict(file=name, url=url, requested_at=datetime.now(timezone.utc).isoformat())
    try:
        response = requests.get(url, timeout=70, headers={'User-Agent': 'MetaACDC-construct-audit/1.0'})
        path.write_bytes(response.content)
        row.update(status=response.status_code, final_url=response.url,
                   content_type=response.headers.get('Content-Type'), bytes=len(response.content),
                   sha256=hashlib.sha256(response.content).hexdigest())
    except requests.RequestException as error:
        row.update(status='request_failed', error=str(error))
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--requests', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    items = json.loads(args.requests.read_text())
    if len({x[0] for x in items}) != len(items):
        raise ValueError('duplicate output names')
    for name, url in items:
        if Path(name).name != name or not url.startswith('https://'):
            raise ValueError('expected basename and HTTPS source')
    args.out.mkdir(parents=True, exist_ok=False)
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(lambda item: fetch(args.out, item), items))
    (args.out / 'manifest.json').write_text(json.dumps(rows, indent=2)+'\n')
    print(json.dumps([{k: r.get(k) for k in ['file', 'status', 'bytes', 'error']} for r in rows], indent=2))


if __name__ == '__main__':
    main()
