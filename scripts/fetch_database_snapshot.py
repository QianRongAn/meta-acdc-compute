"""Replay frozen public source URLs with exact checksums; no login or installation.

Mutable upstream payloads may change. A mismatch stops ingestion; never replace
expected hashes to make a changed upstream dataset pass the frozen release.
"""
import argparse,concurrent.futures,datetime,hashlib,json,urllib.request
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-index',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();index=json.loads(a.source_index.read_text());jobs=[]
    for folder,entries in index['snapshots'].items():
        root=(a.out/folder).resolve()
        if not root.is_relative_to(a.out.resolve()):raise ValueError('snapshot path outside output')
        root.mkdir(parents=True,exist_ok=True)
        for entry in entries:jobs.append((root,entry))
    def one(job):
        root,entry=job;path=(root/entry['file']).resolve()
        if not path.is_relative_to(root):raise ValueError('source path outside snapshot')
        if path.exists():raw=path.read_bytes()
        else:
            with urllib.request.urlopen(urllib.request.Request(entry['url'],headers={'User-Agent':'meta-acdc-database-research/1.0'}),timeout=45) as response:raw=response.read(entry['bytes']+1)
        if len(raw)!=entry['bytes'] or hashlib.sha256(raw).hexdigest()!=entry['sha256']:raise ValueError('changed upstream source: '+entry['url'])
        if not path.exists():path.write_bytes(raw)
        return entry['file']
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for i,name in enumerate(pool.map(one,jobs),1):
            if i%20==0:print('verified',i,'/',len(jobs),flush=True)
    for folder,entries in index['snapshots'].items():(a.out/folder/'retrieval_manifest.json').write_text(json.dumps(entries,indent=2)+'\n')
    print('All frozen payloads verified; source dates in copied manifests refer to the original snapshot acquisition.')
if __name__=='__main__':main()
