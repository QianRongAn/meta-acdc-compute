"""Fetch public ENA FASTQs sequentially into a new task directory; verify complete MD5."""
from pathlib import Path
import argparse,hashlib,json,time,requests
FILES=[(1,3122244336,'16c9ef1a2f171c247eef4b8254d6d6eb'),(2,3154927322,'c11dd383daf5c93f49d6b44bff671429')]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    state=dict(run='SRR34543498',finished=False,files=[]);receipt=a.out/'download_receipt.json'
    for mate,expected_size,expected_md5 in FILES:
        url=f'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR345/098/SRR34543498/SRR34543498_{mate}.fastq.gz';path=a.out/f'SRR34543498_{mate}.fastq.gz';row=dict(mate=mate,url=url,expected_bytes=expected_size,expected_md5=expected_md5,bytes_received=0,complete=False);state['files'].append(row);receipt.write_text(json.dumps(state,indent=2)+'\n');started=time.monotonic();md5=hashlib.md5();sha=hashlib.sha256()
        try:
            with requests.get(url,timeout=(30,60),stream=True) as response:
                row['http_status']=response.status_code;row['content_length']=response.headers.get('Content-Length');response.raise_for_status()
                if response.status_code!=200:raise ValueError('expected complete HTTP response')
                with path.open('xb') as handle:
                    for chunk in response.iter_content(1024*1024):
                        if not chunk:continue
                        if row['bytes_received']+len(chunk)>expected_size:raise ValueError('response larger than ENA declaration')
                        handle.write(chunk);md5.update(chunk);sha.update(chunk);row['bytes_received']+=len(chunk)
                        if row['bytes_received']%(32*1024*1024)==0:
                            row['seconds']=round(time.monotonic()-started,2);receipt.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(row),flush=True)
            row.update(md5=md5.hexdigest(),sha256=sha.hexdigest(),seconds=round(time.monotonic()-started,2));assert row['bytes_received']==expected_size and row['md5']==expected_md5;row['complete']=True
        except Exception as e:
            row.update(error=type(e).__name__+': '+str(e),seconds=round(time.monotonic()-started,2));receipt.write_text(json.dumps(state,indent=2)+'\n');raise
        receipt.write_text(json.dumps(state,indent=2)+'\n')
    state['finished']=True;receipt.write_text(json.dumps(state,indent=2)+'\n');print('COMPLETE_CHECKSUMS_VERIFIED',flush=True)
if __name__=='__main__':main()
