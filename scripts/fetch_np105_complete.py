"""Sequential full-file ENA download via the standard library; no parallel or random writes."""
from pathlib import Path
import argparse,urllib.request,hashlib,json,time
FILES=[(1,3122244336,'16c9ef1a2f171c247eef4b8254d6d6eb'),(2,3154927322,'c11dd383daf5c93f49d6b44bff671429')]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--source',choices=['ena-fastq','ncbi-sra'],default='ena-fastq');a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 state=dict(run='SRR34543498',source=a.source,complete=False,files=[]);receipt=a.out/'download_receipt.json'
 def save():receipt.write_bytes((json.dumps(state,indent=2)+'\n').encode())
 files=[(str(m),n,h,f'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR345/098/SRR34543498/SRR34543498_{m}.fastq.gz',f'SRR34543498_{m}.fastq.gz') for m,n,h in FILES] if a.source=='ena-fastq' else [('sra',6040866339,'bb8d681b918233fada8e0dbe875e4887','https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR34543498/SRR34543498','SRR34543498.sra')]
 for mate,size,expected,url,filename in files:
  row=dict(mate=mate,url=url,expected_bytes=size,expected_md5=expected,bytes_received=0,complete=False);state['files'].append(row);save();start=time.monotonic();md5=hashlib.md5();sha=hashlib.sha256()
  try:
   req=urllib.request.Request(url,headers={'Range':f'bytes=0-{size-1}','User-Agent':'MetaACDC-public-read-audit/1.0'})
   with urllib.request.urlopen(req,timeout=45) as response:
    row.update(http_status=response.status,content_length=response.headers.get('Content-Length'),content_range=response.headers.get('Content-Range'));save()
    if response.status==206 and row['content_range']!=f'bytes 0-{size-1}/{size}':raise ValueError('unexpected full-file range')
    if response.status not in [200,206]:raise ValueError('unexpected response')
    with (a.out/filename).open('xb') as handle:
     while True:
      block=response.read(1024*1024)
      if not block:break
      if row['bytes_received']+len(block)>size:raise ValueError('response exceeds declared size')
      handle.write(block);md5.update(block);sha.update(block);row['bytes_received']+=len(block)
      if row['bytes_received']%(32*1024*1024)==0:
       row['seconds']=round(time.monotonic()-start,2);save();print(json.dumps(dict(mate=mate,bytes=row['bytes_received'],seconds=row['seconds'])),flush=True)
   row.update(md5=md5.hexdigest(),sha256=sha.hexdigest());assert row['bytes_received']==size and row['md5']==expected;row['complete']=True;save()
  except Exception as e:row.update(error=type(e).__name__+': '+str(e),seconds=round(time.monotonic()-start,2));save();raise
 state['complete']=True;save();print('WHOLE_RUN_MD5_VERIFIED',flush=True)
if __name__=='__main__':main()
