"""Resume a preserved SRA prefix with native curl; never changes biological sequences."""
from pathlib import Path
import argparse,hashlib,json,subprocess,time,urllib.request,shutil,re,os
SIZE=6040866339;MD5='bb8d681b918233fada8e0dbe875e4887';ETAG='"a02f10aa61040532a4fdb6806df6039c-721"'
URL='https://sra-pub-run-odp.s3.us-east-1.amazonaws.com/sra/SRR34543498/SRR34543498'
def fingerprint(path,limit=None):
 size=path.stat().st_size if limit is None else limit;left=size;md5=hashlib.md5();sha=hashlib.sha256()
 with path.open('rb') as stream:
  while left:
   block=stream.read(min(left,4*1024*1024))
   if not block:raise ValueError('file shorter than requested fingerprint interval')
   md5.update(block);sha.update(block);left-=len(block)
 return dict(bytes=size,md5=md5.hexdigest(),sha256=sha.hexdigest())
def assess(path,prefix_bytes,prefix_sha,expected_size=SIZE,expected_md5=MD5):
 row=dict(actual_bytes=path.stat().st_size,expected_bytes=expected_size,expected_md5=expected_md5,complete=False)
 if row['actual_bytes']<prefix_bytes:row['error']='preserved prefix truncated';return row
 row['prefix_unchanged']=fingerprint(path,prefix_bytes)['sha256']==prefix_sha
 if not row['prefix_unchanged']:row['error']='preserved prefix changed';return row
 if row['actual_bytes']!=expected_size:row['error']='whole file size not reached or exceeded';return row
 row.update(fingerprint(path));row['complete']=row['md5']==expected_md5
 if not row['complete']:row['error']='whole file MD5 mismatch'
 return row

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();out=args.out;checkpoint=json.loads((out/'prefix_checkpoint.json').read_text());path=out/'SRR34543498.sra';prefix=checkpoint['copied_prefix'];assert fingerprint(path)==prefix
 with (out/'transfer_lease.json').open('x') as stream:json.dump(dict(pid=os.getpid(),source=URL,initial_bytes=prefix['bytes']),stream)
 proxy=urllib.request.getproxies().get('https');curl=shutil.which('curl.exe');assert curl,'native Windows curl required';headers=out/'curl_headers.txt';receipt=out/'download_receipt.json';start=time.monotonic()
 command=[curl,'--fail','--location','--silent','--show-error','--continue-at','-','--connect-timeout','20','--speed-time','60','--speed-limit','1024','--header','If-Match: '+ETAG,'--dump-header',str(headers),'--output',str(path)]
 if proxy:command.extend(['--proxy',proxy])
 command.append(URL);state=dict(source=URL,expected_bytes=SIZE,expected_md5=MD5,preserved_prefix=prefix,system_proxy_used=bool(proxy),complete=False,process_running=True)
 def save():receipt.write_bytes((json.dumps(state,indent=2)+'\n').encode())
 save();last=0
 with (out/'curl.log').open('wb') as log:
  process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
  while process.poll() is None:
   state.update(bytes_received=path.stat().st_size,seconds=round(time.monotonic()-start,2),curl_pid=process.pid);save()
   if time.monotonic()-last>=60:print(json.dumps(dict(bytes=state['bytes_received'],expected_bytes=SIZE,seconds=state['seconds'])),flush=True);last=time.monotonic()
   time.sleep(2)
  state['curl_returncode']=process.returncode
 state['process_running']=False;state['assessment']=assess(path,prefix['bytes'],prefix['sha256']);text=headers.read_text(encoding='iso8859-1') if headers.exists() else '';ranges=re.findall(r'(?im)^Content-Range:\s*(.*?)\r?$',text);tags=re.findall(r'(?im)^ETag:\s*(.*?)\r?$',text);state['response_range']=ranges[-1].strip() if ranges else None;state['response_ETag']=tags[-1].strip() if tags else None
 state['response_range_valid']=state['response_range']==f"bytes {prefix['bytes']}-{SIZE-1}/{SIZE}";state['response_ETag_valid']=state['response_ETag']==ETAG
 state['bytes_received']=state['assessment']['actual_bytes']
 state['complete']=state['curl_returncode']==0 and state['assessment']['complete'] and state['response_range_valid'] and state['response_ETag_valid'];state['seconds']=round(time.monotonic()-start,2);save();print(json.dumps(state),flush=True)
 if not state['complete']:raise SystemExit(2)
if __name__=='__main__':main()
