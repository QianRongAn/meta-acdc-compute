"""Stream full-quality SRA biological read pairs into the unchanged barcode/Q20 selection rule."""
from pathlib import Path
import argparse,collections,hashlib,json,os,re,subprocess,sys,time,shutil
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.resume_np105_native import fingerprint,SIZE,MD5
from scripts.select_np105_raw_reads import trim_ends
ROOT=Path(__file__).resolve().parents[1];PRIOR=ROOT/'docs/fingerprint_np105_raw_20261029'

def read_fastq(stream):
 h=stream.readline()
 if not h:return None
 seq=stream.readline();plus=stream.readline();qual=stream.readline()
 if not all([seq,plus,qual]) or not qual.endswith('\n'):raise ValueError('incomplete exported FASTQ record')
 h,seq,plus,qual=[x.rstrip('\r\n') for x in [h,seq,plus,qual]]
 if not h.startswith('@') or not plus.startswith('+') or len(seq)!=len(qual):raise ValueError('malformed exported FASTQ')
 return h,seq,qual

def selected_pair(one,two,barcodes):
 h,s1,q1=one;h2,s2,q2=two
 if h!=h2:raise ValueError('SRA biological mate spot identifiers differ')
 bc=s1[:16]
 if bc not in barcodes:return None,'not_whitelisted'
 if len(s1)<86 or len(s2)<60:return None,'short_input'
 if min(ord(c)-33 for c in q1[:26])<20:return None,'barcode_or_UMI_quality_below_20'
 umi=s1[16:26]
 if set(umi)-set('ACGT'):return None,'noncanonical_UMI'
 s1,q1=trim_ends(s1[26:],q1[26:]);s2,q2=trim_ends(s2,q2)
 if min(len(s1),len(s2))<60:return None,'short_after_Q20_end_trim'
 if set(s1+s2)-set('ACGTN'):raise ValueError('invalid exported raw DNA base')
 return (bc,umi,[(s1,q1),(s2,q2)]),None

def check_destination(path):
 destination=path.parent.resolve()
 if not destination.is_relative_to(Path('/mnt/d')):raise ValueError('full-run outputs and temporary files must use the Windows D volume')
 free=shutil.disk_usage(destination).free
 if free<8*SIZE:raise ValueError('destination needs at least eight SRA file sizes of free space')
 return free

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--sra',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();destination_free=check_destination(args.out)
 source=fingerprint(args.sra);assert source['bytes']==SIZE and source['md5']==MD5,'full SRA size and MD5 must pass before selection';args.out.mkdir(exist_ok=False);temp=args.out/'sdk_temp';temp.mkdir();tool=args.sra.parent.parent/'sra-tools';env=dict(os.environ,LD_LIBRARY_PATH=str(tool/'usr/lib/x86_64-linux-gnu'),VDB_CONFIG=str(tool/'etc/ncbi'),NCBI_SETTINGS=str(args.sra.parent.parent/'sra-user-settings.mkfg'))
 p=subprocess.run([str(tool/'usr/bin/vdb-dump'),str(args.sra),'-T','SEQUENCE','-r'],env=env,capture_output=True,text=True,check=True);(args.out/'source_range.txt').write_text(p.stdout);m=re.fullmatch(r'id-range: first-row = ([0-9,]+), row-count = ([0-9,]+)\s*',p.stdout);assert m and int(m[1].replace(',',''))==1;expected_spots=int(m[2].replace(',',''));assert expected_spots==50653512
 command=[str(tool/'usr/bin/fasterq-dump'),str(args.sra),'--split-spot','--skip-technical','--stdout','--seq-defline','@SRR34543498.$si','--qual-defline','+','--threads','2','--temp',str(temp)]
 whitelist=json.loads((PRIOR/'barcode_whitelist.json').read_text());barcodes={x[:-2] for x in whitelist['barcodes']};assert len(barcodes)==92;old=json.loads((PRIOR/'selection_receipt.json').read_text());fixture=json.loads((ROOT/'docs/fingerprint_np105_recovery_20261029/sra_first3_raw.json').read_text());first3={}
 for row in fixture:
  offset=8;reads=[]
  for length in row['READ_LEN'][1:]:reads.append((row['READ'][offset:offset+length],''.join(chr(q+33) for q in row['QUALITY'][offset:offset+length])));offset+=length
  first3[row['row_id']]=reads
 handles={};digests={};nbytes=collections.Counter();counts=collections.Counter();umis=collections.defaultdict(set);excluded=collections.Counter();paired=0;selected=0;prefix_exact=False;started=time.monotonic();receipt=args.out/'selection_receipt.json';state=dict(destination_free_bytes_at_start=destination_free,source_file=str(args.sra),source_fingerprint=source,source_whole_MD5_verified=True,expected_spots=expected_spots,complete=False,command=command,biological_reference_used=False,bases_substituted=0,independent_acceptance_accessed=False)
 def write(name,data):
  if name not in handles:handles[name]=(args.out/name).open('xb');digests[name]=hashlib.sha256()
  handles[name].write(data);digests[name].update(data);nbytes[name]+=len(data)
 def save():
  state.update(n_source_pairs=paired,n_selected_pairs=selected,n_selected_barcodes=len(counts),seconds=round(time.monotonic()-started,2),prefix_all186_selected_outputs_byte_exact=prefix_exact);receipt.write_text(json.dumps(state,indent=2)+'\n')
 save()
 with (args.out/'sra_export.stderr').open('w') as log:
  process=subprocess.Popen(command,env=env,cwd=args.out,stdout=subprocess.PIPE,stderr=log,text=True,encoding='ascii',bufsize=1024*1024)
  try:
   while True:
    one=read_fastq(process.stdout)
    if one is None:break
    two=read_fastq(process.stdout)
    if two is None:raise ValueError('orphan exported biological read')
    spot=paired+1
    if one[0]!=f'@SRR34543498.{spot}' or two[0]!=one[0]:raise ValueError('source spot order or count disagrees')
    if len(one[1])!=150 or len(two[1])!=150:raise ValueError('unexpected biological read length; inspect SRA metadata')
    if spot in first3 and [(one[1],one[2]),(two[1],two[2])]!=first3[spot]:raise ValueError('export differs from independently read SRA sequence/quality fixture')
    paired+=1;result,reason=selected_pair(one,two,barcodes)
    if result is None:excluded[reason]+=1
    else:
     bc,umi,reads=result
     for mate,(seq,qual) in enumerate(reads,1):
      data=(one[0]+'.BC_'+bc+'.UMI_'+umi+'/'+str(mate)+'\n'+seq+'\n+\n'+qual+'\n').encode();write('pooled_R'+str(mate)+'.fastq',data);write(bc+'_R'+str(mate)+'.fastq',data)
     selected+=1;counts[bc]+=1;umis[bc].add(umi)
    if paired==old['n_complete_synchronized_prefix_pairs']:
     assert set(digests)==set(old['outputs'])
     for name,x in old['outputs'].items():assert digests[name].hexdigest()==x['sha256'] and nbytes[name]==x['bytes'],name
     prefix_exact=True;save();print('PREFIX_186_OUTPUT_HASHES_EXACT',flush=True)
    if paired%1000000==0:save();print(json.dumps(dict(source_pairs=paired,selected_pairs=selected,seconds=round(time.monotonic()-started,2))),flush=True)
   returncode=process.wait();assert returncode==0,'SRA exporter failed';assert paired==expected_spots and prefix_exact
  except BaseException as error:
   if process.poll() is None:process.terminate();process.wait(timeout=30)
   state['error']=type(error).__name__+': '+str(error);save();raise
  finally:
   process.stdout.close()
   for f in handles.values():f.close()
 state.update(complete=True,read_pairs_by_barcode=dict(counts),UMIs_by_barcode={bc:len(s) for bc,s in umis.items()},excluded_reasons=dict(excluded),technical_SRA_read1_excluded=True,outputs={name:dict(bytes=nbytes[name],sha256=d.hexdigest()) for name,d in digests.items()});save();print('FULL_RUN_SELECTION_COMPLETE',paired,selected,flush=True)
if __name__=='__main__':main()
