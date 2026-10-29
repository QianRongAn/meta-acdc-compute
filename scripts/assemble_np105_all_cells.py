"""Assemble every whitelisted cell independently; no candidate/reference sequence enters assembly."""
import os,subprocess,json,time,hashlib,argparse
from pathlib import Path
from annotate_np105_raw_assembly import annotate
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'docs/fingerprint_np105_validation_20261029';raw=Path('/home/administrator/docs/np105-raw-20261029');selection=raw/'pilot32_selected'
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--selected',type=Path,default=selection);ap.add_argument('--out',type=Path,default=raw/'all92_cell_assemblies');ap.add_argument('--report',type=Path,default=B/'all92_assemblies.json');a=ap.parse_args();out=a.out;out.mkdir(exist_ok=False);selection_path=a.selected;s=json.loads((selection_path/'selection_receipt.json').read_text());env=dict(os.environ,NP105_SPADES_ROOT=str(raw/'spades-tools'),LD_LIBRARY_PATH=str(raw/'spades-tools/usr/lib/x86_64-linux-gnu'));rows=[]
 for bc in sorted(s['read_pairs_by_barcode']):
  start=time.monotonic();dest=out/bc;cmd=['/home/administrator/docs/meta-acdc/.venv/bin/python','scripts/run_np105_spades.py','--rna','--only-assembler','-t','2','-m','4','-1',str(selection_path/(bc+'_R1.fastq')),'-2',str(selection_path/(bc+'_R2.fastq')),'-o',str(dest)]
  with (out/(bc+'.log')).open('w') as log:p=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
  row=dict(barcode=bc,selected_read_pairs=s['read_pairs_by_barcode'][bc],exit_code=p.returncode,seconds=round(time.monotonic()-start,2),command=cmd)
  if p.returncode==0 and (dest/'transcripts.fasta').is_file():row['annotation']=annotate(dest/'transcripts.fasta')
  else:row['error']='assembly_failed_or_no_transcripts'
  rows.append(row);a.report.write_text(json.dumps(dict(scope='same bounded prefix; every whitelisted cell, not independent receptors',biological_reference_used_for_assembly=False,cells=rows,finished=len(rows)==92),indent=2)+'\n')
  if len(rows)%12==0:print('ASSEMBLED',len(rows),'OF',92,flush=True)
 print('ALL92_COMPLETE',flush=True)
if __name__=='__main__':main()
