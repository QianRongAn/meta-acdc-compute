import os,subprocess,json
from pathlib import Path
root=Path('/home/administrator/docs/np105-raw-20261029'); selected=root/'pilot32_selected'
r=json.loads((selected/'selection_receipt.json').read_text());env=dict(os.environ,NP105_SPADES_ROOT=str(root/'spades-tools'),LD_LIBRARY_PATH=str(root/'spades-tools/usr/lib/x86_64-linux-gnu'))
receipts=[]
for bc in r['top_three_by_preassembly_read_count']:
 out=root/('pilot32_cell_'+bc);log=root/(out.name+'.console.log')
 cmd=['/home/administrator/docs/meta-acdc/.venv/bin/python','scripts/run_np105_spades.py','--rna','--only-assembler','-t','2','-m','4','-1',str(selected/(bc+'_R1.fastq')),'-2',str(selected/(bc+'_R2.fastq')),'-o',str(out)]
 with log.open('xb') as f:p=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
 receipts.append(dict(barcode=bc,preassembly_pairs=r['read_pairs_by_barcode'][bc],command=cmd,returncode=p.returncode,log=str(log),output=str(out)))
 print(bc,p.returncode,flush=True)
Path('docs/fingerprint_np105_raw_20261029/cell_assembly_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
