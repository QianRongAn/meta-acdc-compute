from pathlib import Path
import sys,json,importlib.metadata,hashlib
sys.path.insert(0,'/tmp/np-classic-deps-20261029')
import anarci
b=Path('docs/fingerprint_np105_raw_20261029');inputs=[]
for f in sorted(b.glob('pilot32_cell_*_annotation.json')):
 bc=f.name.split('_')[2]
 r=json.loads(f.read_text())
 for i,x in enumerate(r['matches']):inputs.append(dict(name=bc+'_'+x['role']+'_'+str(i),**x))
numbered,details,hits=anarci.anarci([(x['name'],x['ORF']) for x in inputs],scheme='imgt',allow={'A','B'},hmmerpath='/tmp/np-hmmer-bin-20261029/usr/bin',ncpu=2,assign_germline=True,allowed_species=['mouse'],bit_score_threshold=80)
r=dict(inputs=inputs,numbered=numbered,details=details,hits=hits,version=importlib.metadata.version('anarci'),annotation_only=True,reference_used_for_assembly=False,identity_admitted=False)
(b/'raw_classic_numbering.json').write_text(json.dumps(r,indent=2)+'\n')
for x,n,d in zip(inputs,numbered,details):print(x['name'],[(a[1],a[2],len(a[0])) for a in n] if n else None,d)
