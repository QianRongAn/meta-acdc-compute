import sys,json,importlib.metadata
from pathlib import Path
sys.path.insert(0,'/tmp/np-numbering-deps-20261029')
from anarcii import Anarcii
b=Path('docs/fingerprint_np105_raw_20261029');r=json.loads((b/'raw_classic_numbering.json').read_text());inputs=r['inputs'];f=b/'raw_ORFs.fasta';f.write_text(''.join('>'+x['name']+'\n'+x['ORF']+'\n' for x in inputs))
m=Anarcii(seq_type='tcr',mode='accuracy',cpu=True,ncpu=2,batch_size=8);n=m.number(str(f));(b/'raw_neural_numbering.json').write_text(json.dumps(dict(inputs=inputs,results=n,version=importlib.metadata.version('anarcii'),annotation_only=True),indent=2)+'\n')
print({k:dict(chain=v.get('chain_type'),start=v.get('query_start'),end=v.get('query_end'),error=v.get('error')) for k,v in n.items()})
