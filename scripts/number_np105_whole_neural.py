"""Independent neural numbering of observed whole-run ORFs; retain every disagreement."""
import sys,json,argparse,hashlib
from pathlib import Path
sys.path.insert(0,'/tmp/np-numbering-deps-20261029')
from anarcii import Anarcii

def occupied_numbering(result):
 return {tuple(p):aa for p,aa in (result.get('numbering') or []) if aa!='-'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--classic',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();classic=json.loads(a.classic.read_text());inputs={x['name']:x for c in classic['cells'] for role in ['TRA','TRB'] for x in c[role]};unique={}
 for x in inputs.values():unique.setdefault(x['ORF'],x['name'])
 model=Anarcii(seq_type='tcr',mode='accuracy',cpu=True,ncpu=2,batch_size=8);result=model.number([(name,seq) for seq,name in unique.items()]);rows=[]
 for name,x in inputs.items():
  n=result[unique[x['ORF']]];old=x.get('numbered');cmap={tuple(p):aa for p,aa in old[0][0] if aa!='-'} if old and len(old)==1 else None;nmap=occupied_numbering(n);rows.append(dict(name=name,barcode=x['barcode'],role=x['role'],classic_complete=x['complete'],neural=n,non_gap_position_map_identical=cmap==nmap,neural_role_matches_source=n.get('chain_type')==('A' if x['role']=='TRA' else 'B'),ORF_sha256=hashlib.sha256(x['ORF'].encode()).hexdigest()))
 record=dict(n_ORFs=len(inputs),n_unique_ORFs=len(unique),n_classic_complete=sum(x['classic_complete'] for x in rows),n_classic_complete_same_neural_mapping=sum(x['classic_complete'] and x['non_gap_position_map_identical'] and x['neural_role_matches_source'] for x in rows),n_classic_complete_disagreements=sum(x['classic_complete'] and not (x['non_gap_position_map_identical'] and x['neural_role_matches_source']) for x in rows),rows=rows,version='2.0.8',annotation_only=True,identity_admitted=False,protocol_frozen=False,independent_acceptance_accessed=False);a.out.write_text(json.dumps(record,indent=2)+'\n');print({k:v for k,v in record.items() if k!='rows'})
if __name__=='__main__':main()
