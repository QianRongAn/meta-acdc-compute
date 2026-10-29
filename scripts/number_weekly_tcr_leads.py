"""Independent numbering of two newly released public TCR leads; discovery only."""
import sys,json,urllib.parse
from pathlib import Path
sys.path[:0]=['.','src','/tmp/np-classic-deps-20261029']
import anarci
from scripts.audit_np_boundaries import checked_numbering
B=Path('docs/fingerprint_np105_validation_20261029')
def main():
 es=json.loads((B/'sources/weekly_metadata/weekly_entities.json').read_text())['data']['entries'];inputs=[]
 for e in es:
  if e['rcsb_id'] not in ['9XC3','9X9I']:continue
  for entity in e['polymer_entities']:
   desc=entity['rcsb_polymer_entity']['pdbx_description']
   if not desc.startswith('TCR1414_'):continue
   inputs.append(dict(name=entity['rcsb_id'],pdb=e['rcsb_id'],role='alpha' if desc.endswith('alpha') else 'beta',sequence=entity['entity_poly']['pdbx_seq_one_letter_code_can'],description=desc))
 n,d,h=anarci.anarci([(x['name'],x['sequence']) for x in inputs],scheme='imgt',allow={'A','B'},hmmerpath='/tmp/np-hmmer-bin-20261029/usr/bin',ncpu=2,assign_germline=True,allowed_species=['human'],bit_score_threshold=80)
 rows=[];requests=[]
 for x,nn,dd,hh in zip(inputs,n,d,h):
  row=dict(**x,numbered=nn,details=dd,hits=hh,identity_admitted=False)
  try:
   assert nn and len(nn)==1 and dd[0]['species']=='human';mapping,start,end=nn[0];boundary=checked_numbering(x['sequence'],mapping,start,end,dd[0]['chain_type'],'A' if x['role']=='alpha' else 'B');row['boundary']=boundary
   q=dict(query=dict(type='terminal',service='sequence',parameters=dict(evalue_cutoff=.1,identity_cutoff=1.,sequence_type='protein',value=boundary['sequence'])),return_type='polymer_entity',request_options=dict(return_all_hits=True));requests.append([x['name']+'.json','https://search.rcsb.org/rcsbsearch/v2/query?json='+urllib.parse.quote(json.dumps(q,separators=(',',':')))])
  except (AssertionError,ValueError,KeyError) as e:row['unresolved_reason']=str(e) or type(e).__name__
  rows.append(row)
 (B/'weekly_numbering.json').write_text(json.dumps(dict(rows=rows,identity_admitted=False,protocol_frozen=False),indent=2)+'\n');(B/'weekly_chain_requests.json').write_text(json.dumps(requests,indent=2)+'\n');print([(x['name'],x.get('boundary',{}).get('length'),x.get('details'),x.get('unresolved_reason')) for x in rows])
if __name__=='__main__':main()
