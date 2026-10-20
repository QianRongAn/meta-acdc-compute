"""Compare retained per-structure raw tables with the bounded replay."""
from pathlib import Path
import json,hashlib,csv,io
r=Path('/home/administrator/docs/meta-acdc-readiness');d=r/'docs/forensic_followup_20261021';out=d/'recovered_per_structure';out.mkdir(exist_ok=True)
if (d/'per_structure_reconciliation.json').exists():raise FileExistsError('Retain existing evidence; choose a new audit directory for a new run')
full=json.loads((d/'legacy_full_replay.json').read_text());rows=[];invent={};mismatch=[];n=0
for seed in range(5):
 p=Path('/home/administrator/docs/meta-acdc/data/processed')/f'da_scores_seed{seed}.tsv';raw=p.read_bytes();dest=out/p.name
 if dest.exists() and dest.read_bytes()!=raw:raise ValueError('Refusing overwrite')
 dest.write_bytes(raw);table=list(csv.DictReader(io.StringIO(raw.decode()),delimiter='\t'));by={r['job_id']:r for r in table}
 assert len(by)==len(table)
 invent[p.name]=dict(local_source=str(p),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),rows=len(table),record_status='WITHDRAWN_HISTORICAL_RAW')
 for row in full['results']:
  if row['seed']!=seed:continue
  paths=sorted(k for k in full['input_sha256'] if Path(k).stem.startswith(row['job_id']+'_model_'))
  assert len(paths)==len(row['per_af3_model_scores'])
  for path,value in zip(paths,row['per_af3_model_scores']):
   original=by[Path(path).stem];rounded=f'{value:.4f}';n+=1
   if rounded!=original['score']:mismatch.append(dict(job_model=Path(path).stem,seed=seed,replay=rounded,historical=original['score'],absolute_difference=abs(float(rounded)-float(original['score']))))
summary=dict(record_status='FORENSIC_WITHDRAWN',clinical_interpretation_allowed=False,tables=invent,comparisons=n,matches=n-len(mismatch),mismatches=mismatch,max_rounded_absolute_difference=max([v['absolute_difference'] for v in mismatch] or [0]),uncovered_historical_rows=sum(v['rows'] for v in invent.values())-n)
(d/'per_structure_reconciliation.json').write_text(json.dumps(summary,indent=2)+'\n');print(n,summary['matches'],len(mismatch),summary['max_rounded_absolute_difference'],summary['uncovered_historical_rows'])
