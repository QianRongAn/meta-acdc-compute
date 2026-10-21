"""Finite fingerprint search against all retained eligible historical checkpoints."""
from pathlib import Path
import csv,hashlib,io,json,sys,subprocess,tarfile,tempfile
ROOT=Path(__file__).resolve().parents[2];D=Path(__file__).resolve().parent
out=D/'legacy_21_search_results.json'
if out.exists():raise FileExistsError(out)
plan_path=D/'legacy_21_search_plan.json';plan=json.loads(plan_path.read_text())
with (ROOT/'docs/delivery_audit_20261018/historical/data/processed/prediction_scores.relabeled.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
original=Path('/home/administrator/docs/meta-acdc');inputs={};results=[]
with tempfile.TemporaryDirectory(prefix='meta-acdc-legacy-fingerprint-') as td:
 archive=subprocess.check_output(['git','-C',str(ROOT),'archive',plan['source_commit'],'src'])
 with tarfile.open(fileobj=io.BytesIO(archive)) as tf:
  for member in tf:
   if not member.isfile():continue
   if not member.name.startswith('src/') or '..' in Path(member.name).parts:raise ValueError(member.name)
   path=Path(td)/member.name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(tf.extractfile(member).read())
 sys.path.insert(0,str(Path(td)/'src'))
 import torch
 from meta_acdc.models.egnn import EGNN
 from meta_acdc.structure.cif import parse_cif
 from meta_acdc.structure.graph import build_interface_graph_from_residues
 from meta_acdc.structure.score_predictions import score_graph
 torch.set_num_threads(2);graphs={}
 for row in rows:
  p=original/'data/raw/af3_predictions'/(row['job_id']+'.cif')
  if not p.exists():graphs[row['job_id']]=None;continue
  inputs[row['job_id']]=hashlib.sha256(p.read_bytes()).hexdigest()
  graphs[row['job_id']]=build_interface_graph_from_residues(parse_cif(p))
 for candidate in plan['eligible_candidates']:
  p=original/'data/processed'/candidate['file']
  if hashlib.sha256(p.read_bytes()).hexdigest()!=candidate['sha256']:raise ValueError('checkpoint changed')
  c=torch.load(p,map_location='cpu',weights_only=False);net=EGNN(node_dim=c['n_node_dim'],edge_dim=c['n_edge_dim'],depth=6,hidden=128)
  net.load_state_dict(c['model']);net.eval();pred=[]
  for row in rows:
   g=graphs[row['job_id']]
   if g is None:pred.append(dict(job_model=row['job_id'],status='missing_input'));continue
   value=score_graph(net,g,'cpu');score=f'{value:.4f}'
   pred.append(dict(job_model=row['job_id'],replay=score,historical=row['score'],score_matches=score==row['score'],graph_counts_match=(g.n_nodes==int(row['n_nodes']) and len(g.edge_index)==int(row['n_edges']))))
  results.append(dict(checkpoint=candidate['file'],checkpoint_sha256=candidate['sha256'],matches=sum(x.get('score_matches',False) for x in pred),graph_counts_match=sum(x.get('graph_counts_match',False) for x in pred),comparisons=len(rows),predictions=pred))
  print(candidate['file'],results[-1]['matches'],'/',len(rows),flush=True)
 receipt=dict(record_status='FORENSIC_WITHDRAWN',plan_sha256=hashlib.sha256(plan_path.read_bytes()).hexdigest(),input_sha256=inputs,source_archive_sha256=hashlib.sha256(archive).hexdigest(),torch_version=torch.__version__,device='cpu',results=results,complete_matches=[x['checkpoint'] for x in results if x['matches']==len(rows) and x['graph_counts_match']==len(rows)],historical_causal_attestation=False,scientific_model_promotion_allowed=False)
 out.write_text(json.dumps(receipt,indent=2)+'\n')
