"""Test one historically documented candidate; failure is retained as evidence."""
from pathlib import Path
import subprocess,tarfile,io,sys,json,hashlib,csv,tempfile
output=Path(sys.argv[1])
if output.exists():raise FileExistsError(output)
root=Path('/home/administrator/docs/meta-acdc-readiness');temporary=tempfile.TemporaryDirectory(prefix='meta-acdc-21-forensic-');out=Path(temporary.name)
b=subprocess.check_output(['git','-C',str(root),'archive','6797863','src'])
with tarfile.open(fileobj=io.BytesIO(b)) as tf:
 for m in tf:
  if m.isfile() and m.name.startswith('src/'):
   p=out/m.name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(tf.extractfile(m).read())
sys.path.insert(0,str(out/'src'))
import torch
from meta_acdc.models.egnn import EGNN
from meta_acdc.structure.cif import parse_cif
from meta_acdc.structure.graph import build_interface_graph_from_residues
from meta_acdc.structure.score_predictions import score_graph
torch.set_num_threads(2)
p=Path('/home/administrator/docs/meta-acdc/data/processed/egnn_dataset.model.pt');c=torch.load(p,map_location='cpu',weights_only=False)
net=EGNN(node_dim=c['n_node_dim'],edge_dim=c['n_edge_dim'],depth=6,hidden=128);net.load_state_dict(c['model']);net.eval()
with (root/'docs/delivery_audit_20261018/historical/data/processed/prediction_scores.relabeled.tsv').open() as f:
 old={r['job_id']:r for r in csv.DictReader(f,delimiter='\t')}
results=[]
for job in ['5brz_EVDPIGHLY','5brz_ESDPIVAQY','5brz_ILAKFLHWL']:
 path=Path('/home/administrator/docs/meta-acdc/data/raw/af3_predictions')/(job+'_model_0.cif')
 g=build_interface_graph_from_residues(parse_cif(path));value=score_graph(net,g,'cpu')
 results.append(dict(job_id=path.stem,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),replay_score=value,historical_score=old[path.stem]['score'],matches_legacy_four_decimals=f'{value:.4f}'==old[path.stem]['score']))

receipt=dict(record_status='FORENSIC_WITHDRAWN_FAILED_CANDIDATE',clinical_interpretation_allowed=False,historical_causal_attestation=False,source_commit=subprocess.check_output(['git','-C',str(root),'rev-parse','6797863'],text=True).strip(),source_archive_sha256=hashlib.sha256(b).hexdigest(),checkpoint_path=str(p),checkpoint_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),torch_version=torch.__version__,device='cpu',mask_plddt=False,temperature=1.0,selection='Historical scorer default checkpoint, three prespecified disputed jobs model_0; bounded diagnostic, not a model selection search',results=results,matches=sum(r['matches_legacy_four_decimals'] for r in results),comparisons=len(results))
output.write_text(json.dumps(receipt,indent=2)+'\n')
temporary.cleanup()
