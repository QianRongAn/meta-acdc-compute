"""A fixed official ERGO model on an already consumed functional association panel."""
from pathlib import Path
import csv,hashlib,io,json,os,pickle,subprocess,sys,tempfile
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'docs/six_blockers_20261021';S=ROOT/'data/external_evidence/20261021-six-blockers'
PANEL=ROOT/'data/external_evidence/20261017-functional-benchmark-v1/functional_observations.json'
PRIOR=ROOT/'data/external_evidence/20261017-functional-evaluation-v1/predictions.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class DataOnlyUnpickler(pickle.Unpickler):
 def find_class(self,module,name):raise ValueError('Executable pickle global rejected')
def main():
 out=D/'ergo_comparison'
 if out.exists():raise FileExistsError(out)
 plan_path=D/'ergo_comparison_plan_v2.json';plan=json.loads(plan_path.read_text())
 for name,h in plan['source_sha256'].items():
  if sha(S/name)!=h:raise ValueError('Official artifact changed: '+name)
 if sha(PANEL)!=plan['panel_sha256'] or sha(PRIOR)!=plan['prior_baseline_predictions_sha256']:raise ValueError('Panel or baseline changed')
 records=json.loads(PANEL.read_text());prior={r['evidence_id']:r for r in json.loads(PRIOR.read_text())}
 if len(records)!=10000 or set(prior)!={r['evidence_id'] for r in records}:raise ValueError('Observation IDs differ')
 out.mkdir();queries=[(r['cdr3_beta'],r['peptide']) for r in records];execution=queries+[queries[0]]
 with (out/'queries.csv').open('w',newline='') as f:csv.writer(f,lineterminator='\n').writerows(execution)
 with tempfile.TemporaryDirectory(prefix='meta-acdc-ergo-official-') as td:
  td=Path(td)
  for name in ['ERGO.py','ERGO_models.py','ae_utils.py','lstm_utils.py','ergo_data_loader.py']:(td/name).write_bytes((S/('ergo_'+name)).read_bytes())
  command=[sys.executable,str(td/'ERGO.py'),'predict','lstm','mcpas','specific','cpu','--model_file='+str(S/'ergo_lstm_mcpas1.pt'),'--test_data_file='+str(out/'queries.csv')]
  env=dict(os.environ,OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
  run=subprocess.run(command,cwd=td,env=env,text=True,capture_output=True)
  (out/'official_stdout.tsv').write_text(run.stdout);(out/'official_stderr.log').write_text(run.stderr)
  run.check_returncode()
 scored=[]
 for line in run.stdout.splitlines():
  fields=line.split('\t')
  if len(fields)!=3:raise ValueError('Unexpected CLI output')
  scored.append((fields[0],fields[1],float(fields[2])))
 if len(scored)!=10001 or [(a,b) for a,b,_ in scored]!=execution:raise ValueError('CLI dropped or reordered queries')
 if not np.isfinite([v for _,_,v in scored]).all() or any(not 0<=v<=1 for _,_,v in scored):raise ValueError('Invalid score')
 if abs(scored[0][2]-scored[-1][2])>1e-5:raise ValueError('Repeated control inconsistent')
 training=DataOnlyUnpickler(io.BytesIO((S/'ergo_lstm_mcpas_train1.pickle').read_bytes())).load()
 # Author file contains (beta, peptide tuple, p/n) records; use identity fields only.
 train_pairs={(r[0],r[1][0] if isinstance(r[1],(list,tuple)) else r[1]) for r in training}
 train_beta={a for a,b in train_pairs};pred=[]
 for r,(_,_,value) in zip(records,scored[:10000]):
  before=prior[r['evidence_id']]
  if before['peptide']!=r['peptide'] or before['observed_activation']!=r['observed_activation']:raise ValueError('Baseline target mismatch')
  pred.append(dict(evidence_id=r['evidence_id'],receptor_pair_id=r['receptor_pair_id'],peptide=r['peptide'],hla=r['mhc_allele'],observed_activation=r['observed_activation'],ergo_score=value,spr_reference_pkd_score=before['spr_reference_pkd_score'],official_training_exact_beta=r['cdr3_beta'] in train_beta,official_training_exact_beta_peptide=(r['cdr3_beta'],r['peptide']) in train_pairs))
 metrics=[]
 for peptide in sorted({r['peptide'] for r in pred}):
  group=[r for r in pred if r['peptide']==peptide];y=[r['observed_activation'] for r in group]
  if set(y)!={0,1}:raise ValueError('One-class panel')
  row=dict(peptide=peptide,n=len(group),n_positive=sum(y))
  for col in ['ergo_score','spr_reference_pkd_score']:
   scores=[r[col] for r in group];row[col]=dict(auroc=float(roc_auc_score(y,scores)),average_precision=float(average_precision_score(y,scores)))
  metrics.append(row)
 report=dict(task=plan['task'],plan_sha256=sha(plan_path),checkpoint_sha256=sha(S/'ergo_lstm_mcpas1.pt'),tool_commit=plan['tool_commit'],n_observations=len(pred),n_receptors=len({r['receptor_pair_id'] for r in pred}),n_peptides=len(metrics),n_execution_controls=1,control_excluded_from_metrics=True,all_input_pairs_verified_in_original_order=True,repeated_control_absolute_difference=abs(scored[0][2]-scored[-1][2]),official_source_modified=False,official_training_exact_beta_observations=sum(r['official_training_exact_beta'] for r in pred),official_training_exact_beta_peptide_observations=sum(r['official_training_exact_beta_peptide'] for r in pred),metrics={c:dict(macro_auroc=float(np.mean([r[c]['auroc'] for r in metrics])),macro_average_precision=float(np.mean([r[c]['average_precision'] for r in metrics]))) for c in ['ergo_score','spr_reference_pkd_score']},panel_previously_consumed=True,model_fit_or_tuning_performed=False,new_independent_acceptance=False,physical_binding_validation_passed=False,clinical_validation_passed=False,limitations=['Common MIRA activation ranking endpoint only; classifiers/SPR scores are not calibrated activation probabilities','No new independent holdout; alpha-chain, family and study overlap for external model not fully excluded','ERGO uses beta CDR3 and peptide; SPR score uses a different training endpoint and more features','No claim that synthetic decoys are real negatives or that a clinical model is validated'])
 for name,value in [('predictions.json',pred),('per_peptide_metrics.json',metrics),('report.json',report)]: (out/name).write_text(json.dumps(value,indent=2)+'\n')
 (out/'manifest.json').write_text(json.dumps({p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},indent=2)+'\n')
 print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
