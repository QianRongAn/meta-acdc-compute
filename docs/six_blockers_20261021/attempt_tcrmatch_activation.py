"""Locked, retrospective TCRMatch comparison on the already consumed MIRA panel.

This tests activation ranking, not SPR affinity or clinical risk. No model fitting.
"""
from pathlib import Path
import collections,csv,difflib,hashlib,io,json,subprocess,tempfile
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score
ROOT=Path(__file__).resolve().parents[2]
D=ROOT/'docs/six_blockers_20261021';SOURCE=ROOT/'data/external_evidence/20261021-six-blockers'
RAW=Path('/home/administrator/docs/meta-acdc/data/raw/vdjdb/vdjdb-2026-06-03/vdjdb.txt')
PANEL=ROOT/'data/external_evidence/20261017-functional-benchmark-v1/functional_observations.json'
PRIOR=ROOT/'data/external_evidence/20261017-functional-evaluation-v1/predictions.json'
AA=set('ACDEFGHIKLMNPQRSTVWY')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def trim(s):return s[1:-1] if s.startswith('C') and s[-1:] in ('F','W') else s
def near(a,b):
 if min(len(a),len(b))/max(len(a),len(b))<.9:return False
 m=difflib.SequenceMatcher(None,a,b,autojunk=False)
 return m.quick_ratio()>=.9 and m.ratio()>=.9
def eligible(row,peptides):
 s=row['cdr3']
 return row['gene']=='TRB' and row['species']=='HomoSapiens' and int(row['vdjdb.score'])>=2 and row['antigen.epitope'] in peptides and len(s)>=6 and set(s)<=AA

def main():
 out=D/'tcrmatch_comparison'
 if out.exists():raise FileExistsError(out)
 planpath=D/'tcrmatch_comparison_plan.json';plan=json.loads(planpath.read_text())
 for p,key in [(RAW,'reference_sha256'),(PANEL,'panel_sha256'),(PRIOR,'prior_baseline_predictions_sha256')]:
  if sha(p)!=plan[key]:raise ValueError('Locked input changed: '+str(p))
 records=json.loads(PANEL.read_text());prior={r['evidence_id']:r for r in json.loads(PRIOR.read_text())}
 if set(prior)!={r['evidence_id'] for r in records}:raise ValueError('Baseline evidence IDs differ')
 queries=sorted({r['cdr3_beta'] for r in records});peptides={r['peptide'] for r in records}
 if any(not q or set(q)-AA for q in queries):raise ValueError('Noncanonical query')
 with RAW.open() as f:references=[r for r in csv.DictReader(f,delimiter='\t') if eligible(r,peptides)]
 unique={r['cdr3'] for r in references};exclude={s for s in unique if any(near(s,q) for q in queries)}
 references=[r for r in references if r['cdr3'] not in exclude]
 reference_rows=sorted({(trim(r['cdr3']),r['cdr3'],r['reference.id'].replace('\t',' '),r['antigen.epitope'],r['antigen.species'].replace(' ','_'),r['antigen.gene'].replace(' ','_') or 'unknown') for r in references})
 if not reference_rows:raise ValueError('No eligible reference remains')
 print('reference rows',len(reference_rows),'query CDR3s',len(queries),'excluded',len(exclude),flush=True)
 out.mkdir()
 with (out/'reference.tsv').open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['trimmed_seq','original_seq','receptor_group','epitope','organism','antigen']);w.writerows(reference_rows)
 (out/'queries.txt').write_text('\n'.join(queries)+'\n')
 with tempfile.TemporaryDirectory(prefix='meta-acdc-tcrmatch-') as td:
  td=Path(td)
  for src,dst in [('tcrmatch_main.cpp','main.cpp'),('tcrmatch.cpp','tcrmatch.cpp')]: (td/dst).write_bytes((SOURCE/src).read_bytes())
  binary=td/'tcrmatch';command=['g++','-std=c++11','-fopenmp','-O3','-g','-o',str(binary),str(td/'main.cpp')]
  compile_run=subprocess.run(command,text=True,capture_output=True,check=True)
  # Check the unmodified executable returns an exact self-match before evaluating.
  (td/'control.txt').write_text('CASSIRSSYEQYF\n');(td/'control.tsv').write_text('trimmed_seq\toriginal_seq\treceptor_group\tepitope\torganism\tantigen\nASSIRSSYEQY\tCASSIRSSYEQYF\tcontrol\tPEPTIDE\thuman\tcontrol\n')
  control=subprocess.run([str(binary),'-i',str(td/'control.txt'),'-d',str(td/'control.tsv'),'-t','1','-s','.97'],text=True,capture_output=True,check=True)
  selfrows=list(csv.DictReader(io.StringIO(control.stdout),delimiter='\t'))
  if len(selfrows)!=1 or float(selfrows[0]['score'])<.999:raise ValueError('TCRMatch self-test failed')
  command=[str(binary),'-i',str(out/'queries.txt'),'-d',str(out/'reference.tsv'),'-t','2','-s',str(plan['threshold'])]
  run=subprocess.run(command,text=True,capture_output=True,check=True)
  if run.stderr.strip():raise ValueError('TCRMatch reported a warning or incomplete execution: '+run.stderr)
  (out/'raw_matches.tsv').write_text(run.stdout);(out/'execution.stderr').write_text(compile_run.stderr+run.stderr)
  binary_sha=sha(binary);compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0]
 best={};trimmed_queries={trim(q) for q in queries}
 for row in csv.DictReader(io.StringIO(run.stdout),delimiter='\t'):
  q=row['trimmed_input_sequence'];pep=row['epitope'];score=float(row['score'])
  if q not in trimmed_queries or pep not in peptides or not 0<=score<=1.00001:raise ValueError('Unexpected tool output')
  best[q,pep]=max(best.get((q,pep),0),score)
 predictions=[]
 for r in records:
  old=prior[r['evidence_id']]
  if old['peptide']!=r['peptide'] or old['observed_activation']!=r['observed_activation']:raise ValueError('Baseline target mismatch')
  predictions.append(dict(evidence_id=r['evidence_id'],peptide=r['peptide'],hla=r['mhc_allele'],receptor_pair_id=r['receptor_pair_id'],observed_activation=r['observed_activation'],tcrmatch_score=best.get((trim(r['cdr3_beta']),r['peptide']),0),spr_reference_pkd_score=old['spr_reference_pkd_score']))
 metrics=[]
 for peptide in sorted(peptides):
  group=[r for r in predictions if r['peptide']==peptide];y=[r['observed_activation'] for r in group]
  if set(y)!={0,1}:raise ValueError('One-class peptide panel')
  row=dict(peptide=peptide,n=len(group),n_positive=sum(y),n_tcrmatch_nonzero=sum(r['tcrmatch_score']>0 for r in group))
  for col in ['tcrmatch_score','spr_reference_pkd_score']:
   scores=[r[col] for r in group];row[col]=dict(auroc=float(roc_auc_score(y,scores)),average_precision=float(average_precision_score(y,scores)))
  metrics.append(row)
 report=dict(task=plan['task'],plan_sha256=sha(planpath),n_observations=len(records),n_receptors=len({r['receptor_pair_id'] for r in records}),n_peptides=len(peptides),reference_snapshot_sha256=sha(RAW),reference_rows=len(reference_rows),excluded_near_query_reference_cdr3s=len(exclude),query_cdr3s=len(queries),tool_commit=plan['tool_commit'],tool_source_sha256={p:sha(SOURCE/p) for p in ['tcrmatch_main.cpp','tcrmatch.cpp','tcrmatch_LICENSE']},compiler=compiler,binary_sha256=binary_sha,self_match_control_passed=True,no_hit_pairs=sum(r['tcrmatch_score']==0 for r in predictions),metrics={col:dict(macro_auroc=float(np.mean([r[col]['auroc'] for r in metrics])),macro_average_precision=float(np.mean([r[col]['average_precision'] for r in metrics]))) for col in ['tcrmatch_score','spr_reference_pkd_score']},panel_previously_consumed=True,new_independent_acceptance=False,physical_binding_validation_passed=False,clinical_validation_passed=False,model_fit_or_tuning_performed=False,limitations=['Retrospective common MIRA activation ranking endpoint; neither score is a calibrated activation probability','CDR3beta-only reference exclusion does not certify full receptor family or study independence','TCRMatch uses only beta similarity/reference peptide labels; the SPR model was trained on another endpoint','No-hit values are zero ties, not experimentally proven nonbinding','No clinical outcomes or fresh holdout consumed'])
 for name,value in [('predictions.json',predictions),('per_peptide_metrics.json',metrics),('report.json',report)]: (out/name).write_text(json.dumps(value,indent=2)+'\n')
 (out/'manifest.json').write_text(json.dumps({p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},indent=2)+'\n')
 print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
