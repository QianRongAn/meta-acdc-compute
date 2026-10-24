"""Export index-exact legacy metadata and verify cached metrics; no graph arrays exported."""
from pathlib import Path
from collections import defaultdict,Counter
import json,hashlib,argparse
import numpy as np
import torch
from sklearn.metrics import roc_auc_score
from meta_acdc.models.train_egnn import normalize_node_features
from meta_acdc.structure.feature_schema import CONFIDENCE_INDEX
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out-dir',type=Path,required=True);a=ap.parse_args()
 if a.out_dir.exists():raise FileExistsError(a.out_dir)
 a.out_dir.mkdir(parents=True)
 old=ROOT/'docs/fingerprint_retrain_20261023';legacy=ROOT/'docs/fingerprint_retrieval_20261021';path=old/'input/fingerprint_dataset.pt';data=torch.load(path,map_location='cpu',weights_only=False);records=data['records'];n=len(records);fold_path=legacy/'fold_assignment.json';fold=json.loads(fold_path.read_text())['per_graph_fold'];assert len(fold)==n
 protocol=json.loads((old/'protocol.json').read_text());assert sha(path)==protocol['dataset_sha256']
 registry=json.loads((ROOT/'docs/fingerprint_paired_identity_20261024/curation_v3/identity_registry.json').read_text());canonical={x['index']:x for x in registry['rows'] if x['status']=='admitted'}
 exported=[];rid_folds=defaultdict(set);family_folds=defaultdict(set);pdb_rows=defaultdict(list);canonical_folds=defaultdict(set)
 for i,record in enumerate(records):
  assert record['receptor_id']==hashlib.sha256((record['alpha']+'|'+record['beta']).encode()).hexdigest()
  item={k:record.get(k) for k in ['complex_id','pdb','receptor_id','family_id','domain','peptide','source']};item.update(graph_index=i,matrix_row=i,matrix_column=i,fold=fold[i],source_graph_hash=data['graph_hashes'][i],MHC=record.get('MHC',record.get('mhc')),MHC_status='not_recorded_in_legacy_records; no PDB-name inference',identity_semantics='legacy observed full-chain string hash; not validated complete paired variable identity',canonical_identity_status='admitted' if i in canonical else 'quarantined',canonical_paired_variable_id=canonical.get(i,{}).get('receptor_id'))
  if i in canonical:canonical_folds[canonical[i]['receptor_id']].add(fold[i])
  exported.append(item);rid_folds[item['receptor_id']].add(fold[i]);family_folds[item['family_id']].add(fold[i]);pdb_rows[item['pdb']].append(item)
 confidence={}
 for domain in sorted({r['domain'] for r in records}):
  selected=[np.asarray(g[0]) for g,record in zip(data['graphs'],records) if record['domain']==domain];raw=np.concatenate([h[:,CONFIDENCE_INDEX] for h in selected]);norm=np.concatenate([normalize_node_features(torch.as_tensor(h,dtype=torch.float32))[:,CONFIDENCE_INDEX].numpy() for h in selected]);confidence[domain]=dict(raw_unique=np.unique(raw).tolist(),normalized_min=float(norm.min()),normalized_max=float(norm.max()),n_graphs=len(selected))
 binding=[]
 for f in range(5):
  for seed in [20261021,1101,4242]:
   tag=f'fold{f}_seed{seed}';sp=legacy/'followup/model_run_v2'/f'sim_{tag}.npy';ep=sp.with_name(f'emb_{tag}.npy');mp=legacy/'followup'/f'metrics_{tag}.json';S=np.load(sp,allow_pickle=False);E=np.load(ep,allow_pickle=False);assert S.shape==(n,n) and E.shape[0]==n
   scores=[];labels=[]
   for q in range(n):
    if fold[q]!=f:continue
    gallery=[j for j in range(n) if j!=q];truth=[records[q]['receptor_id']==records[j]['receptor_id'] for j in gallery]
    if not any(truth):continue
    scores.extend(S[q,gallery].tolist());labels.extend(truth)
   auc=float(roc_auc_score(labels,scores));reported=json.loads(mp.read_text())['retrieval']['auroc'];assert abs(auc-reported)<1e-12,(tag,auc,reported)
   binding.append(dict(instance=tag,similarity_file=sp.relative_to(ROOT).as_posix(),similarity_sha256=sha(sp),embedding_file=ep.relative_to(ROOT).as_posix(),embedding_sha256=sha(ep),n_graphs=n,recomputed_original_pooled_pair_AUROC=auc,reported_original_AUROC=reported,exact_original_metric_replay=True))
 mixed=[dict(pdb=p,rows=rows) for p,rows in pdb_rows.items() if len({x['fold'] for x in rows})>1]
 result=dict(n_graphs=n,n_observed_sequence_keys=len(rid_folds),n_legacy_families=len(family_folds),observed_keys_with_ge2=sum(c>=2 for c in Counter(x['receptor_id'] for x in records).values()),n_legacy_observed_keys_crossing_folds=sum(len(v)>1 for v in rid_folds.values()),n_legacy_families_crossing_folds=sum(len(v)>1 for v in family_folds.values()),n_pdb_names_crossing_folds=len(mixed),n_graphs_in_cross_fold_pdbs=sum(len(x['rows']) for x in mixed),n_admitted_canonical_ids_crossing_legacy_folds=sum(len(v)>1 for v in canonical_folds.values()),n_canonical_identities_checked=len(canonical_folds),n_original_instances_metric_replayed=len(binding),seed_20261021_mean_fold_AUROC=float(np.mean([x['reported_original_AUROC'] for x in binding if x['instance'].endswith('seed20261021')])),dataset_sha256=sha(path),fold_assignment_sha256=sha(fold_path),binding_scope='retrospective index binding replayed against original metrics; no contemporaneous per-instance dataset hash in legacy metric files',independent_acceptance_accessed=False,no_graph_arrays_exported=True,confidence_column_by_domain=confidence)
 for name,obj in [('records_index_exact.json',dict(summary=result,records=exported)),('legacy_matrix_binding.json',dict(summary=result,instances=binding)),('cross_fold_cases.json',dict(summary=result,pdb_cases=mixed))]:(a.out_dir/name).write_text(json.dumps(obj,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
