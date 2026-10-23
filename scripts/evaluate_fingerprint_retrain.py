"""Frozen endpoint analysis of corrected training, baselines, and permutation control."""
import hashlib
import json
from pathlib import Path
from itertools import combinations
import numpy as np
import torch
from scripts.fingerprint_shortcut_audit import query_auc, cosine_rows
from scripts.rashomon_fingerprint import bootstrap_ci, finite_json, spearman_rho

ROOT=Path(__file__).resolve().parents[1]
STRATA=('all','same_family','same_peptide','same_domain')

def evaluate_run(protocol, out):
    result_path=out/'evaluation.json'
    if result_path.exists():
        raise FileExistsError(result_path)
    data=torch.load(protocol['dataset'],map_location='cpu',weights_only=False)
    graphs,records=data['graphs'],data['records']
    assignment=json.loads((ROOT/protocol['fold_assignment']).read_text())
    fold_of=assignment['per_graph_fold']
    if assignment['graph_records'] != [r['complex_id'] for r in records]:
        raise ValueError('graph identities do not match protocol')
    n=len(records);queries=list(range(n));groups=[r['family_id'] for r in records]
    comp=np.array([np.asarray(g[0])[:,4:24].mean(0) for g in graphs])
    size=np.array([len(g[0]) for g in graphs],dtype=float)
    radius=np.array([np.linalg.norm(np.asarray(g[1])-np.asarray(g[1]).mean(0),axis=1).mean() for g in graphs])
    baselines={'composition':cosine_rows(comp),'size_cosine_constant':cosine_rows(size[:,None]),
               'size_log_distance':-abs(np.log(size[:,None])-np.log(size[None,:])),
               'radius_distance':-abs(radius[:,None]-radius[None,:])}
    baseline_results,baseline_rows={},{}
    for name,S in baselines.items():
        baseline_results[name]={};baseline_rows[name]={}
        for st in STRATA:
            baseline_results[name][st],baseline_rows[name][st]=query_auc(S,records,queries,st)
    all_rows={};model_results={};source_hashes={};seed_fold_metrics={};seed_fold_sim={};permutation=[]
    for arm in protocol['arms']:
        model_results[arm]={};all_rows[arm]={};seed_fold_metrics[arm]={};seed_fold_sim[arm]={}
        for seed in protocol['seeds']:
            assembled=np.zeros((n,n));fragments=[];seed_fold_sim[arm][seed]={}
            for fold in range(protocol['n_folds']):
                tag=f'fold{fold}_seed{seed}'
                base=out/arm
                metric_path=base/f'metrics_{tag}.json'
                metric=json.loads(metric_path.read_text());fragments.append(metric)
                for rel,expected in metric['artifact_sha256'].items():
                    p=base/rel;actual=hashlib.sha256(p.read_bytes()).hexdigest()
                    if actual != expected:
                        raise ValueError('artifact binding mismatch: '+str(p))
                    source_hashes[str(p.relative_to(out))]=actual
                source_hashes[str(metric_path.relative_to(out))]=hashlib.sha256(metric_path.read_bytes()).hexdigest()
                S=np.load(base/'model_run_v2'/f'sim_{tag}.npy',allow_pickle=False)
                E=np.load(base/'model_run_v2'/f'emb_{tag}.npy',allow_pickle=False).astype(float)
                if not np.allclose(S,E@E.T,atol=2e-6,rtol=0):
                    raise ValueError('similarity/embedding mismatch')
                idx=[i for i,f in enumerate(fold_of) if f==fold]
                assembled[idx]=S[idx];seed_fold_sim[arm][seed][fold]=S
                lm=json.loads((base/'model_run_v2'/f'labels_{tag}.json').read_text())
                if any(m['target_graph'] in idx or m['source_graph'] in idx for m in lm['mapping']):
                    raise ValueError('permutation touched test fold')
                permutation.append({'arm':arm,'fold':fold,'seed':seed,'n_train':len(lm['mapping']),
                                    'n_changed_receptor_labels':lm['n_changed_receptor_labels']})
            model_results[arm][str(seed)]={};all_rows[arm][seed]={}
            for st in STRATA:
                model_results[arm][str(seed)][st],all_rows[arm][seed][st]=query_auc(assembled,records,queries,st)
            seed_fold_metrics[arm][str(seed)]={'fold_mean_pooled_auroc':float(np.mean([m['retrieval']['auroc'] for m in fragments])),
                                              'fold_mean_rank1':float(np.mean([m['retrieval']['rank1'] for m in fragments])),
                                              'fold_cosine_medians':[m['median_pairwise_cosine'] for m in fragments]}
    summaries={};contrasts={}
    for arm in protocol['arms']:
        summaries[arm]={};contrasts[arm]={}
        for st in STRATA:
            a=np.stack([all_rows[arm][seed][st] for seed in protocol['seeds']],axis=1)
            summaries[arm][st]=bootstrap_ci(a,groups)
            contrasts[arm][st]={name:bootstrap_ci(a-baseline_rows[name][st][:,None],groups) for name in baselines}
    true_minus_null={st:bootstrap_ci(np.stack([all_rows['true'][seed][st]-all_rows['shuffled'][seed][st]
                                for seed in protocol['seeds']],axis=1),groups) for st in STRATA}
    stability={}
    for arm in protocol['arms']:
        stability[arm]={}
        for st in ('all','same_family'):
            vals=[]
            for q in queries:
                f=fold_of[q]
                gallery=[j for j in queries if j!=q and (st=='all' or
                         records[j]['receptor_id']==records[q]['receptor_id'] or records[j]['family_id']==records[q]['family_id'])]
                valid=st=='all' or np.isfinite(all_rows[arm][protocol['seeds'][0]][st][q])
                vals.append([spearman_rho(seed_fold_sim[arm][a][f][q,gallery],seed_fold_sim[arm][b][f][q,gallery])
                             if valid else np.nan for a,b in combinations(protocol['seeds'],2)])
            stability[arm][st]=bootstrap_ci(np.asarray(vals),groups)
    primary=contrasts['true']['same_family']['composition']
    individual=[bootstrap_ci(all_rows['true'][s]['same_family']-baseline_rows['composition']['same_family'],groups)['mean']
                for s in protocol['seeds']]
    positive=(primary['ci95'][0] is not None and primary['ci95'][0]>0 and
              all(v is not None and v>0 for v in individual) and true_minus_null['same_family']['ci95'][0] is not None
              and true_minus_null['same_family']['ci95'][0]>0)
    result=dict(scope='previously exposed developmental identity-retrieval panel; no independent binding acceptance',
                protocol=protocol,primary=primary,primary_per_seed_family_mean=individual,
                preregistered_developmental_evidence_criterion_met=bool(positive),
                practical_or_biological_gain_claim_allowed=False,
                baselines=baseline_results,models=model_results,seed_fold_metrics=seed_fold_metrics,
                family_equal_seed_averaged=summaries,model_minus_baselines=contrasts,
                true_minus_shuffled=true_minus_null,ranking_stability=stability,permutation_checks=permutation,
                artifact_sha256=source_hashes,replay=json.loads((out/'replay_check.json').read_text()),
                uncomputed=['same MHC/pMHC: audited metadata unavailable','new independent binding/clinical acceptance'],
                uncertainty='conditional family-cluster bootstrap, seeds averaged within query; does not resample model training or sources')
    result_path.write_text(json.dumps(finite_json(result),indent=2,allow_nan=False)+'\n')
    lines=['# Corrected fingerprint training: frozen developmental experiment','',
           'This panel was exposed previously. No hyperparameter tuning, no early stopping, no best-seed selection; no independent measured acceptance set was read.',
           '', '| Representation | All-negative family-equal AUROC | Same-family family-equal AUROC |',
           '|---|---:|---:|']
    def fmt(x):return 'undefined' if x is None else f'{x:.4f}'
    for name,res in baseline_results.items():
        lines.append(f"| {name} | {fmt(res['all']['mean'])} | {fmt(res['same_family']['mean'])} |")
    for name,res in summaries.items():
        lines.append(f"| {name} (3 seeds averaged per query) | {fmt(res['all']['mean'])} | {fmt(res['same_family']['mean'])} |")
    lines+=['',f"Primary corrected model minus composition (same-family): {primary['mean']}, conditional CI95 {primary['ci95']}.",
            f"Positive developmental evidence criterion met: {positive}. No claim of practical or biological improvement is allowed.",
            f"Full replay: {result['replay']}",
            'All checkpoints, training-only label mappings, embeddings, similarities, metric fragments and logs are retained; evaluation.json binds hashes.',
            'Only 3 label permutations per fold: a descriptive control, not a calibrated permutation p-value. Same-MHC analysis unavailable.']
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'primary':primary,'criterion_met':positive},indent=2),flush=True)
