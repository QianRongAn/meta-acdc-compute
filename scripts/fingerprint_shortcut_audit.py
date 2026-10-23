"""Retrospective E12/E13/E14 diagnostics; never trains or overwrites inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import torch
from meta_acdc.models.train_fingerprint import auroc
from scripts.rashomon_fingerprint import bootstrap_ci, finite_json, spearman_rho

STRATA = ('all', 'same_family', 'same_peptide', 'same_domain')

def cosine_rows(a):
    norms = np.linalg.norm(a, axis=1, keepdims=True)
    return (a / np.maximum(norms, 1e-12)) @ (a / np.maximum(norms, 1e-12)).T

def query_auc(scores, records, queries, stratum):
    values, pair_scores, pair_labels = [], [], []
    field = {'same_family':'family_id','same_peptide':'peptide','same_domain':'domain'}.get(stratum)
    for q in queries:
        gallery = [j for j in range(len(records)) if j != q and
                   (records[j]['receptor_id'] == records[q]['receptor_id'] or
                    field is None or records[j][field] == records[q][field])]
        labels = np.array([records[j]['receptor_id'] == records[q]['receptor_id'] for j in gallery])
        s = scores[q, gallery]
        valid = labels.any() and (~labels).any()
        values.append(auroc(s, labels.astype(float)) if valid else np.nan)
        if valid:
            pair_scores.extend(s.tolist()); pair_labels.extend(labels.astype(float).tolist())
    vals = np.asarray(values)
    result = bootstrap_ci(vals, [records[q]['family_id'] for q in queries])
    result['per_query_macro_auroc'] = float(np.nanmean(vals)) if np.isfinite(vals).any() else None
    result['pooled_pair_auroc_defined_queries'] = auroc(np.array(pair_scores),np.array(pair_labels)) if pair_labels else None
    result['n_pairs_defined_queries'] = len(pair_labels)
    result['n_undefined_queries'] = int(np.sum(~np.isfinite(vals)))
    return result, vals

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset',type=Path,required=True)
    ap.add_argument('--fold-assignment',type=Path,required=True)
    ap.add_argument('--followup-dir',type=Path,required=True)
    ap.add_argument('--plan',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True)
    args = ap.parse_args()
    if args.out_dir.exists():
        raise FileExistsError(args.out_dir)
    files = [args.dataset,args.fold_assignment,args.plan]
    files += [args.followup_dir/'model_run_v2'/f'sim_fold{f}_seed20261021.npy' for f in range(5)]
    def hashes():
        return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    before = hashes()
    d = torch.load(args.dataset,map_location='cpu',weights_only=False)
    graphs, records = d['graphs'], d['records']
    assignment = json.loads(args.fold_assignment.read_text())
    fold = assignment['per_graph_fold']
    if assignment['graph_records'] != [r['complex_id'] for r in records]:
        raise ValueError('graph identities differ from frozen assignment')
    if len(fold) != len(records):
        raise ValueError('fold/dataset length mismatch')
    comp = np.array([np.asarray(g[0])[:,4:24].mean(0) for g in graphs])
    size = np.array([len(g[0]) for g in graphs],dtype=float)
    radius = np.array([np.linalg.norm(np.asarray(g[1])-np.asarray(g[1]).mean(0),axis=1).mean() for g in graphs])
    domains = sorted(set(r['domain'] for r in records))
    domain = np.eye(len(domains))[[domains.index(r['domain']) for r in records]]
    matrices = dict(composition_cosine=cosine_rows(comp), size_scalar_cosine=cosine_rows(size[:,None]),
                    size_log_distance=-abs(np.log(size[:,None])-np.log(size[None,:])),
                    radius_distance=-abs(radius[:,None]-radius[None,:]),domain_metadata=cosine_rows(domain))
    model = np.zeros((len(records),len(records)))
    for f in range(5):
        s = np.load(files[3+f],allow_pickle=False)
        if s.shape != model.shape or not np.isfinite(s).all():
            raise ValueError('invalid model matrix')
        idx = [i for i,v in enumerate(fold) if v == f]
        model[idx] = s[idx]
    matrices['historical_v2_seed20261021'] = model
    queries = list(range(len(records)))
    results, rows = {}, {}
    for name,s in matrices.items():
        results[name], rows[name] = {}, {}
        for stratum in STRATA:
            results[name][stratum], rows[name][stratum] = query_auc(s,records,queries,stratum)
            per_fold = [query_auc(s, records, [q for q in queries if fold[q] == f], stratum)[0]
                        for f in range(5)]
            defined = [v['pooled_pair_auroc_defined_queries'] for v in per_fold
                       if v['pooled_pair_auroc_defined_queries'] is not None]
            results[name][stratum]['fold_mean_pooled_pair_auroc'] = float(np.mean(defined)) if defined else None
            results[name][stratum]['per_fold_pooled_pair_auroc'] = [v['pooled_pair_auroc_defined_queries'] for v in per_fold]
            results[name][stratum]['cross_fold_pooling_caveat'] = 'Different fold models have uncalibrated score scales; use fold mean or query macro, not cross-fold pooled pair AUROC.'
    paired = {}
    for name in matrices:
        if name == 'historical_v2_seed20261021':
            continue
        paired[name] = {}
        for st in STRATA:
            a,b = rows['historical_v2_seed20261021'][st],rows[name][st]
            paired[name][st] = bootstrap_ci(a-b,[r['family_id'] for r in records])
    rho = np.array([spearman_rho(matrices['composition_cosine'][q, np.arange(len(records))!=q],
                               matrices['composition_cosine'][q, np.arange(len(records))!=q]) for q in queries])
    confidence = {}
    for domain_name in domains:
        hs = [np.asarray(g[0])[:,3] for g,r in zip(graphs,records) if r['domain']==domain_name]
        confidence[domain_name] = dict(n_graphs=len(hs),minimum=float(min(h.min() for h in hs)),
                                      maximum=float(max(h.max() for h in hs)),
                                      variable_graphs=sum(bool(np.ptp(h)>0) for h in hs))
    if before != hashes():
        raise RuntimeError('source inputs changed during analysis')
    result = dict(scope='retrospective defective cached training, not new validation',source_sha256=before,
                  input_hashes_unchanged=True,confidence_audit=confidence,results=results,
                  paired_model_minus_baseline=paired,
                  composition_identical_instance_rho=bootstrap_ci(rho,[r['family_id'] for r in records]),
                  composition_null_scope='deterministic identical copies, NOT shuffled-label training',
                  missing=['audited MHC for same-MHC stratum','shuffled-label trained null','corrected full training E15'])
    args.out_dir.mkdir(parents=True)
    (args.out_dir/'results.json').write_text(json.dumps(finite_json(result),indent=2,allow_nan=False)+'\n')
    lines = ['# Retrospective shortcut diagnostics','',
             'Historical v2 weights used incorrect batch edges. All numbers below diagnose artifacts; none restore a validated gain.',
             '', '| Representation | Mean fold pooled AUROC | Query macro | Family equal (conditional 95% CI) | Same-family macro |',
             '|---|---:|---:|---|---:|']
    def fmt(x): return 'undefined' if x is None else f'{x:.4f}'
    for name,rs in results.items():
        a=rs['all'];h=rs['same_family']
        lines.append(f"| {name} | {fmt(a['fold_mean_pooled_pair_auroc'])} | {fmt(a['per_query_macro_auroc'])} | {fmt(a['mean'])} {a['ci95']} | {fmt(h['per_query_macro_auroc'])} |")
    lines += ['', 'All positives are retained; only negatives are restricted in each stratum. Undefined queries are counted in JSON.',
              'Scalar positive size cosine is degenerate by design; log-size distance is the predeclared sensitivity analysis.',
              'Domain metadata is a confound diagnostic, not evidence the confidence feature varies (it does not in this pool).',
              'Identical composition copies have rho=1 by construction where defined; this is not a trained permutation control.',
              'No audited MHC field: same-MHC/pMHC analysis is unavailable. No fresh binding acceptance data were read.']
    (args.out_dir/'README.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'confidence':confidence,'macro_all':{k:v['all']['per_query_macro_auroc'] for k,v in results.items()}},indent=2))

if __name__ == '__main__':
    main()
