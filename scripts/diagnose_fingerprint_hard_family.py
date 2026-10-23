"""Retrospective error analysis only: exposed identity panel and old frozen outputs."""
import argparse
from collections import Counter,defaultdict
import hashlib
import json
from pathlib import Path
import random
import numpy as np
import torch
from meta_acdc.models.train_egnn import normalize_node_features
from scripts.fingerprint_shortcut_audit import cosine_rows,query_auc
from scripts.rashomon_fingerprint import finite_json

ROOT=Path(__file__).resolve().parents[1]

def edit_distance(a,b):
    row=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        new=[i]
        for j,y in enumerate(b,1):new.append(min(new[-1]+1,row[j]+1,row[j-1]+(x!=y)))
        row=new
    return row[-1]

def coverage(records,order,batch_size=64):
    counts=Counter()
    for s in range(0,len(order),batch_size):
        idx=order[s:s+batch_size]
        for q in idx:
            pos=[j for j in idx if j!=q and records[j]['receptor_id']==records[q]['receptor_id']]
            hard=[j for j in idx if records[j]['family_id']==records[q]['family_id'] and records[j]['receptor_id']!=records[q]['receptor_id']]
            counts['anchors']+=1;counts['positive_pairs']+=len(pos);counts['hard_negative_pairs']+=len(hard)
            counts['anchors_with_positive']+=bool(pos)
            counts['active_anchors_with_hard_negative']+=bool(pos and hard)
    return dict(counts)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out-dir',type=Path,required=True);a=ap.parse_args()
    if a.out_dir.exists():raise FileExistsError(a.out_dir)
    base=ROOT/'docs/fingerprint_retrain_20261023'
    data_path=base/'input/fingerprint_dataset.pt';d=torch.load(data_path,map_location='cpu',weights_only=False)
    graphs,records=d['graphs'],d['records'];n=len(records)
    folds=json.loads((base/'fold_assignment.json').read_text())['per_graph_fold']
    comp=np.array([np.asarray(g[0])[:,4:24].mean(0) for g in graphs]);Scomp=cosine_rows(comp)
    _,comp_auc=query_auc(Scomp,records,list(range(n)),'same_family')
    seed_rows=[];scores=[];source={str(data_path.relative_to(ROOT)):hashlib.sha256(data_path.read_bytes()).hexdigest()}
    for seed in [20261021,1101,4242]:
        S=np.zeros((n,n))
        for f in range(5):
            p=base/'run/true/model_run_v2'/f'sim_fold{f}_seed{seed}.npy'
            arr=np.load(p,allow_pickle=False);idx=[i for i,v in enumerate(folds) if v==f];S[idx]=arr[idx]
            source[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
        _,row=query_auc(S,records,list(range(n)),'same_family');seed_rows.append(row);scores.append(S)
    model_auc=np.nanmean(np.stack(seed_rows),axis=0,where=np.isfinite(np.stack(seed_rows)))
    queries=[];family=defaultdict(list)
    for q in range(n):
        if not np.isfinite(comp_auc[q]):continue
        pos=[j for j in range(n) if j!=q and records[j]['receptor_id']==records[q]['receptor_id']]
        neg=[j for j in range(n) if records[j]['family_id']==records[q]['family_id'] and records[j]['receptor_id']!=records[q]['receptor_id']]
        meanS=np.mean([S[q] for S in scores],axis=0)
        wrong=max(neg,key=lambda j:meanS[j]);bestpos=max(pos,key=lambda j:meanS[j])
        row=dict(query_index=q,complex_id=records[q]['complex_id'],family_id=records[q]['family_id'],
                 receptor_id=records[q]['receptor_id'],domain=records[q]['domain'],peptide=records[q]['peptide'],
                 n_positive=len(pos),n_hard_negative=len(neg),composition_auc=float(comp_auc[q]),
                 model_auc=float(model_auc[q]),model_minus_composition=float(model_auc[q]-comp_auc[q]),
                 per_seed_auc=[float(x[q]) for x in seed_rows],
                 highest_scoring_wrong_index=wrong,highest_scoring_wrong_complex=records[wrong]['complex_id'],
                 highest_scoring_positive_index=bestpos,
                 wrong_alpha_edit=edit_distance(records[q]['alpha'],records[wrong]['alpha']),
                 wrong_beta_edit=edit_distance(records[q]['beta'],records[wrong]['beta']))
        queries.append(row);family[row['family_id']].append(row)
    families=[]
    for fam,rows in sorted(family.items()):
        idx=[i for i,r in enumerate(records) if r['family_id']==fam]
        reps={records[i]['receptor_id']:records[i] for i in idx}
        families.append(dict(family_id=fam,fold=folds[idx[0]],n_graphs=len(idx),n_receptors=len(reps),
                             n_eligible_queries=len(rows),domain_counts=dict(Counter(records[i]['domain'] for i in idx)),
                             peptide_counts=dict(Counter(records[i]['peptide'] for i in idx)),
                             composition_auc=float(np.mean([r['composition_auc'] for r in rows])),
                             model_auc=float(np.mean([r['model_auc'] for r in rows])),
                             delta=float(np.mean([r['model_minus_composition'] for r in rows]))))
    norm=[];hashgroups=defaultdict(list)
    for i,g in enumerate(graphs):
        h=torch.tensor(g[0],dtype=torch.float32);z=normalize_node_features(h)
        norm.append(dict(index=i,raw_identity_mean_l2=float(h[:,4:24].mean(0).norm()),
                         normalized_identity_mean_max_abs=float(z[:,4:24].mean(0).abs().max())))
        b=b''.join(np.ascontiguousarray(x).tobytes() for x in g)
        hashgroups[hashlib.sha256(b).hexdigest()].append(i)
    collisions=[idx for idx in hashgroups.values() if len({records[i]['receptor_id'] for i in idx})>1]
    batch_coverage=[]
    for f in range(5):
        train=[i for i,v in enumerate(folds) if v!=f]
        for seed in [20261021,1101,4242]:
            rng=random.Random(seed);total=Counter()
            for epoch in range(120):
                order=train[:];rng.shuffle(order);total.update(coverage(records,order))
            batch_coverage.append(dict(fold=f,seed=seed,**total))
    result=dict(scope='retrospective diagnostics only; no fitting or new independent data',source_sha256=source,
                families=families,queries=queries,node_normalization=norm,
                exact_graph_collisions_across_receptors=collisions,random_training_pair_coverage=batch_coverage,
                developmental_partition=dict(fit_original_folds=[2,3,4],validation_original_folds=[1],quarantine_original_folds=[0],
                    counts={str(f):sum(v==f for v in folds) for f in range(5)},
                    eligible_families_by_fold=dict(Counter(x['fold'] for x in families))))
    a.out_dir.mkdir(parents=True)
    (a.out_dir/'diagnostics.json').write_text(json.dumps(finite_json(result),indent=2,allow_nan=False)+'\n')
    lines=['# Retrospective hard-family error analysis','',
           'All 331 graphs and model outputs were exposed previously. Findings motivate developmental ablations; no independent evaluation is claimed.',
           '', '| Family | Fold | Eligible queries | Composition | Model | Delta |', '|---|---:|---:|---:|---:|---:|']
    for x in sorted(families,key=lambda x:x['delta']):
        lines.append(f"| {x['family_id']} | {x['fold']} | {x['n_eligible_queries']} | {x['composition_auc']:.4f} | {x['model_auc']:.4f} | {x['delta']:.4f} |")
    total=Counter()
    for x in batch_coverage:
        total.update({k:v for k,v in x.items() if k not in ('fold','seed')})
    lines+=['',f"Active positive anchors / all: {total['anchors_with_positive']}/{total['anchors']}.",
            f"Active anchors also seeing same-family negatives / active: {total['active_anchors_with_hard_negative']}/{total['anchors_with_positive']}.",
            f"Maximum identity-column mean after z-score: {max(x['normalized_identity_mean_max_abs'] for x in norm):.9g}.",
            f"Exact identical graph tuples across different receptor labels: {len(collisions)} groups.",
            'Zeroing one-hot means removes the direct composition mean, but does not prove all identity information is erased.',
            'Edit distances refer to complete recorded chains; interface observability requires node-to-residue mapping and is not inferred here.']
    (a.out_dir/'README.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'families':families,'coverage':dict(total),'partition':result['developmental_partition']},indent=2))

if __name__=='__main__':main()
