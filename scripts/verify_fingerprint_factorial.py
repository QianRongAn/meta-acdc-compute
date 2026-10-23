"""Artifact and boundary acceptance for the frozen diagnostic factorial experiment."""
import argparse
from collections import Counter,defaultdict
import json
from pathlib import Path
import random
import sys
import numpy as np
import torch
from meta_acdc.models.train_fingerprint_v2 import EmbedderV2
from meta_acdc.models.fingerprint_factorial import ARMS,encode,auc_rows,epoch_order
from scripts.diagnose_fingerprint_hard_family import coverage
from scripts.run_fingerprint_factorial import validate,ROOT,sha
from scripts.rashomon_fingerprint import finite_json,bootstrap_ci
from scripts.fingerprint_shortcut_audit import cosine_rows

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--protocol',type=Path,required=True);ap.add_argument('--run',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    p=validate(a.protocol);torch.set_num_threads(p['recipe']['threads']);torch.use_deterministic_algorithms(True)
    data=torch.load(ROOT/p['dataset'],map_location='cpu',weights_only=False);records=data['records']
    if json.loads((a.run/'status.json').read_text())['state']!='complete':raise ValueError('incomplete run')
    if set(data['source_indices'])&set(p['quarantined_source_indices']):raise ValueError('quarantined source in development')
    checks=[];rows={};traces={}
    instances=[(arm,s,a.run/arm/str(s)) for s in p['seeds'] for arm in ARMS]+[('reference',p['seeds'][0],a.run/'replay')]
    for arm,seed,out in instances:
        metric=json.loads((out/'metrics.json').read_text())
        if metric['protocol_sha256']!=p['protocol_sha256'] or metric['arm']!=arm or metric['seed']!=seed:raise ValueError('instance binding failed')
        for name,digest in metric['artifact_sha256'].items():
            if sha(out/name)!=digest:raise ValueError('artifact changed')
        checkpoint=torch.load(out/'model.pt',map_location='cpu',weights_only=False)
        if not checkpoint['diagnostic_only'] or checkpoint['identity_truth_validated']:raise ValueError('unjustified biological promotion')
        model=EmbedderV2(26,16,emb_dim=p['recipe']['emb_dim']);model.load_state_dict(checkpoint['model'],strict=True)
        E=encode(model,data['graphs'],records,list(range(len(records))),ARMS[arm]['raw_identity'])
        S=(E@E.T).numpy()
        if not np.array_equal(E.numpy(),np.load(out/'embeddings.npy')) or not np.array_equal(S,np.load(out/'similarities.npy')):
            raise ValueError('weight to score reproduction failed')
        for partition,idx in [('fit',p['fit_indices']),('validation',p['validation_indices'])]:
            for name,mask,hard in [('raw_hard',[],True),('known_same_excluded',p['known_same_pairs'],True),('all',[],False)]:
                summary,values=auc_rows(S,records,idx,mask,hard)
                if finite_json(dict(summary=summary,query_auc=values.tolist()))!=metric['partitions'][partition][name]:
                    raise ValueError('metric mismatch')
        for intervention in p['interventions']:
            idx=p['validation_indices'];e=encode(model,data['graphs'],records,idx,ARMS[arm]['raw_identity'],intervention)
            if not np.array_equal(e.numpy(),np.load(out/(intervention+'.npy'))):raise ValueError('intervention reproduction failed')
            summary,values=auc_rows((e@e.T).numpy(),[records[i] for i in idx],list(range(len(idx))))
            if finite_json(dict(summary=summary,query_auc=values.tolist()))!=metric['interventions'][intervention]:
                raise ValueError('intervention metric mismatch')
        rng=random.Random(seed);total=Counter();trace=json.loads((out/'trace.json').read_text())
        for epoch in range(p['recipe']['epochs']):
            order=epoch_order(records,p['fit_indices'],rng,ARMS[arm]['structured'])
            if sorted(order)!=sorted(p['fit_indices']):raise ValueError('sampling loss/duplication or validation contamination')
            counts=coverage(records,order,p['recipe']['batch_size']);total.update(counts)
            if counts!=trace[epoch]['coverage']:raise ValueError('coverage mismatch')
        if dict(total)!=metric['coverage']:raise ValueError('coverage totals mismatch')
        checks.append(dict(arm=arm,seed=seed,location=out.relative_to(a.run).as_posix(),
                           weight_to_embedding_exact=True,similarity_exact=True,metrics_exact=True,
                           interventions_exact=True,fit_sampling_exact=True,diagnostic_gate_preserved=True))
        if out.name!='replay':rows[(arm,seed)]=metric
        print('VERIFIED',len(checks),out.relative_to(a.run),flush=True)
    comparisons={};family_errors=[]
    comp=cosine_rows(np.array([np.asarray(g[0])[:,4:24].mean(0) for g in data['graphs']]))
    for part,idx in [('fit',p['fit_indices']),('validation',p['validation_indices'])]:
        groups=[records[i]['family_id'] for i in idx];comparisons[part]={}
        br=auc_rows(comp,records,idx)[1];bc=auc_rows(comp,records,idx,p['known_same_pairs'])[1]
        comparisons[part]['composition_mask_minus_raw_common']=bootstrap_ci(bc-br,groups)
        for arm in ARMS:
            raw=np.array([rows[(arm,s)]['partitions'][part]['raw_hard']['query_auc'] for s in p['seeds']],dtype=float).T
            clean=np.array([rows[(arm,s)]['partitions'][part]['known_same_excluded']['query_auc'] for s in p['seeds']],dtype=float).T
            comparisons[part][arm]=bootstrap_ci(clean-raw,groups)
            for family in sorted(set(groups)):
                positions=[j for j,g in enumerate(groups) if g==family]
                vals=raw[positions];valid=np.isfinite(vals)
                if not valid.any():continue
                family_errors.append(dict(partition=part,arm=arm,family_id=family,
                                          n_queries=int(np.isfinite(vals).any(1).sum()),
                                          mean=float(vals[valid].mean()),
                                          per_seed=[float(np.mean(v[np.isfinite(v)])) if np.isfinite(v).any() else None for v in vals.T]))
    validate(a.protocol)
    a.out.write_text(json.dumps(finite_json(dict(status='pass',n_models=len(checks),checks=checks,
                        common_query_pair_mask_sensitivity=comparisons,per_family_errors=family_errors,
                        source_quarantine_count=len(p['quarantined_source_indices']),
                        protocol_sha256=p['protocol_sha256'],measured_acceptance_sets_read=False,
                        biological_candidate_promotable=False,verifier_scope='post-training computational acceptance only',
                        environment=dict(python=sys.version,torch=torch.__version__,numpy=np.__version__))),indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
