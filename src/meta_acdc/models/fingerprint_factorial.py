"""Controlled diagnostics of observed-sequence labels; no biological identity promotion."""
from collections import defaultdict,Counter
import hashlib
import math
import random
from pathlib import Path
import numpy as np
import torch
from meta_acdc.models.train_fingerprint import batch_graphs
from meta_acdc.models.train_fingerprint_v2 import EmbedderV2,set_seed,supcon_loss_v2
from meta_acdc.structure.feature_schema import AA_START,AA_STOP
from scripts.diagnose_fingerprint_hard_family import coverage
from scripts.rashomon_fingerprint import bootstrap_ci
from sklearn.metrics import roc_auc_score

ARMS={'reference':dict(raw_identity=False,structured=False),
      'identity_preserved':dict(raw_identity=True,structured=False),
      'family_batches':dict(raw_identity=False,structured=True),
      'identity_and_family':dict(raw_identity=True,structured=True)}

def assemble(graphs,records,indices,raw_identity=False,intervention=None):
    h,x,ei,ea,bt,pm,rids=batch_graphs(graphs,records,indices,'cpu',None)
    if raw_identity:
        h[:,AA_START:AA_STOP]=torch.cat([torch.as_tensor(graphs[i][0],dtype=torch.float32)[:,AA_START:AA_STOP] for i in indices])
    if intervention=='identity_zero':h[:,AA_START:AA_STOP]=0
    if intervention=='numeric_geometry_zero':
        x=torch.zeros_like(x);ea=torch.zeros_like(ea);h[:,25]=0
    return h,x,ei,ea,bt,pm,rids

def epoch_order(records,indices,rng,structured):
    if not structured:
        order=list(indices);rng.shuffle(order);return order
    families=defaultdict(lambda:defaultdict(list))
    for i in indices:families[records[i]['family_id']][records[i]['receptor_id']].append(i)
    keys=sorted(families);rng.shuffle(keys);order=[]
    for family in keys:
        rids=sorted(families[family]);rng.shuffle(rids)
        for rid in rids:
            views=list(families[family][rid]);rng.shuffle(views);order.extend(views)
    return order

def tensor_digest(state):
    digest=hashlib.sha256()
    for key,value in sorted(state.items()):
        digest.update(key.encode());digest.update(str(value.shape).encode());digest.update(value.detach().cpu().numpy().tobytes())
    return digest.hexdigest()

def encode(model,graphs,records,indices,raw_identity,intervention=None):
    model.eval();values=[]
    with torch.no_grad():
        for i in indices:
            batch=assemble(graphs,records,[i],raw_identity,intervention)
            if intervention=='peptide_pool_zero':
                pooled=model.egnn(*batch[:5],embed_only=True,peptide_mask=batch[5])
                hidden=model.egnn.embed.out_features;pooled[:,2*hidden:3*hidden]=0
                out=model.proj(pooled)
            else:out=model(*batch[:6])
            values.append(out.squeeze(0).cpu())
    E=torch.stack(values);return E/E.norm(dim=1,keepdim=True).clamp(min=1e-8)

def train_model(graphs,records,fit_idx,seed,arm,recipe):
    set_seed(seed);torch.set_num_threads(recipe['threads'])
    model=EmbedderV2(26,16,emb_dim=recipe['emb_dim']);initial=tensor_digest(model.state_dict())
    opt=torch.optim.Adam(model.parameters(),lr=recipe['lr'],weight_decay=recipe['weight_decay'])
    rng=random.Random(seed);trace=[];total=Counter();options=ARMS[arm]
    for epoch in range(recipe['epochs']):
        model.train();order=epoch_order(records,fit_idx,rng,options['structured']);counts=Counter();losses=[]
        if epoch<recipe['warmup_epochs']:lr=recipe['lr']*(epoch+1)/recipe['warmup_epochs']
        else:
            t=(epoch-recipe['warmup_epochs'])/max(1,recipe['epochs']-recipe['warmup_epochs'])
            lr=recipe['lr_min']+0.5*(recipe['lr']-recipe['lr_min'])*(1+math.cos(math.pi*t))
        for group in opt.param_groups:group['lr']=lr
        for start in range(0,len(order),recipe['batch_size']):
            idx=order[start:start+recipe['batch_size']]
            batch=assemble(graphs,records,idx,options['raw_identity'])
            emb=model(*batch[:6]);emb=emb/emb.norm(dim=1,keepdim=True).clamp(min=1e-8)
            loss=supcon_loss_v2(emb,batch[6],[records[i]['family_id'] for i in idx],recipe['temperature'],1.0)
            if not torch.isfinite(loss):raise RuntimeError('nonfinite loss')
            opt.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),10.0,error_if_nonfinite=True);opt.step()
            losses.append(float(loss.detach()));counts.update(coverage(records,idx,recipe['batch_size']))
        trace.append(dict(epoch=epoch,lr=lr,mean_batch_loss=float(np.mean(losses)),coverage=dict(counts)))
        total.update(counts)
        if epoch%15==0 or epoch==recipe['epochs']-1:print(f'{arm} seed={seed} epoch={epoch} loss={np.mean(losses):.4f}',flush=True)
    return model,initial,trace,dict(total)

def auc_rows(S,records,indices,known_same=(),hard=True):
    excluded=set(map(tuple,known_same));rows=[];removed=0
    for q in indices:
        gallery=[j for j in indices if j!=q and (not hard or records[j]['receptor_id']==records[q]['receptor_id'] or records[j]['family_id']==records[q]['family_id'])]
        clean=[]
        for j in gallery:
            if records[j]['receptor_id']!=records[q]['receptor_id'] and (q,j) in excluded:removed+=1
            else:clean.append(j)
        labels=np.array([records[j]['receptor_id']==records[q]['receptor_id'] for j in clean],dtype=int)
        rows.append(float(roc_auc_score(labels,S[q,clean])) if labels.any() and (labels==0).any() else np.nan)
    array=np.asarray(rows)
    summary=bootstrap_ci(array,[records[q]['family_id'] for q in indices])
    summary['n_excluded_known_same_negative_pairs']=removed
    summary['query_macro_mean']=float(np.mean(array[np.isfinite(array)])) if np.isfinite(array).any() else None
    return summary,array
