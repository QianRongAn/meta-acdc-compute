"""Freeze new development-only leave-one-eligible-V-gene-context-out diagnostics."""
from pathlib import Path
import hashlib,json
from collections import defaultdict,Counter
import torch
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    base=ROOT/'docs/fingerprint_paired_identity_20261024';out=base/'development'
    if out.exists():raise FileExistsError(out)
    out.mkdir();registry=json.loads((base/'curation_v3/identity_registry.json').read_text())
    old=torch.load(ROOT/'docs/fingerprint_retrain_20261023/input/fingerprint_dataset.pt',map_location='cpu',weights_only=False)
    accepted=[x for x in registry['rows'] if x['status']=='admitted' and x['development_partition'] in ['fit','validation']]
    records=[];graphs=[]
    for row in accepted:
        i=row['index'];rec=dict(old['records'][i]);rec.update(receptor_id=row['receptor_id'],construct_id=row['construct_id'],source_index=i,old_family_id=rec['family_id'],curated_family_id=row['family_id'],alpha_variable=row['alpha']['variable'],beta_variable=row['beta']['variable'],alpha_v_gene=row['alpha']['v_allele'].split('*')[0],beta_v_gene=row['beta']['v_allele'].split('*')[0],species=row['alpha'].get('species','unknown'))
        records.append(rec);graphs.append(old['graphs'][i])
    parents=list(range(len(records)))
    def find(i):
        while parents[i]!=i:parents[i]=parents[parents[i]];i=parents[i]
        return i
    seen={}
    for i,r in enumerate(records):
        keys=[('old',r['curated_family_id']),('alpha',r['alpha_variable']),('beta',r['beta_variable']),('V_pair',r['species'],r['alpha_v_gene'],r['beta_v_gene'])]
        for k in keys:
            if k in seen:parents[find(i)]=find(seen[k])
            else:seen[k]=i
    groups=defaultdict(list)
    for i in range(len(records)):groups[find(i)].append(i)
    contexts=[]
    for indices in groups.values():
        key='G'+hashlib.sha256('|'.join(sorted({records[i]['receptor_id'] for i in indices})).encode()).hexdigest()[:16]
        for i in indices:records[i]['family_id']=key
        counts=Counter(records[i]['receptor_id'] for i in indices)
        contexts.append(dict(family_id=key,indices=indices,n_paired_variable_ids=len(counts),n_eligible_hard_queries=sum(n for n in counts.values() if n>1) if len(counts)>1 else 0))
    hard=sorted([g for g in contexts if g['n_eligible_hard_queries']>0],key=lambda g:g['family_id'])
    folds=[]
    for g in hard:
        val=g['indices'];fit=sorted(set(range(len(records)))-set(val))
        assert not {records[i]['receptor_id'] for i in fit}&{records[i]['receptor_id'] for i in val}
        assert not {records[i]['alpha_variable'] for i in fit}&{records[i]['alpha_variable'] for i in val}
        assert not {records[i]['beta_variable'] for i in fit}&{records[i]['beta_variable'] for i in val}
        folds.append(dict(fold_id=g['family_id'],fit_indices=fit,validation_indices=val))
    torch.save(dict(graphs=graphs,records=records,source_indices=[r['source_index'] for r in records],complete_variable_identity_admitted=True,scope='deposited engineered paired-variable sequence retrieval, not binding or independent generalization'),out/'dataset.pt')
    pold=json.loads((ROOT/'docs/fingerprint_hard_family_20261024/protocol.json').read_text())
    sources=['src/meta_acdc/structure/paired_variable_identity.py','scripts/curate_fingerprint_paired_identity.py','scripts/freeze_paired_identity_development.py','scripts/run_paired_identity_development.py','src/meta_acdc/models/fingerprint_factorial.py','src/meta_acdc/models/train_fingerprint.py','src/meta_acdc/models/train_fingerprint_v2.py','src/meta_acdc/models/egnn.py','scripts/rashomon_fingerprint.py','scripts/diagnose_fingerprint_hard_family.py','src/meta_acdc/structure/feature_schema.py']
    sources+=['docs/fingerprint_paired_identity_20261024/curation_v3/identity_registry.json','docs/fingerprint_paired_identity_20261024/development/dataset.pt']
    p=dict(protocol_id='PAIRED-VARIABLE-DEV-20261024-v1',scope='historical exposed structural panel only; no binding labels, no independent acceptance; new metadata-defined development cross-validation',dataset=str((out/'dataset.pt').relative_to(ROOT)),registry='docs/fingerprint_paired_identity_20261024/curation_v3/identity_registry.json',n_graphs=len(records),n_paired_variable_ids=len({r['receptor_id'] for r in records}),contexts=contexts,folds=folds,seeds=[20261021,1101,4242],arms=['reference','identity_preserved','identity_preserved_shuffled'],recipe=pold['recipe'],epochs_selected='fixed final 120; no seed, fold, epoch or arm selection',batch_factor='not tested: every fit partition fits in one batch, so structured ordering cannot change pair exposure',partition_change='old fixed fit/validation has zero defined same-family hard queries after identity repair; newly freeze leave-one-eligible-V-gene-context-out DEVELOPMENT folds from old folds1-4 only; not a new test set; no performance-driven reassignment',family_definition='union of original curated-family links, identical full single variable chains and exact species+alpha/beta V-gene pair; no same complete chain/V-pair crossing fit/validation',hard_metric='per-query AUROC against same context, equal query mean within context, then equal context mean; positive=identical complete paired variable amino acids; negative=different paired variable sequence, NOT nonbinding',baselines=['constant 0.5','raw residue composition cosine','graph-size negative Euclidean after training-only standardization'],label_shuffle='permute fit graph identity labels only, preserve counts, deterministic local RNG(seed+9173); validation truth unchanged; one fixed identity-preserved control',promotion=False,limitations=['few contexts; descriptive conditional intervals only','positives may be multiple copies of one crystal, not independent experiments','engineering and constant-region differences tracked separately; no native-clone or clinical identity inferred'],independent_acceptance_accessed=False,frozen_sha256={s:sha(ROOT/s) for s in sources})
    (out/'protocol.json').write_text(json.dumps(p,indent=2)+'\n');(out/'protocol.sha256').write_text(sha(out/'protocol.json')+'  protocol.json\n')
    print('FROZEN',len(records),'graphs',len(hard),'hard contexts',sum(x['n_eligible_hard_queries'] for x in hard),'hard queries',[(len(x['fit_indices']),len(x['validation_indices'])) for x in folds],sha(out/'protocol.json'))
if __name__=='__main__':main()
