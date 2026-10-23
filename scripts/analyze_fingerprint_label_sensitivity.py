"""Predeclared sensitivity to confirmed false negative pairs; no model fitting."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import torch
from meta_acdc.models.fingerprint_factorial import auc_rows
from scripts.fingerprint_shortcut_audit import cosine_rows
from scripts.rashomon_fingerprint import bootstrap_ci,finite_json
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    root=ROOT/'docs/fingerprint_hard_family_20261024'
    old=ROOT/'docs/fingerprint_retrain_20261023'
    d=torch.load(old/'input/fingerprint_dataset.pt',map_location='cpu',weights_only=False)
    records=d['records'];n=len(records);idx=list(range(n));groups=[r['family_id'] for r in records]
    audit=json.loads((root/'identity_audit/identity_audit.json').read_text());known=audit['known_same_crystal_pairs']
    fold=json.loads((old/'fold_assignment.json').read_text())['per_graph_fold']
    Scomp=cosine_rows(np.array([np.asarray(g[0])[:,4:24].mean(0) for g in d['graphs']]))
    raw_base,raw_rows=auc_rows(Scomp,records,idx)
    clean_base,clean_rows=auc_rows(Scomp,records,idx,known)
    raw_models=[];clean_models=[];source_sha={}
    for seed in [20261021,1101,4242]:
        S=np.zeros((n,n))
        for f in range(5):
            path=old/'run/true/model_run_v2'/f'sim_fold{f}_seed{seed}.npy'
            source_sha[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
            scores=np.load(path,allow_pickle=False);rows=[q for q,v in enumerate(fold) if v==f];S[rows]=scores[rows]
        raw_models.append(auc_rows(S,records,idx)[1]);clean_models.append(auc_rows(S,records,idx,known)[1])
    raw=np.array(raw_models).T;clean=np.array(clean_models).T
    result=dict(scope='retrospective pair-mask diagnostic; no biological identity repair, no new test',
                known_same_negative_pairs=len(known),composition_raw=raw_base,composition_masked=clean_base,
                model_raw=bootstrap_ci(raw,groups),model_masked=bootstrap_ci(clean,groups),
                model_mask_minus_raw_on_common_queries=bootstrap_ci(clean-raw,groups),
                composition_mask_minus_raw_on_common_queries=bootstrap_ci(clean_rows-raw_rows,groups),
                raw_model_minus_composition=bootstrap_ci(raw-raw_rows[:,None],groups),
                masked_model_minus_composition=bootstrap_ci(clean-clean_rows[:,None],groups),
                source_sha256=source_sha,
                caveat='Only primary-backed identical full-alpha/beta crystal pairs are masked. AF3 full-construct identities and cross-PDB variable-region equivalence are not imputed. Undefined queries/families are counted; no selected-family deletion.')
    a.out.write_text(json.dumps(finite_json(result),indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['model_raw','model_masked','masked_model_minus_composition']},indent=2))

if __name__=='__main__':main()
