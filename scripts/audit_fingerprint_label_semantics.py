"""Audit existing labels against local depositor evidence; never relabel or merge by similarity."""
import argparse
import hashlib
import json
from pathlib import Path
from collections import Counter
import torch
from meta_acdc.structure.identity import pdb_entities
from meta_acdc.structure.receptor_identity_audit import seqres_sequences,is_subsequence,same_construct_conflicts,IDENTITY_SEMANTICS
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inventory',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True);args=ap.parse_args()
    if args.out_dir.exists():raise FileExistsError(args.out_dir)
    data_path=ROOT/'docs/fingerprint_retrain_20261023/input/fingerprint_dataset.pt'
    d=torch.load(data_path,map_location='cpu',weights_only=False);records=d['records'];by_id={r['complex_id']:i for i,r in enumerate(records)}
    inv=torch.load(args.inventory,map_location='cpu',weights_only=False)
    args.out_dir.mkdir(parents=True);headerdir=args.out_dir/'headers';headerdir.mkdir()
    rows=[];sources={}
    for asset in inv['assets']:
        meta=asset['record'];cid=meta['complex_id']
        if cid not in by_id or meta['domain']!='crystal':continue
        path=Path(meta['source_path']);raw=path.read_bytes();sha=hashlib.sha256(raw).hexdigest()
        if sha!=meta['source_sha256']:raise ValueError('source PDB changed')
        text=raw.decode('utf-8');seq=seqres_sequences(text);entities=pdb_entities(path)
        alpha,beta=meta['alpha_chain'],meta['beta_chain'];r=records[by_id[cid]]
        valid=(alpha in seq and beta in seq and entities.get(alpha,{}).get('kind')=='tra'
               and entities.get(beta,{}).get('kind')=='trb'
               and is_subsequence(r['alpha'],seq[alpha]) and is_subsequence(r['beta'],seq[beta]))
        row=dict(index=by_id[cid],complex_id=cid,pdb=meta['pdb'],family_id=r['family_id'],
                 observed_sequence_id=r['receptor_id'],alpha_chain=alpha,beta_chain=beta,
                 alpha_entity=entities.get(alpha),beta_entity=entities.get(beta),
                 observed_alpha=r['alpha'],observed_beta=r['beta'],seqres_alpha=seq.get(alpha),seqres_beta=seq.get(beta),
                 eligible_construct_evidence=bool(valid),source_sha256=sha)
        rows.append(row)
        if meta['pdb'] not in sources:
            header='\n'.join(line for line in text.splitlines() if line.startswith(('HEADER','TITLE','COMPND','SEQRES','SEQADV','DBREF','SOURCE')))+'\n'
            target=headerdir/(meta['pdb']+'.txt');target.write_text(header)
            sources[meta['pdb']]=dict(original_path=str(path),original_sha256=sha,original_bytes=len(raw),
                                     extracted_header_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    conflicts=same_construct_conflicts(rows)
    pairs=[]
    for group in conflicts:
        for a in group:
            for b in group:
                if a['observed_sequence_id']!=b['observed_sequence_id']:
                    pairs.append([a['index'],b['index']])
    out=dict(identity_semantics=IDENTITY_SEMANTICS,paired_variable_identity_validated=False,
             scope='exact deposited construct equivalence proves some original negatives invalid; not complete variable-region curation',
             source_inventory_sha256=hashlib.sha256(args.inventory.read_bytes()).hexdigest(),
             input_sha256=hashlib.sha256(data_path.read_bytes()).hexdigest(),rows=rows,sources=sources,
             confirmed_same_construct_conflicts=conflicts,known_same_crystal_pairs=pairs,
             n_conflict_groups=len(conflicts),n_conflicting_observed_ids=len({r['observed_sequence_id'] for g in conflicts for r in g}),
             no_automatic_merges=True,no_AF3_full_construct_identity_inferred=True,
             missing='audited paired complete variable sequences, terminal/tag/constant-region boundaries and AF3 construct correspondence')
    (args.out_dir/'identity_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['n_conflict_groups','n_conflicting_observed_ids','paired_variable_identity_validated']},indent=2))

if __name__=='__main__':main()
