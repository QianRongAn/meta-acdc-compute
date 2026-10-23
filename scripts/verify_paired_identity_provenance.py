"""Post-freeze provenance replay and detailed quarantine evidence; does not alter admission."""
from pathlib import Path
import hashlib,json,subprocess,os,tempfile
from collections import Counter,defaultdict
import torch
from meta_acdc.structure.receptor_identity_audit import seqres_sequences
from meta_acdc.structure.paired_variable_identity import cif_declared_sequences
ROOT=Path(__file__).resolve().parents[1];MAIN=Path('/home/administrator/docs/meta-acdc')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    base=ROOT/'docs/fingerprint_paired_identity_20261024';original=base/'curation_v3'
    replay=base/'curation_replay'
    if replay.exists():raise FileExistsError(replay)
    subprocess.run([str(MAIN/'.venv/bin/python'),'scripts/curate_fingerprint_paired_identity.py','--out-dir',str(replay)],cwd=ROOT,check=True)
    assert sha(original/'identity_registry.json')==sha(replay/'identity_registry.json'),'curation replay differs'
    registry=json.loads((original/'identity_registry.json').read_text());old_audit=json.loads((ROOT/'docs/fingerprint_hard_family_20261024/identity_audit/identity_audit.json').read_text())
    old=torch.load(ROOT/'docs/fingerprint_retrain_20261023/input/fingerprint_dataset.pt',map_location='cpu',weights_only=False)
    inv=torch.load(MAIN/'data/releases/20261017-structural-v7-atomic-repair/native_inventory_graphs.pt',map_location='cpu',weights_only=False);assets={x['record']['complex_id']:x['record'] for x in inv['assets']}
    source_checks=[];details=[];crystal={x['index']:x for x in old_audit['rows']}
    for row in registry['rows']:
        i=row['index'];rec=old['records'][i];detail=dict(index=i,complex_id=rec['complex_id'],status=row['status'],reasons=row['reasons'],observed_alpha=rec['alpha'],observed_beta=rec['beta'])
        if i in crystal:
            meta=assets[rec['complex_id']];path=Path(meta['source_path']);assert sha(path)==meta['source_sha256']==crystal[i]['source_sha256']
            seq=seqres_sequences(path.read_text());detail.update(source_sha256=sha(path),alpha_chain=meta['alpha_chain'],beta_chain=meta['beta_chain'],full_alpha_construct=seq.get(meta['alpha_chain']),full_beta_construct=seq.get(meta['beta_chain']))
            for role in ['alpha','beta']:
                if row.get(role,{}).get('boundary_status')=='pass':
                    v=row[role];assert v['full_construct'][v['variable_start0']:v['variable_end0']]==v['variable']
        else:
            if rec['source']=='v7-inventory':path=Path(assets[rec['complex_id']]['source_path']);expected=assets[rec['complex_id']]['source_sha256']
            else:path=MAIN/'data/raw/af3_predictions'/(rec['complex_id']+'_model_0.cif');expected=old['graph_hashes'][i]
            assert sha(path)==expected;declared=cif_declared_sequences(path.read_text());detail.update(source_sha256=expected,declared_cif_chains=declared,matching_alpha_chains=[c for c,s in declared.items() if s==rec['alpha']],matching_beta_chains=[c for c,s in declared.items() if s==rec['beta']])
        details.append(detail);source_checks.append(i)
    admitted=[x for x in registry['rows'] if x['status']=='admitted'];groups=defaultdict(list)
    for row in admitted:groups[row['receptor_id']].append(row)
    merges=[dict(corrected_receptor_id=rid,old_observed_ids=sorted({x['old_receptor_id'] for x in rows}),indices=[x['index'] for x in rows]) for rid,rows in groups.items() if len({x['old_receptor_id'] for x in rows})>1]
    (base/'QUARANTINE_AND_CONSTRUCT_DETAILS.json').write_text(json.dumps(dict(scope='source evidence only; no admission or training changes',rows=details),indent=2)+'\n')
    result=dict(status='pass',curation_registry_byte_replay_exact=True,registry_sha256=sha(original/'identity_registry.json'),n_source_hashes_checked=len(source_checks),n_admitted_observed_keys=len({x['old_receptor_id'] for x in admitted}),n_corrected_variable_keys=len(groups),n_confirmed_merge_groups=len(merges),confirmed_merges=merges,independent_acceptance_accessed=False)
    (base/'provenance_verification.json').write_text(json.dumps(result,indent=2)+'\n');print('PROVENANCE_PASS',len(source_checks),'MERGES',len(merges))
if __name__=='__main__':main()
