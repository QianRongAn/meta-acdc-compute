"""Recheck every old isolated item against new exact full-variable evidence."""
from pathlib import Path
import json,hashlib
from collections import defaultdict,Counter
import torch
from meta_acdc.structure.receptor_identity_audit import seqres_sequences,is_subsequence
from meta_acdc.structure.paired_variable_identity import fasta_records,cif_declared_sequences
ROOT=Path(__file__).resolve().parents[1];MAIN=Path('/home/administrator/docs/meta-acdc')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    base=ROOT/'docs/fingerprint_cross_source_20261024';oldbase=ROOT/'docs/fingerprint_paired_identity_20261024';old=json.loads((oldbase/'curation_v3/identity_registry.json').read_text());new=json.loads((base/'source_balanced_design.json').read_text());native=json.loads((ROOT/'docs/fingerprint_hard_family_20261024/identity_audit/identity_audit.json').read_text());old_native={x['index']:x for x in native['rows']};records=json.loads((oldbase/'record_index.json').read_text());inventory=torch.load(MAIN/'data/releases/20261017-structural-v7-atomic-repair/native_inventory_graphs.pt',map_location='cpu',weights_only=False);assets={x['record']['complex_id']:x['record'] for x in inventory['assets']};details={x['index']:x for x in json.loads((oldbase/'QUARANTINE_AND_CONSTRUCT_DETAILS.json').read_text())['rows']}
    admitted=[x for x in new['rows'] if x['status']=='identity_admitted_structure_QC_pending'];bychain={(x['pdb'],x['alpha_chain'],x['beta_chain']):x for x in admitted};pairs={x['receptor_id']:x for x in admitted};role_refs=defaultdict(list)
    for p in (ROOT/'af3_pending_upload').rglob('*.fasta'):
        rec=dict(fasta_records(p.read_text()))
        if set(rec)=={'alpha','beta','peptide','mhc','b2m'}:role_refs[(rec['alpha'],rec['beta'])].append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
    rows=[]
    for item in old['rows']:
        if item['status']!='quarantined':continue
        i=item['index'];record=records[i];row=dict(index=i,complex_id=item['complex_id'],domain=item['domain'],previous_reasons=item['reasons'],status='still_isolated',reasons=[],source_hash_rechecked=False)
        try:
            if i in old_native:
                original=old_native[i];meta=assets[record['complex_id']];path=Path(meta['source_path']);assert sha(path)==meta['source_sha256']==original['source_sha256'];row['source_hash_rechecked']=True
                proof=bychain.get((original['pdb'],original['alpha_chain'],original['beta_chain']))
                if proof is None:raise ValueError('no new complete paired-domain evidence for this exact primary chain pair')
                if not original['eligible_construct_evidence']:raise ValueError('old role/construct mapping was invalid; source graph must be rebuilt and verified, not relabeled')
                seq=seqres_sequences(path.read_text())
                for role in ['alpha','beta']:
                    if seq.get(original[role+'_chain'])!=proof[role]['full_construct']:raise ValueError('old and new deposited construct versions differ')
                    if not is_subsequence(record[role],proof[role]['full_construct']):raise ValueError('old observed chain does not agree with full construct')
                row['proof_primary_pdb']=proof['pdb']
            else:
                if record['source']=='v7-inventory':path=Path(assets[record['complex_id']]['source_path'])
                else:path=MAIN/'data/raw/af3_predictions'/(record['complex_id']+'_model_0.cif')
                assert sha(path)==details[i]['source_sha256'];row['source_hash_rechecked']=True;declared=cif_declared_sequences(path.read_text());pair=(record['alpha'],record['beta'])
                if pair not in role_refs:raise ValueError('no exact role FASTA; no role inference from nearest sequence')
                if not any(sa==pair[0] and sb==pair[1] and a!=b for a,sa in declared.items() for b,sb in declared.items()):raise ValueError('declared CIF complete chains do not match record')
                matches=[x for x in pairs.values() if pair[0].count(x['alpha']['variable'])==1 and pair[1].count(x['beta']['variable'])==1]
                if len(matches)!=1:raise ValueError('actual AF3 input still lacks a unique complete paired-domain match; no residue filling')
                proof=matches[0];row['role_fasta_evidence']=role_refs[pair]
            row.update(status='identity_evidence_resolved_release_still_blocked',receptor_id=proof['receptor_id'],alpha_variable=proof['alpha']['variable'],beta_variable=proof['beta']['variable'],V_context=proof['V_context'],identity_proof=dict(pdb=proof['pdb'],alpha_chain=proof['alpha_chain'],beta_chain=proof['beta_chain'],primary_header_sha256=proof['primary_header_sha256']),source_sha256=sha(path),old_observed_id=item['old_receptor_id'],not_a_training_release=True)
        except (ValueError,KeyError,OSError) as error:row['reasons'].append(str(error))
        rows.append(row)
    assert len(rows)==285 and all(x['source_hash_rechecked'] for x in rows)
    summary=dict(n_previous_isolated=285,n_identity_evidence_resolved=sum(x['status']=='identity_evidence_resolved_release_still_blocked' for x in rows),n_still_isolated=sum(x['status']=='still_isolated' for x in rows),n_original_source_hashes_rechecked=sum(x['source_hash_rechecked'] for x in rows),counts_by_domain_status=dict(Counter(x['domain']+'/'+x['status'] for x in rows)))
    (base/'isolation_followup.json').write_text(json.dumps(dict(summary=summary,rows=rows,old_release_unchanged=True,no_automatic_training_admission=True,independent_acceptance_accessed=False),indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
