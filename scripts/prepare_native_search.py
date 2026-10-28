"""Prepare metadata-only paired-chain discovery; never complete truncated sequences."""
from pathlib import Path
import argparse
import hashlib
import json
from urllib.parse import quote

from meta_acdc.structure.paired_variable_identity import numbered_domains

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    previous = ROOT/'docs/fingerprint_positive_controls_20261027/source_balanced_design.json'
    followup = ROOT/'docs/fingerprint_construct_followup_20261027'
    rows = json.loads(previous.read_text())['rows']
    updates = json.loads((followup/'construct_adjudication_final.json').read_text())['rows']
    by_key = {(r['pdb'],r['alpha_chain'],r['beta_chain']):r for r in rows}
    by_key.update({(r['pdb'],r['alpha_chain'],r['beta_chain']):r for r in updates})
    readiness = json.loads((followup/'joint_readiness_final.json').read_text())
    representatives = {}
    for unit in readiness['QC_pass_units']:
        row = by_key[(unit['pdb'],unit['alpha_chain'],unit['beta_chain'])]
        representatives.setdefault(row['receptor_id'],row)
    assert len(representatives) == 29
    targets = []
    for rid,row in sorted(representatives.items()):
        targets.append(dict(target_id=rid,anchor=row['pdb'],scope='complete paired variable identity with an existing fixed-QC-passing anchor',
                            sequences={role:row[role]['variable'] for role in ['alpha','beta']},
                            full_boundary_validated=True,anchor_study=row['source_study_key']))
    # Natural AS8.4 passed identity, but not atomic QC. Search partners without waiving its QC failure.
    as84 = next(r for r in updates if r['pdb']=='8cx4')
    targets.append(dict(target_id='AS8.4-discovery',anchor='8cx4',scope='complete natural variable sequences; current coordinate fails fixed atomic gate',
                        sequences={role:as84[role]['variable'] for role in ['alpha','beta']},
                        full_boundary_validated=True,anchor_study=as84['source_study_key']))
    annotation = ROOT/'docs/fingerprint_cross_source_extension_20261024/sources/stcrdab_all'
    for pdb,name in [('7n2s','AS3.1'),('7n2n','AS4.2'),('7n2p','AS4.3')]:
        row = next(r for r in rows if r['pdb']==pdb)
        domains = numbered_domains((annotation/(pdb+'.html')).read_text())
        seqs = {role:domains[row[role+'_chain']]['sequence'] for role in ['alpha','beta']}
        for role,seq in seqs.items():
            assert row[role]['full_construct'].count(seq) == 1
        targets.append(dict(target_id=name+'-truncated-discovery',anchor=pdb,
                            scope='existing deposited numbered fragment used for discovery only; NOT complete identity',
                            sequences=seqs,full_boundary_validated=False,anchor_study=row['source_study_key'],
                            unchanged_quarantine_reasons=row['reasons']))
    unique = sorted({(role,seq) for t in targets for role,seq in t['sequences'].items()})
    queries,requests = [],[]
    for i,(role,sequence) in enumerate(unique):
        name=f'q{i:03d}.json'
        query=dict(query=dict(type='terminal',service='sequence',parameters=dict(
            evalue_cutoff=.1,identity_cutoff=1.,sequence_type='protein',value=sequence)),
            return_type='polymer_entity',request_options=dict(return_all_hits=True))
        queries.append(dict(file=name,role=role,sequence=sequence,sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),query=query))
        requests.append([name,'https://search.rcsb.org/rcsbsearch/v2/query?json='+quote(json.dumps(query,separators=(',',':')))])
    plan=dict(scope='public structural-development discovery; not an evaluation protocol freeze',
              targets=targets,queries=queries,n_QC_anchor_identities=29,n_AS_discovery_targets=4,
              interpretation='100% local alignment hits do not establish complete-domain identity, chain pairing, independent study or atomic QC',
              source_files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in [previous,followup/'construct_adjudication_final.json',followup/'joint_readiness_final.json']},
              evaluation_protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    (args.out/'search_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    (args.out/'requests.json').write_text(json.dumps(requests,indent=2)+'\n')
    print(json.dumps(dict(n_targets=len(targets),n_unique_chain_queries=len(queries)),indent=2))


if __name__=='__main__':
    main()
