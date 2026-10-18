"""Fail-closed validation of a new, unscored, currently incomplete SPR reserve.

No scoring or training entry point. Old consumed test records are never loaded.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from difflib import SequenceMatcher
from .release import sha256
from .experimental_evidence import reserved_hits
from .independent_spr import author_keys,source_components,connected_test_blocks,validate_partition
ROOT=Path(__file__).resolve().parents[3]
RESERVE=ROOT/'docs/fresh_spr_reserve_20261018'

def read(name):return json.loads((RESERVE/name).read_text())

def assert_development_safe(rows, reserved=None, studies=None):
    reserved=read('reserved_measurements.json') if reserved is None else reserved
    studies=read('study_sources.json') if studies is None else studies
    by_dev={s['pmid']:s for s in json.loads((ROOT/'docs/independent_spr_20261018/development_study_sources.json').read_text())}
    keys=set().union(*(author_keys(s) for s in studies))
    for row in rows:
        if row['pmid'] not in by_dev or author_keys(by_dev[row['pmid']])&keys:raise ValueError('reserved study source in development')
        if reserved_hits(row,reserved):raise ValueError('reserved receptor family in development')
        if any(SequenceMatcher(None,row['peptide'],r['peptide'],autojunk=False).ratio()>=.8 for r in reserved):raise ValueError('reserved peptide family in development')
    return True

def verify(*,check_raw=True):
    seal=read('manifest.json')
    for n,h in seal['artifacts'].items():
        if sha256(ROOT/n)!=h:raise ValueError('reserve artifact changed: '+n)
    for n,h in seal['reference_files'].items():
        if sha256(ROOT/n)!=h:raise ValueError('reserve reference changed: '+n)
    for n,h in seal['code_files'].items():
        if sha256(ROOT/n)!=h:raise ValueError('reserve code changed: '+n)
    if check_raw:
        for n,h in read('required_sources.json').items():
            if sha256(ROOT/n)!=h:raise ValueError('primary source changed: '+n)
        for n,meta in read('raw_sources.json').items():
            if sha256(ROOT/n)!=meta['sha256']:raise ValueError('search artifact changed: '+n)
    test=read('reserved_measurements.json');studies=read('study_sources.json');protocol=read('protocol.json');access=read('data_access.json')
    dev_path=ROOT/protocol['development_path']
    if sha256(dev_path)!=seal['development_sha256']:raise ValueError('source-safe development changed')
    dev=json.loads(dev_path.read_text())
    assert_development_safe(dev,test,studies)
    parent=json.loads((ROOT/'docs/sequence_followup_20261018/review_ledger.json').read_text())
    # The previous source-disjoint20 is the authorized parent. Older excluded
    # studies stay excluded; do not reinstate them because this is a new source.
    plan=json.loads((ROOT/'docs/affinity_development_20261018/plan_v2.json').read_text())
    parent=[x for x in parent if x['evidence_id'] in plan['parent_record_sha256']]
    dev_sources=json.loads((ROOT/'docs/independent_spr_20261018/development_study_sources.json').read_text())
    # Reuse every existing biological, source, parent-membership and range gate.
    # Only the final cohort-size failure is expected for this partial reservation.
    try:validate_partition(test,dev,studies,dev_sources,parent)
    except ValueError as exc:
        if str(exc)!='insufficient independent blocks':raise
    else:raise ValueError('partial protocol must be replaced before formal acceptance')
    counts=dict(observations=len(test),primary_observations=sum(x['metric_role']=='primary' for x in test),sensitivity_observations=sum(x['metric_role']=='sensitivity_only' for x in test),lineages=len({x['lineage'] for x in test}),studies=len(studies),research_source_components=len(source_components(studies)),connected_blocks=len(connected_test_blocks(test)))
    if counts!=protocol['counts']:raise ValueError('reserve counts changed')
    if check_raw:
        from .evidence_completion import structure_registry
        refs=[json.loads((ROOT/'docs/database_binding/affinity_dataset.json').read_text()),json.loads((ROOT/'docs/clinical_primary_20261018/clinical_reservations.json').read_text())['receptors'],json.loads((ROOT/'docs/independent_spr_20261018/reservation_registry.json').read_text())['receptors']]
        refs.extend(structure_registry(json.loads(p.read_text())) for p in (ROOT/'data/releases').glob('*/structure_records.json'))
        for x in test:
            if any(reserved_hits(x,reg) for reg in refs):raise ValueError('historically exposed reserve family')
        identity=read('sequence_identity_audit.json');chains={}
        raw=ROOT/'data/external_evidence/20261018-fresh-spr-v1/9K2R.fasta'
        for rec in raw.read_text().split('>')[1:]:
            header,*lines=rec.splitlines()
            for c in ('alpha','beta'):
                if 'receptor' in header.lower() and c in header.lower():chains[c]=''.join(lines)
        for c,motif in [('alpha','DIQNP'),('beta','EDLNK')]:
            if chains[c]!=identity['full_'+c] or chains[c].split(motif)[0]!=identity[c+'_variable']:raise ValueError('deposited variable sequence mismatch')
            if any(x[c+'_variable']!=identity[c+'_variable'] for x in test):raise ValueError('reserve variable mismatch')
    if protocol['evaluation_allowed'] is not False or access['model_scored'] is not False or access['number_of_model_predictions']!=0:raise ValueError('partial reserve cannot be scored')
    if protocol['status']!='partial_reserved_not_formal_acceptance':raise ValueError('unexpected reserve state')
    return dict(status=protocol['status'],evaluation_allowed=False,training_allowed=False,model_scored=False,counts=counts,source_safe_development_rows=len(dev),inherited_anchor_has_configuration_exposure=True)

def require_acceptance_ready():
    result=verify()
    if not result['evaluation_allowed']:raise ValueError('new acceptance blocked: fewer than three source components/connected blocks; retained anchor also has inherited configuration exposure')
    return result

def main():print(json.dumps(verify(),indent=2))
if __name__=='__main__':main()
