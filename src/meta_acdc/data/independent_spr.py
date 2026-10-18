"""Validate a frozen external SPR panel and its source-disjoint development subset.

This is a data contract, not a training/scoring entry point. Old datasets/folds are
never reassigned. Repeated measurements, engineered variants, and shared labs
must not inflate the independent-sample count.
"""
from __future__ import annotations
import argparse, json, math, re, unicodedata
from collections import Counter
from difflib import SequenceMatcher
from .build_release import similar
from pathlib import Path
from .release import sha256
from .experimental_evidence import reserved_hits
from .evidence_completion import structure_registry

ROOT=Path(__file__).resolve().parents[3]
FREEZE=ROOT/'docs/independent_spr_20261018'
AA=set('ACDEFGHIKLMNPQRSTVWY')


def author_keys(study):
    authors=study.get('authors',{}).get('author',[])
    if not authors or any(not x.get('fullName') for x in authors):
        raise ValueError('missing source authors')
    # Full normalized surname+initials; never match surnames alone. Ambiguous
    # homonyms conservatively exclude rather than justify independence.
    keys=set()
    for author in authors:
        name=unicodedata.normalize('NFKD',author['fullName'].casefold())
        keys.add(re.sub('[^a-z0-9]','',name))
        identifier=author.get('authorId') or {}
        if identifier.get('type')=='ORCID' and identifier.get('value'):
            keys.add('orcid:'+identifier['value'].rstrip('/').split('/')[-1].casefold())
    return keys


def source_components(studies):
    if len({s['pmid'] for s in studies})!=len(studies):raise ValueError('duplicate source study')
    groups=[]
    for s in studies:
        keys=author_keys(s);joined=[g for g in groups if g[1]&keys]
        pmids={s['pmid']}
        for g in joined:
            groups.remove(g);pmids|=g[0];keys|=g[1]
        groups.append((pmids,keys))
    return sorted(sorted(x[0]) for x in groups)


def connected_test_blocks(rows):
    """Original study/near-receptor/near-peptide block rule, without fold assignment."""
    parent=list(range(len(rows)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,a in enumerate(rows):
        for j,b in enumerate(rows[:i]):
            if (a['pmid']==b['pmid'] or similar(a['alpha_variable'],b['alpha_variable'])
                or similar(a['beta_variable'],b['beta_variable'])
                or SequenceMatcher(None,a['peptide'],b['peptide'],autojunk=False).ratio()>=.8):
                parent[find(i)]=find(j)
    groups={}
    for i,x in enumerate(rows):groups.setdefault(find(i),[]).append(x['evidence_id'])
    return sorted(sorted(g) for g in groups.values())


def validate_partition(test, development, test_sources, development_sources, parent_review):
    """Recompute source/family disjointness and immutable parent membership."""
    by_test={x['pmid']:x for x in test_sources};by_dev={x['pmid']:x for x in development_sources}
    if {x['pmid'] for x in test}!=set(by_test):raise ValueError('test/source coverage')
    test_keys=set().union(*(author_keys(x) for x in test_sources))
    parent={x['evidence_id']:x for x in parent_review if x['training_allowed']}
    if len({x['evidence_id'] for x in test+development})!=len(test)+len(development):raise ValueError('duplicate evidence or fold overlap')
    for x in test:
        if (x.get('fold')!='external_test' or x.get('training_allowed') is not False
            or x.get('model_scored') is not False or any(x.get(k) is not None for k in ('training_label','binary_binding_label','clinical_risk_label'))):
            raise ValueError('external endpoint promotion')
        if (x.get('endpoint')!='tcr_binding_spr' or x.get('species')!='Homo sapiens'
            or x.get('receptor_type')!='alphabeta' or not re.fullmatch(r'HLA-[ABC]\*\d{2,3}:\d{2,3}',x.get('mhc_allele',''))):
            raise ValueError('wrong biological endpoint')
        if x.get('raw_unit') not in ('uM','nM') or x.get('qualifier')!='=':raise ValueError('unknown units or censored value')
        v=float(x['raw_value'])*(1000 if x['raw_unit']=='uM' else 1)
        if not math.isfinite(v) or v<=0 or not math.isclose(v,x['kd_nm'],rel_tol=1e-12):raise ValueError('KD conversion')
        if not math.isclose(x['pkd'],9-math.log10(v),rel_tol=1e-12):raise ValueError('pKD conversion')
        if not 273<=x.get('temperature_K',0)<=320 or x.get('fit') not in ('equilibrium','kinetic') or 'Biacore' not in x.get('instrument',''):raise ValueError('SPR conditions')
        for chain in ('alpha','beta'):
            seq=x.get(chain+'_variable','');cdr=x.get('cdr3_'+chain,'')
            if len(seq)<70 or len(cdr)<8 or cdr not in seq or set(seq)-AA:raise ValueError('paired variable identity')
        if not x.get('peptide') or set(x['peptide'])-AA:raise ValueError('noncanonical peptide')
        for k in ('sources','sequence_sources','conditions','primary_locator','lineage','research_source'):
            if not x.get(k):raise ValueError('primary evidence missing '+k)
        if x.get('metric_role') not in ('primary','sensitivity_only'):raise ValueError('undefined metric subset')
        top=x.get('max_tested_um')
        if top is not None and v/1000>top and x['metric_role']=='primary':raise ValueError('range extrapolation in primary metric')
    for x in development:
        p=parent.get(x['evidence_id'])
        if p is None or p['parent_fold'] not in ('train','validation'):raise ValueError('unapproved development record')
        if (x['fold'],x['split_block'])!=(p['parent_fold'],p['parent_block']):raise ValueError('old fold or block reassigned')
        # Changes to annotations cannot silently alter identity or original labels.
        for k,v in p['original_record'].items():
            if k not in ('primary_review','temperature_K') and x.get(k)!=v:raise ValueError('immutable parent field '+k)
        if x['pmid'] not in by_dev or author_keys(by_dev[x['pmid']])&test_keys:raise ValueError('research source overlap')
        if reserved_hits(x,test):raise ValueError('test receptor family in development')
        if any(SequenceMatcher(None,x['peptide'],y['peptide'],autojunk=False).ratio()>=.8 for y in test):raise ValueError('test peptide family in development')
    # Exact coverage: no outcome/score-dependent dropping of otherwise eligible rows.
    expected=set()
    for eid,p in parent.items():
        if p['parent_fold'] not in ('train','validation'):continue
        if p['pmid'] not in by_dev:raise ValueError('missing development source')
        if not author_keys(by_dev[p['pmid']])&test_keys and not reserved_hits(p['original_record'],test):expected.add(eid)
    if expected!={x['evidence_id'] for x in development}:raise ValueError('development filter coverage')
    groups=source_components(test_sources)
    counts=Counter(x['fold'] for x in development)
    blocks={f:len({x['split_block'] for x in development if x['fold']==f}) for f in ('train','validation')}
    test_blocks=connected_test_blocks(test)
    if min(blocks.values())<3 or len(groups)<3 or len(test_blocks)<3:raise ValueError('insufficient independent blocks')
    return dict(observations=len(test),primary=sum(x['metric_role']=='primary' for x in test),lineages=len({x['lineage'] for x in test}),studies=len(test_sources),research_source_components=len(groups),development=dict(counts),development_blocks=blocks,test_connected_blocks=len(test_blocks))


def verify_freeze(path=FREEZE, *, check_raw=True):
    path=Path(path)
    read=lambda n:json.loads((path/n).read_text())
    seal=read('manifest.json')
    if seal['model_scored'] is not False or seal['optimizer_started'] is not False:raise ValueError('freeze exposure state')
    for group in ('artifacts','code_files','previous_freezes','reference_files'):
        for n,digest in seal[group].items():
            p=(ROOT/n).resolve()
            if not p.is_relative_to(ROOT) or sha256(p)!=digest:raise ValueError('freeze hash mismatch: '+n)
    if check_raw:
        for n,digest in read('required_sources.json').items():
            p=(ROOT/n).resolve()
            if not p.is_relative_to(ROOT) or sha256(p)!=digest:raise ValueError('primary source changed: '+n)
        for n,meta in read('raw_sources.json').items():
            if sha256(ROOT/n)!=meta['sha256']:raise ValueError('raw search source changed: '+n)
    test=read('frozen_test.json');dev=read('development_records.json')
    result=validate_partition(test,dev,read('test_study_sources.json'),read('development_study_sources.json'),json.loads((ROOT/'docs/sequence_followup_20261018/review_ledger.json').read_text()))
    # Independently recompute historical family exposure, not trust a curated zero.
    old=json.loads((ROOT/'docs/database_binding/affinity_dataset.json').read_text())
    historical=[old]+[structure_registry(json.loads(p.read_text())) for p in sorted((ROOT/'data/releases').glob('*/structure_records.json'))]
    clinical=json.loads((ROOT/'docs/clinical_primary_20261018/clinical_reservations.json').read_text())['receptors']
    for x in test:
        if any(reserved_hits(x,refs) for refs in historical+[clinical]):raise ValueError('historically exposed receptor in frozen test')
    # Bind transcribed sequence and deposited chains, never infer a paired chain.
    for x in read('primary_sequence_audit.json'):
        if not x['no_imputation']:raise ValueError('imputed sequence')
        if x['pdb']:
            raw=ROOT/'data/external_evidence/20261018-independent-spr-v1'/(x['pdb']+'.fasta')
            chains={}
            if check_raw:
                for record in raw.read_text().split('>')[1:]:
                    h,*lines=record.splitlines()
                    for chain in ('alpha','beta'):
                        if 'TCR '+chain in h and 'Homo sapiens (9606)' in h:chains[chain]=''.join(lines)
                if chains!={c:x['full_'+c] for c in ('alpha','beta')}:raise ValueError('primary/deposited sequence conflict')
    counts=read('protocol.json')['counts']
    if (counts['test_observations'],counts['primary_metric_observations'],counts['lineages'],counts['research_source_components'])!=(result['observations'],result['primary'],result['lineages'],result['research_source_components']):raise ValueError('protocol count mismatch')
    result.update(training_allowed=False,model_scored=False,optimizer_started=False)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--freeze',type=Path,default=FREEZE);a=p.parse_args()
    print(json.dumps(verify_freeze(a.freeze),indent=2))

if __name__=='__main__':main()
