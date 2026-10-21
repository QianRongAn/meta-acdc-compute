"""Audit candidate receptor/source connections without assigning folds or scoring models."""
from pathlib import Path
import argparse,collections,hashlib,json,re,sys,unicodedata
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from meta_acdc.data.experimental_evidence import reserved_hits

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def author_key(surname,given):
    clean=lambda x:re.sub('[^a-z]','',unicodedata.normalize('NFKD',x.casefold()))
    a,b=clean(surname or ''),clean(given or '')
    return a+':'+b[0] if a and b else None

def components(nodes,edges):
    parent={x:x for x in nodes}
    def find(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    for e in edges:
        a,b=find(e['a']),find(e['b']);parent[max(a,b)]=min(a,b)
    groups=collections.defaultdict(list)
    for x in sorted(nodes):groups[find(x)].append(x)
    return sorted(groups.values())

def audit(root):
    d=root/'docs/measured_intake_20261021';primary=load(d/'primary_observations.json');secondary=load(root/'docs/six_blockers_20261021/batcave_candidate_index.json')
    source_paths=['docs/measured_intake_20261021/primary_observations.json','docs/measured_intake_20261021/secondary_reconciliation.json','docs/six_blockers_20261021/batcave_candidate_index.json','docs/independent_spr_20261018/development_study_sources.json']
    nodes=[]
    for x in primary['receptors']:
        y=dict(x);y.update(origin='primary_PTEAM',sequence_pair_review_pass=x['DNA_pair_passes'] and not x['author_excluded']);nodes.append(y)
    for x in secondary['constructs']:
        nodes.append(dict(receptor_id='BATCAVE:'+x['construct_id'],receptor_name=x['tcr_name'],pmid=x['pmid'],origin='secondary_BATCAVE',alpha_variable=x['alpha_published'] or '',beta_variable=x['beta_published'] or '',cdr3_alpha=x['cdr3_alpha'] or '',cdr3_beta=x['cdr3_beta'] or '',sequence_pair_review_pass=x['canonical_pair_with_CDR3'],source_file=x['source_file'],source_rows=x['source_rows'],assay=x['assay'],training_allowed=False,model_scored=False))
    studies=load(root/'docs/independent_spr_20261018/development_study_sources.json');authors={}
    for s in studies:
        authors[s['pmid']]={k for a in s['authors']['author'] if (k:=author_key(a.get('lastName'),a.get('firstName') or a.get('initials')))}
    for path in sorted((d/'source_metadata').glob('*.json')):
        raw=load(path)['resultList']['result'][0]
        if raw['id']!=path.stem or raw['source']!='MED':raise ValueError('source PMID identity mismatch')
        rows=raw.get('authorList',{}).get('author',[])
        keys={k for a in rows if (k:=author_key(a.get('lastName'),a.get('firstName') or a.get('initials')))}
        if not keys:raise ValueError('source author metadata missing')
        authors[raw['id']]=keys;source_paths.append(str(path.relative_to(root)))
    source_paths.append('docs/measured_intake_20261021/source_metadata_retrieval.json')
    authors['39151427']={k for a in primary['authors'] if (k:=author_key(a['surname'],a['given_names']))}
    assert authors['39151427'],'primary authors missing'
    registries={}
    for label,path in [('old_affinity_development','docs/database_binding/affinity_dataset.json'),('previous_reservations_including_MIRA','docs/database_binding/combined_reservation_registry.json'),('clinical_reservation','docs/clinical_primary_20261018/clinical_reservations.json'),('consumed_SPR_test_identity','docs/independent_spr_20261018/frozen_test.json'),('fresh_SPR_reservation','docs/fresh_spr_reserve_20261018/reservation_registry.json')]:
        source_paths.append(path);raw=load(root/path);raw=raw if isinstance(raw,list) else raw['receptors']
        registries[label]=[{k:x.get(k,'') for k in ['receptor_name','alpha_variable','beta_variable','cdr3_alpha','cdr3_beta','pmid']} for x in raw]
    exposure=[]
    for x in nodes:
        hits={label:reserved_hits(x,refs) for label,refs in registries.items()}
        source_hits=[dict(development_pmid=s['pmid'],shared_author_keys=sorted(authors.get(x['pmid'],set())&authors[s['pmid']])) for s in studies if authors.get(x['pmid'],set())&authors[s['pmid']]]
        exposure.append(dict(receptor_id=x['receptor_id'],registry_hits=hits,development_author_overlap=source_hits,source_author_metadata_complete=x['pmid'] in authors,formal_independent_acceptance_allowed=False,training_allowed=False))
    edges=[];near=[]
    for i,a in enumerate(nodes):
        for b in nodes[:i]:
            reasons=[]
            if a['pmid']==b['pmid']:reasons.append('same_primary_study')
            shared=sorted(authors.get(a['pmid'],set())&authors.get(b['pmid'],set()))
            if a['pmid']!=b['pmid'] and shared:reasons.append('shared_author_conservative')
            if reserved_hits(a,[b]):
                reasons.append('corresponding_chain_or_CDR3_near_match')
                near.append(dict(a=a['receptor_id'],b=b['receptor_id']))
            if reasons:edges.append(dict(a=a['receptor_id'],b=b['receptor_id'],reasons=reasons,shared_author_keys=shared if a['pmid']!=b['pmid'] else []))
    graph=components([x['receptor_id'] for x in nodes],edges);families=components([x['receptor_id'] for x in nodes],near)
    primary_ids={x['receptor_id'] for x in primary['receptors']}
    pblocks=[g for g in graph if primary_ids&set(g)]
    summary=dict(nodes=len(nodes),primary_receptors=len(primary_ids),secondary_constructs=len(secondary['constructs']),primary_sequence_pass=sum(x['DNA_pair_passes'] and not x['author_excluded'] for x in primary['receptors']),primary_observations=len(primary['observations']),primary_sequence_pass_observations=sum(x['sequence_review_pass'] for x in primary['observations']),primary_studies=len({x['pmid'] for x in primary['receptors']}),primary_joint_blocks=len(pblocks),provisional_joint_blocks=len(graph),provisional_sequence_components=len(families),primary_receptors_with_prior_reservation_hits=sum(any(e['registry_hits'].values()) for e in exposure if e['receptor_id'] in primary_ids),primary_receptors_with_development_author_overlap=sum(bool(e['development_author_overlap']) for e in exposure if e['receptor_id'] in primary_ids),unresolved_secondary_author_pmids=sorted({x['pmid'] for x in nodes if x['pmid'] not in authors}),new_frozen_test_created=False,training_release_created=False,model_scored=False)
    exposures={x['receptor_id']:x for x in exposure}
    block_review=[dict(members=g,prior_family_exposure=any(any(exposures[x]['registry_hits'].values()) for x in g),development_author_overlap=any(exposures[x]['development_author_overlap'] for x in g),unresolved_source_metadata=any(not exposures[x]['source_author_metadata_complete'] for x in g),independent_acceptance_allowed=False,training_allowed=False) for g in graph]
    return dict(block_review=block_review,rule='Existing reserved_hits: conservative corresponding-chain or CDR3 near-match, plus same primary PMID or shared surname/first-initial author key; transitive connected blocks. Sequence components are exclusion screens, not validated biological families. Missing author metadata never establishes independence.',source_sha256={n:sha(root/n) for n in source_paths},summary=summary,nodes=nodes,edges=edges,provisional_sequence_components=families,provisional_joint_blocks=graph,exposure=exposure,frozen_test_labels_used_for_tuning=False,partition_assignment=None)

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    result=audit(ROOT);a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['summary'],indent=2))
if __name__=='__main__':main()
