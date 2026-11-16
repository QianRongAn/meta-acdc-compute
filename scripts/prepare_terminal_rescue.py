"""Plan exact observed-domain rescue; cached query reuse is not identity admission."""
from pathlib import Path
import hashlib,json,re
from collections import defaultdict
from urllib.parse import parse_qs,urlsplit,urlencode
from meta_acdc.structure.paired_variable_identity import numbered_domains
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_terminal_rescue_20261101'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cache_records():
    roots=['fingerprint_np_boundaries_20261029/broad_search','fingerprint_native_search_20261029','fingerprint_control_closure_20261101/fragment_rescue'];cache={}
    for folder in roots:
        for manifest in sorted((ROOT/'docs'/folder).rglob('*manifest*.json')):
            records=json.loads(manifest.read_text())
            if not isinstance(records,list):continue
            for record in records:
                if not isinstance(record,dict) or str(record.get('status')) not in ['200','204']:continue
                url=record.get('url','')
                if 'search.rcsb.org/rcsbsearch/v2/query?' not in url:continue
                q=json.loads(parse_qs(urlsplit(url).query)['json'][0]);terminal=q.get('query',{});params=terminal.get('parameters',{})
                if terminal.get('service')!='sequence' or params.get('identity_cutoff')!=1.0 or params.get('evalue_cutoff')!=.1 or params.get('sequence_type')!='protein' or q.get('return_type')!='polymer_entity' or q.get('request_options',{}).get('return_all_hits') is not True:continue
                path=manifest.parent/record['file'];raw=path.read_bytes() if path.exists() else b''
                assert hashlib.sha256(raw).hexdigest()==record['sha256']
                if str(record['status'])=='200':
                    result=json.loads(raw);assert result['result_type']=='polymer_entity';assert isinstance(result['result_set'],list)
                else:assert not raw
                cache[params['value']]=dict(file=str(path.relative_to(ROOT)),manifest=str(manifest.relative_to(ROOT)),sha256=record['sha256'],status=str(record['status']),URL=url)
    return cache

def main():
    BASE.mkdir(exist_ok=False);source=ROOT/'docs/fingerprint_positive_controls_20261027/source_balanced_design.json';rows=json.loads(source.read_text())['rows'];groups=defaultdict(list);skips=[];selected=0
    for row in rows:
        if row['status']!='quarantined' or not any(re.search(r'N terminus|J/FR4|terminal numbering',s,re.I) for s in row['reasons']):continue
        selected+=1;pdb=row['pdb'];p=ROOT/'docs/fingerprint_cross_source_extension_20261024/sources/stcrdab_all'/(pdb+'.html')
        if not p.exists():p=ROOT/'docs/fingerprint_cross_source_20261024/sources/stcrdab'/(pdb+'.html')
        try:
            ann=numbered_domains(p.read_text());seqs={}
            for role,kind in [('alpha','VA'),('beta','VB')]:
                domain=ann[row[role+'_chain']];assert domain['kind']==kind;seq=domain['sequence'];assert 80<=len(seq)<=140 and row[role]['full_construct'].count(seq)==1;seqs[role]=seq
            evidence=dict(pdb=pdb,alpha_chain=row['alpha_chain'],beta_chain=row['beta_chain'],annotation_file=str(p.relative_to(ROOT)),annotation_sha256=sha(p),reasons=row['reasons'],source_study=row.get('source_study_key'))
            groups[(seqs['alpha'],seqs['beta'])].append(evidence)
        except (OSError,ValueError,KeyError,AssertionError) as error:skips.append(dict(pdb=pdb,chains=[row['alpha_chain'],row['beta_chain']],reason=str(error) or type(error).__name__))
    cached=cache_records();targets=[];queries=[];requests=[]
    for i,((a,b),evidence) in enumerate(sorted(groups.items())):targets.append(dict(target_id='fragment_'+hashlib.sha256((a+'|'+b).encode()).hexdigest(),sequences=dict(alpha=a,beta=b),evidence=evidence,complete_identity_claimed=False))
    for i,seq in enumerate(sorted({seq for t in targets for seq in t['sequences'].values()})):
        query=dict(query=dict(type='terminal',service='sequence',parameters=dict(evalue_cutoff=.1,identity_cutoff=1.,sequence_type='protein',value=seq)),return_type='polymer_entity',request_options=dict(return_all_hits=True));name=f'q{i:03d}.json';row=dict(file=name,sequence=seq,query=query,cached_response=cached.get(seq))
        if not row['cached_response']:requests.append([name,'https://search.rcsb.org/rcsbsearch/v2/query?'+urlencode({'json':json.dumps(query,separators=(',',':'))})])
        queries.append(row)
    plan=dict(date='2026-11-01',selection='all terminal-quarantined identity rows with both direct external numbered observed domains uniquely present in the declared full construct; 80-140 aa each; no truncation/imputation',source_file=str(source.relative_to(ROOT)),source_sha256=sha(source),n_terminal_quarantined_rows=selected,n_unique_observed_fragment_pairs=len(targets),n_unique_query_sequences=len(queries),n_cached_queries=sum(bool(x['cached_response']) for x in queries),n_new_queries=len(requests),targets=targets,queries=queries,excluded_rows=skips,interpretation='Grouping identical fragment queries only deduplicates discovery; it does not merge biological receptors or establish complete boundaries.',admission_gates_changed=False,protocol_frozen=False,training_started=False,independent_acceptance_accessed=False)
    (BASE/'search_plan.json').write_text(json.dumps(plan,indent=2)+'\n');(BASE/'new_requests.json').write_text(json.dumps(requests,indent=2)+'\n');print({k:v for k,v in plan.items() if k.startswith('n_')});print('EXCLUDED_ROWS',len(skips))
if __name__=='__main__':main()
