"""Download a bounded IEDB experimental-evidence discovery snapshot.

Public API only; no installation or credentials. This query intentionally
ascertains receptor groups with negative calls, so it is not population sampling.
All downloaded responses are preserved with URL, UTC retrieval time and SHA256.
"""
from __future__ import annotations
import argparse,datetime,hashlib,json,time,urllib.parse,urllib.request
from pathlib import Path
from .experimental_evidence import assay_ids

BASE='https://query-api.iedb.org/'

class Snapshot:
    def __init__(self,folder):
        self.folder=folder;folder.mkdir(parents=True,exist_ok=True)
        self.path=folder/'retrieval_manifest.json'
        self.entries=json.loads(self.path.read_text()) if self.path.exists() else []
    def get(self,name,url,*,json_data=True):
        p=self.folder/name;old=next((x for x in self.entries if x['file']==name),None)
        if p.exists():
            raw=p.read_bytes()
            if not old or old['url']!=url or old['sha256']!=hashlib.sha256(raw).hexdigest():raise ValueError('untracked/mismatched cached source '+name)
        else:
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'meta-acdc-evidence-audit/1.0'}),timeout=60) as response:raw=response.read()
                    break
                except Exception:
                    if attempt==2:raise
                    time.sleep(2**attempt)
            if json_data:json.loads(raw)
            else:
                from xml.etree.ElementTree import fromstring
                fromstring(raw)
            p.write_bytes(raw)
            self.entries.append({'file':name,'url':url,'retrieved_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)})
            temp=self.path.with_suffix('.tmp');temp.write_text(json.dumps(self.entries,indent=2));temp.replace(self.path)
        return json.loads(raw) if json_data else raw
    def page(self,name,endpoint,params,limit=1000):
        result=[];offset=0
        while True:
            url=BASE+endpoint+'?'+urllib.parse.urlencode(dict(params,limit=str(limit),offset=str(offset)))
            rows=self.get(f'{name}_page_{offset//limit}.json',url)
            if not isinstance(rows,list):raise ValueError('API did not return records')
            result.extend(rows)
            if len(rows)<limit:break
            offset+=limit
        return result


def acquire(output):
    s=Snapshot(output)
    candidates=s.page('negative_paired_candidates','tcr_search',{'receptor_type':'eq.alphabeta','chain1_cdr3_seq':'not.is.null','chain2_cdr3_seq':'not.is.null','qualitative_measures':'cs.{Negative}','receptor_species_names':'eq.Homo sapiens (human)',
        'select':'receptor_group_id,pdb_ids,iedb_assay_ids,reference_ids,pubmed_ids,qualitative_measures,receptor_names,chain1_cdr3_seq,chain2_cdr3_seq,assay_names','order':'receptor_group_id'})
    ids=sorted(x['receptor_group_id'] for x in candidates);exports=[]
    for i in range(0,len(ids),20):
        url=BASE+'tcr_export?'+urllib.parse.urlencode({'receptor_group_id':'in.('+','.join(map(str,ids[i:i+20]))+')','limit':'10000','order':'receptor_group_id,receptor__iedb_receptor_id,epitope__iedb_iri,assay__iedb_ids'})
        rows=s.get(f'tcr_export_batch_{i//20}.json',url)
        if len(rows)>=10000:raise ValueError('receptor response may be truncated; reduce batch size')
        exports.extend(rows)
    assays=sorted({aid for row in exports if row.get('assay__type')=='T cell' for aid in assay_ids(row.get('assay__iedb_ids'))})
    for endpoint,key,select in [('tcell_export','assay_id',None),('tcell_search','tcell_id','tcell_id,linear_sequence,structure_id,reference_id,pubmed_id,mhc_allele_name,mhc_restriction,mhc_class,qualitative_measure,assay_names,receptor_ids,tcr_receptor_group_ids')]:
        for i in range(0,len(assays),80):
            params={key:'in.('+','.join(map(str,assays[i:i+80]))+')','limit':'10000','order':key}
            if select:params['select']=select
            rows=s.get(f'{endpoint}_batch_{i//80}.json',BASE+endpoint+'?'+urllib.parse.urlencode(params))
            if {x[key] for x in rows}!=set(assays[i:i+80]):raise ValueError('missing/extra assay exports')
    s.get('PMC7283474.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7283474/fullTextXML',json_data=False)
    for pdb in ['6VRM','6VRN','6VQO']:
        for entity in [3,4]:s.get(f'{pdb}_entity_{entity}.json',f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{entity}')
    print(json.dumps({'candidate_groups':len(ids),'receptor_exports':len(exports),'assays':len(assays),'sources':len(s.entries)}))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args();acquire(a.out)
if __name__=='__main__':main()
