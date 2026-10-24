"""Metadata-only source-balanced retrieval eligibility; no outcomes or shuffled splits."""
import re,unicodedata
from collections import defaultdict

def author_keys(row):
    authors=(row.get('primary_citation') or {}).get('rcsb_authors') or []
    keys={re.sub(r'[^a-z0-9]','',unicodedata.normalize('NFKD',x).casefold()) for x in authors}
    return keys-{''}

def representatives(rows):
    groups=defaultdict(list)
    for row in rows:
        if row.get('source_study_key'):groups[(row['receptor_id'],row['source_study_key'])].append(row)
    return [sorted(items,key=lambda x:(x['pdb'],x['alpha_chain'],x['beta_chain']))[0] for _,items in sorted(groups.items())]

def source_balanced_blocks(rows,disjoint_authors=False):
    records=representatives(rows);blocks=[]
    for query in records:
        studies=sorted({x['source_study_key'] for x in records if x['source_study_key']!=query['source_study_key'] and x['V_context']==query['V_context']})
        for study in studies:
            gallery=[x for x in records if x['source_study_key']==study and x['pdb']!=query['pdb'] and x['V_context']==query['V_context']]
            positive=[x for x in gallery if x['receptor_id']==query['receptor_id']];negative=[x for x in gallery if x['receptor_id']!=query['receptor_id']]
            if not positive or not negative:continue
            qauthors=author_keys(query);gauthors=set().union(*(author_keys(x) for x in gallery));independent=bool(qauthors and gauthors and not qauthors&gauthors)
            if disjoint_authors and not independent:continue
            blocks.append(dict(query_pdb=query['pdb'],query_receptor_id=query['receptor_id'],query_study=query['source_study_key'],gallery_study=study,V_context=query['V_context'],positive_pdbs=[x['pdb'] for x in positive],negative_pdbs=[x['pdb'] for x in negative],negative_receptor_ids=[x['receptor_id'] for x in negative],author_disjoint=independent))
    return blocks

def joint_components(rows):
    parent=list(range(len(rows)));seen={}
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,row in enumerate(rows):
        keys=[('receptor',row['receptor_id']),('V_context',row['V_context'])]
        if row.get('source_study_key'):keys.append(('study',row['source_study_key']))
        for role in ['alpha','beta']:
            if row.get(role,{}).get('variable'):keys.append((role,row[role]['variable']))
        for key in keys:
            if key in seen:parent[find(i)]=find(seen[key])
            else:seen[key]=i
    groups=defaultdict(list)
    for i,row in enumerate(rows):groups[find(i)].append(row)
    return list(groups.values())
