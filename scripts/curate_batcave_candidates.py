from pathlib import Path
import collections,hashlib,json,openpyxl,re,math
ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'data/external_evidence/20261021-six-blockers';OUT=ROOT/'docs/six_blockers_20261021'
AA=set('ACDEFGHIKLMNPQRSTVWY')
def seq(value):
 s=re.sub(r'\s+','',str(value or '')).upper()
 return s if len(s)>=70 and set(s)<=AA else None
def usable_pair(row):
 a,b=seq(row.get('va')),seq(row.get('vb'))
 ca,cb=str(row.get('cdr3a') or ''),str(row.get('cdr3b') or '')
 return bool(a and b and ca and cb and ca in a and cb in b)
def sheet(name):
 w=openpyxl.load_workbook(RAW/name,read_only=True,data_only=True);s=w.active;it=iter(s.values);keys=next(it)
 return [(i,dict(zip(keys,row))) for i,row in enumerate(it,start=2)]
def main():
 output=OUT/'batcave_candidate_index.json'
 if output.exists():raise FileExistsError(output)
 counts={};constructs={};missing=collections.Counter();allrows=[]
 for name in ['batcave_mhci.xlsx','batcave_mhcii.xlsx']:
  rows=sheet(name);hum=[(i,x) for i,x in rows if str(x['tcr_source_organism']).casefold()=='human'];valid=[(i,x) for i,x in hum if usable_pair(x)]
  counts[name]=dict(rows=len(rows),human_rows=len(hum),human_canonical_paired_sequence_rows=len(valid),human_zero_activity_rows=sum(x['peptide_activity']==0 for _,x in hum),human_valid_pair_zero_activity_rows=sum(x['peptide_activity']==0 for _,x in valid),source_sha256=hashlib.sha256((RAW/name).read_bytes()).hexdigest())
  for i,x in rows:
   if str(x['tcr_source_organism']).casefold()!='human':missing['nonhuman_functional_rows']+=1;continue
   if not usable_pair(x):missing['missing_noncanonical_or_CDR3_inconsistent_paired_sequence']+=1
   key=(str(x['pmid']),x['tcr'],x['assay'],x['mhc'],seq(x.get('va')),seq(x.get('vb')))
   ident=hashlib.sha256(json.dumps(key).encode()).hexdigest()[:20]
   if ident not in constructs:
    constructs[ident]=dict(construct_id=ident,tcr_name=x['tcr'],pmid=str(x['pmid']),assay=x['assay'],mhc=x['mhc'],index_peptide=x['index_peptide'],species='human',source_file=name,source_rows=[],alpha_published=seq(x.get('va')),beta_published=seq(x.get('vb')),cdr3_alpha=x['cdr3a'],cdr3_beta=x['cdr3b'],canonical_pair_with_CDR3=usable_pair(x),clinical_alias_requires_reservation=x['tcr'] in {'a3a','c259','TIL1383I','DMF4','DMF5'},physical_binding_label=None,clinical_risk_label=None,full_clinical_engineered_DNA_verified=False,training_allowed=False,model_scored=False,admission='primary_assay_protocol_and_identity_review_required')
   constructs[ident]['source_rows'].append(i)
 kd=sheet('batcave_kd.xlsx');human_names={'c259','TIL1383I','868Z11'}
 kdsummary=dict(rows=len(kd),censored_rows=sum(x['Kd_is_greater']==1 for _,x in kd),censored_by_receptor=dict(collections.Counter(x['tcr'].strip() for _,x in kd if x['Kd_is_greater']==1)),human_rows=sum(x['tcr'].strip() in human_names for _,x in kd),human_censored_rows=sum(x['tcr'].strip() in human_names and x['Kd_is_greater']==1 for _,x in kd),human_identities='868Z11, c259, TIL1383I per deposited table; clinical construct equivalence not assumed',all_censored_rows_nonhuman=True,training_allowed=False,reason='secondary curated/possibly digitized values; original SPR conditions and numeric provenance needed; preserve Kd_is_greater, never turn bound into exact KD')
 data=dict(purpose='source-indexed candidate intake, not a training release',source_repository='meyer-lab-cshl/BATMAN-paper',source_commit='37d8d1ba23921d90736418b7b6cd62a5cdf68a91',counts=counts,quarantine_counts=dict(missing),constructs=list(constructs.values()),affinity_candidate_summary=kdsummary,zero_activity_is_universal_nonbinding=False,clinical_validation_passed=False,all_new_candidates_unscored=True)
 output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps(dict(counts=counts,kd=kdsummary,constructs=len(constructs)),indent=2))
if __name__=='__main__':main()
