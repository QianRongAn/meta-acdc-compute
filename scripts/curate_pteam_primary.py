"""Keep primary research TCR DNA/response data separate from clinical products."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from meta_acdc.data.clinical_primary import CODONS
RAW=ROOT/'data/external_evidence/20261021-six-blockers';OUT=ROOT/'docs/six_blockers_20261021'
def translate(dna):
 dna=re.sub(r'\s+','',str(dna or '')).upper()
 if not dna or set(dna)-set('ACGT'):return dict(dna=dna,valid=False,reason='noncanonical_or_missing_DNA')
 aa=''.join(CODONS[dna[i:i+3]] for i in range(0,len(dna)-2,3));body=aa.removesuffix('*')
 return dict(dna=dna,amino_acids=body,trailing_nt=dna[len(aa)*3:],terminal_stop=aa.endswith('*'),internal_stop='*' in body,valid=(len(dna)%3==0 and '*' not in body))
def main():
 import openpyxl
 target=OUT/'pteam_primary_sequence_index.json'
 if target.exists():raise FileExistsError(target)
 receptors=[];summary={}
 for name in ['pteam_cmv.xlsx','pteam_rnf43.xlsx']:
  w=openpyxl.load_workbook(RAW/name,read_only=True,data_only=True)
  if 'cmv' in name:
   matrix=list(w['Mean'].values);allowed=set(matrix[0][1:]);seqs=list(w['TCR sequences'].values)[1:];items=[]
   for rowno,row in enumerate(seqs,start=2):
    label=str(row[2] or '').removeprefix('TCR').strip().replace('-','_')
    if label not in allowed:continue
    items.append(dict(source_row=rowno,name='TCR'+label.replace('_','-'),matrix_column=label,cdr3_alpha=row[9],cdr3_beta=row[6],alpha=translate(row[10]),beta=translate(row[11]),alpha_published_aa=row[12],beta_published_aa=row[13],scope='published_contiguous_research_chain_DNA',peptide=row[1],mhc=row[0]))
  else:
   matrix=list(w['Normalized by PC'].values);allowed=set(x for x in matrix[0][4:11] if x);seqs=list(w['Sequences'].values)[2:];items=[]
   for rowno,row in enumerate(seqs,start=3):
    if row[0] not in allowed:continue
    items.append(dict(source_row=rowno,name=row[0],matrix_column=row[0],cdr3_alpha=row[1],cdr3_beta=row[2],alpha=translate(row[10]),beta=translate(row[12]),alpha_constant_DNA=row[11],beta_constant_DNA=row[13],alpha_published_aa=None,beta_published_aa=None,scope='published_variable_leader_DNA_with_separate_constant_segments_not_assembled_here',peptide=row[3],mhc=row[4]))
  observations=0
  for x in items:
   for chain in ['alpha','beta']:
    x[chain]['cdr3_matches']=bool(x['cdr3_'+chain]) and x['cdr3_'+chain] in x[chain].get('amino_acids','')
    expected=x[chain+'_published_aa']
    x[chain]['published_translation_matches']=(x[chain].get('amino_acids')==str(expected).removesuffix('*')) if expected else None
   x['DNA_pair_passes']=all(x[c]['valid'] and x[c]['cdr3_matches'] and x[c]['published_translation_matches'] is not False for c in ['alpha','beta'])
   x['author_excluded']=x['name']=='R27';x['source_file']=name;x['pmid']='39151427';x['clinical_administered_product_link']=None;x['clinical_risk_label']=None;x['clinical_training_allowed']=False;x['model_scored']=False
   cols=[i for i,v in enumerate(matrix[0]) if v==x['matrix_column']]
   if len(cols)!=1:raise ValueError('Ambiguous response column')
   col=cols[0];values=[dict(source_excel_row=i,value=float(row[col])) for i,row in enumerate(matrix[1:],start=2) if isinstance(row[col],(float,int))]
   x['functional_observations']=values;x['functional_endpoint']='published_normalized_T_cell_activation_not_SPR_KD';x['molecular_binding_label']=None;x['training_allowed']=False;observations+=len(values)
  summary[name]=dict(receptors=len(items),DNA_pairs_passing=sum(x['DNA_pair_passes'] for x in items),author_excluded_receptors=sum(x['author_excluded'] for x in items),observations=observations,eligible_for_further_primary_review=sum(len(x['functional_observations']) for x in items if x['DNA_pair_passes'] and not x['author_excluded']),source_sha256=hashlib.sha256((RAW/name).read_bytes()).hexdigest())
  receptors.extend(items)
 data=dict(status='unscored_research_DNA_and_functional_evidence_candidates',pmid='39151427',doi='10.1016/j.xgen.2024.100634',source='author-deposited Zenodo 11197942, primary workbook sequence and response sheets',summary=summary,receptors=receptors,clinical_validation_passed=False,independence_not_assumed=True,training_release_created=False)
 target.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
