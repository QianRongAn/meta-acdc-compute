from pathlib import Path
import json,sys,csv
from collections import Counter
import argparse
r=Path(__file__).resolve().parents[1];sys.path.insert(0,str(r/'src'))
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args()
args.out.mkdir(parents=True,exist_ok=False)
from meta_acdc.data.build_release import domain
from meta_acdc.data.release import sha256
from meta_acdc.data.primary_panel import validate_measurements
from meta_acdc.data.experimental_evidence import write_tsv
exp=r/'data/external_evidence/20261016-expansion';old=r/'data/external_evidence/20261016-repro';inputs=r/'docs/evidence';d=args.out
def src(path,url,locator):return {'path':str(path.relative_to(r)),'sha256':sha256(path),'url':url,'locator':locator}
def seq_sources(folder,pdb,ids):return [src(folder/f'{pdb}_entity_{i}.json',f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{i}','entity_poly.pdbx_seq_one_letter_code_can') for i in ids]
rows=[]
for x in json.loads((inputs/'p53_primary_curation.json').read_text()):
 ss=json.loads(x['receptor_sequence_source']);ss=[src(old/a['file'],a['url'],'entity_poly.pdbx_seq_one_letter_code_can') for a in ss]
 rows.append(dict(evidence_id=x['evidence_id'],receptor_name=x['receptor_name'],lineage=x['receptor_name'],species='Homo sapiens',receptor_type='alphabeta',construct_notes=x['construct_notes'],alpha_variable=x['alpha_variable'],beta_variable=x['beta_variable'],sequence_sources=ss,pmid=x['pmid'],primary_locator=x['primary_locator'],sources=[src(old/x['primary_file'],x['primary_url'],x['primary_locator'])],mhc_allele=x['mhc_allele'],peptide=x['peptide'],endpoint='tcr_binding_spr',observed_outcome='measured_binding' if x['label'] else 'not_detected',raw_value=x['raw_value'] if x['label'] else 'NB',raw_unit=x['raw_unit'] if x['label'] else None,qualifier='=' if x['label'] else None,kd_nm=x['kd_nm'],conditions=x['conditions'],max_tested_um=x['negative_max_tcr_um'],condition_limit='SPR non-detection is conditional; not a cellular or clinical safety label',training_label=None,use='reserved_external_only'))
def receptor(name,pdb,ids,pmid,pmc,mhc,lineage,conditions,notes,locator):
 ss=seq_sources(exp,pdb,ids);seqs=[domain(json.loads((r/a['path']).read_text())['entity_poly']['pdbx_seq_one_letter_code_can'])[0] for a in ss]
 return dict(receptor_name=name,lineage=lineage,species='Homo sapiens',receptor_type='alphabeta',construct_notes=notes,alpha_variable=seqs[0],beta_variable=seqs[1],sequence_sources=ss,pmid=pmid,primary_locator=locator,sources=[src(exp/(pmc+'.xml'),f'https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/',locator)],mhc_allele=mhc,conditions=conditions,condition_limit='Assay-specific observation; missing conditions remain unknown; no absolute nonbinding or safety inference',training_label=None,use='reserved_external_only')
def add(base,peptide,kd,unit='uM',max_tested=None,endpoint='tcr_binding_spr',activation=None):
 row=dict(base);row.update(evidence_id=f"PMID{base['pmid']}:{base['receptor_name']}:{endpoint}:{peptide}",peptide=peptide,endpoint=endpoint,observed_outcome=('activation_detected' if activation else 'activation_not_detected') if activation is not None else ('measured_binding' if kd is not None else 'not_detected'),raw_value=str(kd) if kd is not None else ('NB' if activation is None else ('detected' if activation else 'not_detected')),raw_unit=unit if kd is not None else None,qualifier='=' if kd is not None else None,kd_nm=kd*({'uM':1000,'M':1e9}[unit]) if kd is not None else None,max_tested_um=max_tested)
 rows.append(row);return row
j=receptor('JDI','7PB2',[4,5],'36088370','PMC9464187','HLA-A*11:01','JDI','Biacore T200/8K; weak binding steady-state multicycle, 400-500 RU; JDI target concentration series 0.1-50 uM. Exact per-mimetic maximum and baseline temperature not stated. Pools use different loading and are excluded.','Parent JDI; not interchangeable with affinity-enhanced JDI variants','Tables 1, 3, 4; Methods SPR; Supplementary Figure S1a')
add(j,'VVVGADGVGK',63,max_tested=50)
for pep in ['VVVGAGGVGK','VVVGACGVGK','VVVGASGVGK','VVVGGGGVGK','VVVGDGGVGK','VVFGAGGVGK','VVLGSGGVGK']:add(j,pep,None)
a=receptor('JDIa41b1','7OW6',[4,5],'36088370','PMC9464187','HLA-A*11:01','JDI','Biacore T200/8K; mutant strong single-cycle 100-200 RU, injections 2.5,5,10,20,40 nM; WT ten two-fold injections starting 5 uM. Baseline temperature unspecified.','Affinity-enhanced alpha and beta construct; WT binds measurably; NEVER inherits parent NB','Tables 2, 3, 4; Supplementary Figure S1d,e')
add(a,'VVVGADGVGK',7.43e-10,'M');add(a,'VVVGAGGVGK',3e-6,'M')
n=receptor('N17.1.2','8YIV',[4,5],'39576860','PMC11584006','HLA-A*01:01','N17.1.2','BIAcore T100; 25 C; immobilized pMHC about 1000 RU; 50 uL/min; free-biotin blank; 600 s dissociation; positive 0.25-16 uM; WT up to 200 uM.','Human alpha/beta receptor using TRDV1 in alpha chain (not a gamma/delta receptor)','Figure 1; Results; Methods SPR and structure determination')
add(n,'ILDTAGKEEY',1.2,max_tested=16);add(n,'ILDTAGREEY',3.4,max_tested=16);add(n,'ILDTAGQEEY',None,max_tested=200)
for name,pdb,ids,kd in [('TCR4','7L1D',[4,5],49),('TCR3','7RRG',[3,4],200)]:
 p=receptor(name,pdb,ids,'35484264','PMC9117146','HLA-A*03:01',name,'Biacore T200; 25 C; TCR immobilized 900-1300 RU; pMHC analyte; duplicate injections; activated/deactivated reference; HEPES/NaCl/EDTA/P20 pH7.4. No WT SPR non-detection value inferred.','Human variable sequences; exact deposited recombinant construct; cellular experiments use murine constant domains','Figure 4b,c (TCR4, n=4); Extended Data Figure 3e (TCR3, n=3); Methods SPR')
 add(p,'ALHGGWTTK',kd)
 if name=='TCR4':
  cell=dict(p);cell.update(conditions='Human variable TCR4 with murine constant regions in CD8+ live mTCR+ T cells; HLA-A*03:01 targets; peptide 1 uM; TNF-alpha ICS; E:T 1:1; 6 h; 37 C; n=3.',primary_locator='Figure 3d,e; Methods intracellular cytokine staining',construct_notes='Cellular construct differs in constant domains from soluble SPR; separate endpoint')
  add(cell,'ALHGGWTTK',None,endpoint='tcr_cellular_tnf_alpha',activation=True);add(cell,'AHHGGWTTK',None,endpoint='tcr_cellular_tnf_alpha',activation=False)
report=validate_measurements(rows,r)
(d/'primary_measurements.json').write_text(json.dumps(rows,indent=2)+'\n')
# All selected receptors are reserved, including positive-only and cellular-only observations.
registry=json.loads((inputs/'holdout_registry.json').read_text());unique={}
for x in rows:
 unique[x['pmid'],x['receptor_name']]={k:x[k] for k in ['receptor_name','lineage','alpha_variable','beta_variable','pmid','sequence_sources']}
registry['receptors']=list(unique.values());registry['scope']='all selected constructs and corresponding near-chain families; JDI parent/enhanced share a biological lineage and cannot be independent folds';(d/'holdout_registry.json').write_text(json.dumps(registry,indent=2)+'\n')
queue=list(csv.DictReader((inputs/'primary_review_queue.tsv').open(),delimiter='\t'))
for x in queue:
 if x['receptor_name']=='CF34':status='blocked_primary_fulltext_conditions_unavailable';reason='IEDB retained; primary full-text/supplement conditions inaccessible. Multimer and SPR outcomes cannot be pooled; >200 uM is a censored bound, not exact KD or absolute negative.'
 elif x['endpoint']=='tcr_binding_spr_kd':status='provisional_not_in_strict_panel';reason='Original Figure 7 SPR and later published table support both ligands binding; original primary bytes not captured, so no strict provenance-complete admission.'
 else:status='excluded_from_intrinsic_binding_panel';reason='Cellular or multimer endpoint; CD4/CD8, peptide versus whole-protein and multi-construct assay aggregates differ. CD4 nonresponse cannot negate CD8 response or imply intrinsic nonbinding.'
 x.update(adjudication=status,reason=reason,training_label='',review_scope='endpoint/context adjudication; no invented numeric remeasurement')
write_tsv(d/'priority_adjudications.tsv',queue)
report.update(binding_outcomes=dict(Counter(x['observed_outcome'] for x in rows if x['endpoint']=='tcr_binding_spr')),priority_adjudications=dict(Counter(x['adjudication'] for x in queue)),limitations=['Targeted case-selected panel, not representative benchmark','JDI parent and enhanced share lineage; do not count as independent families','Conditional NB not proof of absent binding outside measured conditions','Cellular nonresponse does not imply intrinsic nonbinding','No observed WT SPR result fabricated for PIK3CA','Missing full sequences for JDIa9bwt/JDIa96b35 and workbook Mutant IDs: excluded; pooled peptides not split into negatives'])
(d/'primary_panel_validation.json').write_text(json.dumps(report,indent=2)+'\n')
# Retain raw supplementary replicate values and file/cell provenance, separate from aggregate pairs.
rep=[]
for name,pep,cell,values in [('JDI','VVVGADGVGK','D6:E6',[6.313e-5,6.252e-5]),('JDIa41b1','VVVGADGVGK','D19:E19',[7.25824717153411e-10,7.612904041720185e-10]),('JDIa41b1','VVVGAGGVGK','D23:F23',[3e-6,2.99e-6,3.01e-6])]:
 rep.append(dict(receptor_name=name,peptide=pep,raw_kd_values_M=values,sheet='Figure S1, Table 1, 2 and 4',cells=cell,source=src(exp/'Source Data.xlsx','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9464187/supplementaryFiles',cell),count_as_independent_pairs=False))
(d/'supplementary_replicates.json').write_text(json.dumps(rep,indent=2)+'\n')
manifest={'status':'curated_external_measurements_only','validation':report,'artifacts':{n:sha256(d/n) for n in ['primary_measurements.json','holdout_registry.json','priority_adjudications.tsv','primary_panel_validation.json','supplementary_replicates.json']},'source_retrieval_manifest':{'path':str((exp/'retrieval_manifest.json').relative_to(r)),'sha256':sha256(exp/'retrieval_manifest.json')},'supplementary_archive':{'sha256':sha256(exp/'PMC9464187_supplementary.zip'),'url':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9464187/supplementaryFiles'},'training_allowed':False}
(d/'expanded_evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(report,indent=2))
