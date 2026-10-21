"""Curate a source-backed, unscored human SPR detectability challenge."""
from pathlib import Path
import hashlib,json,openpyxl,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1];RAW=ROOT/'data/external_evidence/20261021-six-blockers';OUT=ROOT/'docs/six_blockers_20261021'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 output=OUT/'tcr4_primary_challenge.json'
 if output.exists():raise FileExistsError(output)
 refs={}
 for rec in (RAW/'rcsb_7L1D.fasta').read_text().split('>')[1:]:
  header,*lines=rec.splitlines()
  for chain in ['alpha','beta']:
   if ('receptor, '+chain+' chain') in header:refs[chain]=''.join(lines)
 if set(refs)!={'alpha','beta'} or 'CAGNTGTASKLTF' not in refs['alpha'] or 'CASSGLAGGPVSGANVLTF' not in refs['beta']:raise ValueError('Deposited receptor identity mismatch')
 workbook=openpyxl.load_workbook(RAW/'pi3kca_source_data.xlsx',read_only=True,data_only=True)
 data=list(workbook['Fig 1c'].values);curves={}
 for name,col in [('WT',0),('WT_repeat',2),('neo',5),('neo_repeat',7),('irrelevant_gp100',12)]:
  curves[name]=[dict(concentration_uM=float(row[col]),response_RU=float(row[col+1])) for row in data[2:] if isinstance(row[col],(float,int)) and isinstance(row[col+1],(float,int))]
 if max(x['concentration_uM'] for x in curves['WT'])!=208.24:raise ValueError('Original response range changed')
 text=ET.fromstring((RAW/'PMC11756396.xml').read_bytes())
 relevant_figures=[''.join(x.itertext()) for x in text.iter('fig') if x.get('id')=='Fig1']
 if not relevant_figures:relevant_figures=[''.join(x.itertext()) for x in text.iter('fig') if 'TCR4' in ''.join(x.itertext()) and '62' in ''.join(x.itertext())]
 if not any('4' in x and '°C' in x and 'six' in x for x in relevant_figures):raise ValueError('Figure-specific assay temperature/replicates not verified')
 construct='PMID:39833157:TCR4:PDB:7L1D'
 common=dict(construct_id=construct,receptor_name='TCR4',alpha_pdb_construct=refs['alpha'],beta_pdb_construct=refs['beta'],cdr3_alpha='CAGNTGTASKLTF',cdr3_beta='CASSGLAGGPVSGANVLTF',species='human',method='SPR',instrument='Biacore T200',temperature_K=277.15,temperature_basis='published Figure 1 caption, panels c/d, 4 C',orientation='TCR immobilized; soluble peptide-HLA analyte',chip='CM5 Series S',immobilization_RU_general_method_range=[400,6000],flow_uL_min=5,buffer='HBS-EP, pH 7.4',fit='1:1 steady-state; OriginPro 2024',clinical_risk_label=None,full_clinical_engineered_DNA_verified=False,model_scored=False,training_allowed=False,formal_independent_acceptance_allowed=False)
 observations=[dict(common,evidence_id=construct+':ALHGGWTTK:HLA-A0301',peptide='ALHGGWTTK',mhc_allele='HLA-A*03:01',endpoint='SPR_affinity',assay_result='binding_detected',reported_KD_uM=62,reported_KD_sd_uM=6,n_published_affinity_replicates=6,point_regression_target_available=True,assay_binding_detected=True,universal_binding_label=None),dict(common,evidence_id=construct+':AHHGGWTTK:HLA-A0301',peptide='AHHGGWTTK',mhc_allele='HLA-A*03:01',endpoint='SPR_detectability',assay_result='no_detectable_binding_under_reported_conditions',reported_KD_uM=None,reported_KD_sd_uM=None,assayed_analyte_max_uM=208.24,quantified_detection_limit_RU=None,valid_KD_lower_bound_uM=None,point_regression_target_available=False,assay_binding_detected=False,universal_binding_label=None)]
 result=dict(status='unscored_primary_assay_challenge_not_training_release',doi='10.1038/s41467-025-56004-8',pmid='39833157',pdb_identity='7L1D TCR4 per main text reference and original RCSB chain definitions; not any other paper named TCR4',observations=observations,measured_source_curves=curves,positive_sample_integrity_control=dict(reagent='s3-4 scTv',WT_KD_uM=1.8,WT_KD_sd_uM=.1,neo_KD_uM=1.7,neo_KD_sd_uM=.2,n_published_affinity_replicates=3),irrelevant_MHC_control_included_in_model_dataset=False,limits=['No numeric KD invented for undetectable binding; maximum analyte concentration is not a KD lower bound without Rmax/LOD assumptions','Full PDB construct proteins are available; no inferred clinical engineering DNA or patient outcome','One receptor lineage and one research source; not a sufficient independent acceptance cohort','Figure 1b source spreadsheet says ug/mL while published axis says uM; cellular numeric dose-response not admitted or converted here'],sources={name:sha(RAW/name) for name in ['PMC11756396.xml','rcsb_7L1D.fasta','rcsb_7L1D.json','pi3kca_source_data.xlsx','pi3kca_figure1.png']})
 output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print('PRIMARY_CHALLENGE',len(observations))
if __name__=='__main__':main()
