"""Read primary P-TEAM workbooks without modifying them; no model or label fitting."""
from pathlib import Path
import argparse,hashlib,json,math,re,xml.etree.ElementTree as ET

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require_unique(rows,key):
    result={}
    for row in rows:
        k=row[key]
        if k in result:raise ValueError('duplicate source join key: '+str(k))
        result[k]=row
    return result

def build(root):
    import openpyxl
    d=root/'docs/six_blockers_20261021';sources=d/'sources'
    parent=json.loads((d/'pteam_primary_sequence_index.json').read_text())
    receptors=[];observations=[];sources_used={}
    for fname in ['pteam_cmv.xlsx','pteam_rnf43.xlsx']:
        path=sources/fname;w=openpyxl.load_workbook(path,read_only=True,data_only=True)
        sources_used[path.relative_to(root).as_posix()]=digest(path)
        cmv='cmv' in fname
        sheet='Mean' if cmv else 'Normalized by PC';matrix=list(w[sheet].values)
        if cmv:
            peps=require_unique([dict(id=x[0],peptide=x[1],position=x[2]) for x in list(w['peptides'].values)[1:]],'id')
            raw=list(w['Replicates'].values)
        else:
            peps=require_unique([dict(id=x[0],peptide=x[3],position=x[1]) for x in matrix[1:]],'id')
            raw=list(w['Individual APL screening'].values)
            pc={x[0]:(i,x[1],x[2]) for i,x in enumerate(w['Positive control'].values,1) if x[0] and isinstance(x[1],(float,int))}
        raw_rows=require_unique([dict(id=x[0],excel_row=i,values=x) for i,x in enumerate(raw[1:],2)],'id')
        for old in parent['receptors']:
            if old['source_file']!=fname:continue
            rid='PTEAM:'+('cmv:' if cmv else 'rnf43:')+old['name']
            seq=list(w['TCR sequences' if cmv else 'Sequences'].values)[old['source_row']-1]
            trav=seq[7] if cmv else seq[5];trbv=seq[3] if cmv else seq[7]
            rec=dict(receptor_id=rid,receptor_name=old['name'],pmid=old['pmid'],source_file=fname,sequence_source_row=old['source_row'],alpha_variable=old['alpha'].get('amino_acids',''),beta_variable=old['beta'].get('amino_acids',''),cdr3_alpha=old['cdr3_alpha'],cdr3_beta=old['cdr3_beta'],trav_raw=trav,trbv_raw=trbv,DNA_pair_passes=old['DNA_pair_passes'],author_excluded=old['author_excluded'],training_allowed=False,model_scored=False)
            rec['sequence_reasons']=[]
            if not rec['DNA_pair_passes']:rec['sequence_reasons'].append('DNA_translation_or_frame_CDR3_conflict')
            if rec['author_excluded']:rec['sequence_reasons'].append('author_excluded')
            receptors.append(rec)
            col=[i for i,x in enumerate(matrix[0]) if x==old['matrix_column']]
            if len(col)!=1:raise ValueError('ambiguous matrix receptor')
            for rowno,row in enumerate(matrix[1:],2):
                if not isinstance(row[col[0]],(float,int)):continue
                pep=peps[row[0]]['peptide']
                if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]{9}',pep):raise ValueError('noncanonical peptide')
                value=float(row[col[0]]);rr=raw_rows[row[0]];reps=[];controls=None
                rawcols=[i for i,x in enumerate(raw[0]) if x==old['matrix_column']]
                if cmv:
                    if len(rawcols)!=2:raise ValueError('expected paired replicate columns')
                    reps=[float(rr['values'][i]) for i in rawcols]
                    expected=sum(reps)/len(reps)
                else:
                    if len(rawcols)!=1:raise ValueError('ambiguous unnormalized source')
                    pr,mean,factor=pc[old['matrix_column']];unnorm=float(rr['values'][rawcols[0]])
                    expected=unnorm*factor;controls=dict(sheet='Positive control',source_row=pr,mean_PMA_ionomycin=mean,normalization_factor=factor,unnormalized_activation=unnorm)
                if not math.isfinite(value) or not math.isclose(expected,value,rel_tol=1e-10,abs_tol=1e-10):raise ValueError('primary response reconstruction mismatch')
                observations.append(dict(evidence_id=rid+':peptide:'+str(row[0]),receptor_id=rid,pmid='39151427',source_file=fname,source_sheet=sheet,source_excel_row=rowno,source_excel_column=col[0]+1,peptide_id=row[0],peptide=pep,is_cognate_parent=(row[0]==1),mhc_allele='HLA-A*02:01' if cmv else 'HLA-B*07:02',endpoint='JTPR_NFAT_flow_cytometry_activation',normalization='TCR_transgenic_fraction' if cmv else 'PMA_ionomycin_positive_control',value=value,normalized_replicates=reps,positive_control=controls,incubation_hours=24,peptide_pulse_concentration=None,pulse_dose_primary_verification='not_established_in_examined_methods',molecular_binding_label=None,clinical_label=None,binary_activation_label=None,training_allowed=False,model_scored=False,sequence_review_pass=not rec['sequence_reasons'],cognate_preselection='none_by_MHC_prediction' if cmv else 'NetMHCpan_4.1_prefiltered_APLs'))
    tree=ET.parse(sources/'PMC11480844.xml');authors=[]
    for node in tree.findall('.//article-meta//contrib'):
        name=node.find('name')
        if node.get('contrib-type') in (None,'author') and name is not None:
            authors.append(dict(surname=name.findtext('surname'),given_names=name.findtext('given-names')))
    sources_used[(sources/'PMC11480844.xml').relative_to(root).as_posix()]=digest(sources/'PMC11480844.xml')
    sources_used[(d/'pteam_primary_sequence_index.json').relative_to(root).as_posix()]=digest(d/'pteam_primary_sequence_index.json')
    return dict(purpose='primary endpoint intake, not a training release or frozen acceptance set',primary_pmid='39151427',primary_doi='10.1016/j.xgen.2024.100634',source_sha256=sources_used,authors=authors,receptors=receptors,observations=observations,normalizations_must_not_be_pooled=True,clinical_validation=False)

def reconcile_secondary(root,primary):
    import openpyxl
    path=root/'docs/six_blockers_20261021/sources/batcave_mhci.xlsx'
    def name(x):return re.sub(r'[^A-Z0-9]','',x.upper().removeprefix('TCR'))
    recs={x['receptor_id']:x for x in primary['receptors']}
    by_key=require_unique([dict(key=(name(recs[x['receptor_id']]['receptor_name']),x['peptide']),observation=x) for x in primary['observations']],'key')
    rows=[]
    for rowno,x in enumerate(openpyxl.load_workbook(path,read_only=True,data_only=True).active.values,1):
        if str(x[14])!='39151427':continue
        key=(name(x[0]),x[16]);target=by_key.get(key)
        if target is None:raise ValueError('secondary row unmatched to primary matrix')
        obs=target['observation'];raw=obs['positive_control']['unnormalized_activation'] if obs['positive_control'] else None
        eq=lambda v: v is not None and math.isclose(float(x[17]),v,abs_tol=1e-8,rel_tol=1e-8)
        rows.append(dict(secondary_excel_row=rowno,primary_evidence_id=obs['evidence_id'],secondary_value=x[17],matches_primary_normalized=eq(obs['value']),matches_primary_unnormalized=eq(raw),secondary_assay=x[10],primary_assay=obs['endpoint'],assay_name_conflict=x[10]=='NFAT luminescence',duplicate_primary_observation=True,admit_as_new_measurement=False))
    if len(rows)!=len(primary['observations']) or len({x['primary_evidence_id'] for x in rows})!=len(rows):raise ValueError('incomplete secondary bijection')
    return dict(source_sha256=digest(path),primary_pmid='39151427',rows=rows,new_independent_measurements=0,warning='RNA/functional response source scales differ; never concatenate secondary rows as new observations or pool normalization.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);p.add_argument('--reconcile-primary',type=Path);a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    root=Path(__file__).resolve().parents[1];result=reconcile_secondary(root,json.loads(a.reconcile_primary.read_text())) if a.reconcile_primary else build(root);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(rows=len(result['rows']),normalized_matches=sum(x['matches_primary_normalized'] for x in result['rows']),unnormalized_matches=sum(x['matches_primary_unnormalized'] for x in result['rows'])) if a.reconcile_primary else dict(receptors=len(result['receptors']),observations=len(result['observations']),sequence_pass_observations=sum(x['sequence_review_pass'] for x in result['observations']))))
if __name__=='__main__':main()
