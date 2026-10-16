"""Assay-level experimental evidence; database calls never imply a universal label.

Offline processing of provenance-tracked IEDB API snapshots plus explicit
primary-paper curation. No random negatives, cross-assay majority voting or
single-chain epitope transfer. The frozen pilot is an external sanity panel,
not a population benchmark or training dataset.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re
from collections import Counter,defaultdict
from pathlib import Path
from .build_release import domain,similar
from .release import sha256

AA=set('ACDEFGHIKLMNPQRSTVWY')

def assay_ids(value):
    if value is None or value=='':return []
    if not isinstance(value,str) or not re.fullmatch(r'[0-9,;\s]+',value):raise ValueError('malformed assay ID list')
    return sorted({int(x) for x in re.split(r'[,;\s]+',value.strip()) if x})


def endpoint(method,response):
    # Retain assay modality and measured response; do not pool endpoint families.
    m=(method or '').casefold();r=(response or '').casefold()
    if m=='surface plasmon resonance (spr)' and r=='dissociation constant kd':return 'tcr_binding_spr_kd'
    if m=='multimer/tetramer' and r=='qualitative binding':return 'tcr_multimer_binding'
    if 'release' in r:return 'tcr_cellular_'+re.sub('[^a-z0-9]+','_',r).strip('_')
    if r in {'cytotoxicity','activation','degranulation'}:return 'tcr_cellular_'+r
    return 'other_assay'


def observed_call(value):
    if value=='Negative':return 'conditional_negative'
    if (value or '').startswith('Positive'):return 'positive'
    return 'unknown'


def finite_nm(value,unit):
    factors={'nM':1.,'uM':1000.,'µM':1000.,'μM':1000.,'mM':1e6,'M':1e9,'pM':.001}
    if unit not in factors or value is None:return None
    try:v=float(value)*factors[unit]
    except (TypeError,ValueError):return None
    return v if math.isfinite(v) and v>0 else None


def paired_sequences(row):
    if row.get('receptor__type')!='alphabeta':raise ValueError('not alpha/beta')
    result=[]
    for n,kind in [(1,'alpha'),(2,'beta')]:
        if row.get(f'chain_{n}__type')!=kind:raise ValueError('chain identity mismatch')
        if not (row.get(f'chain_{n}__organism_iri') or '').endswith('9606'):raise ValueError('not confirmed human chain')
        seq=row.get(f'chain_{n}__protein_sequence') or row.get(f'chain_{n}__v_domain_calculated') or ''
        if len(seq)<70 or not set(seq)<=AA:raise ValueError('complete canonical paired variable sequence unavailable')
        cdr=row.get(f'chain_{n}__cdr3_calculated') or row.get(f'chain_{n}__cdr3_curated') or ''
        if not cdr or cdr not in seq:raise ValueError('CDR3/variable-sequence mismatch')
        result.append(domain(seq)[0])
    return result


def pair_key(alpha,beta):return hashlib.sha256((alpha+'|'+beta).encode()).hexdigest()


def overlaps(alpha,beta,references):
    # Conservative reservation: either homologous chain excludes independent use.
    return sorted({x['family_id'] for x in references if similar(alpha,x['alpha']) or similar(beta,x['beta'])})


def make_ledger(exports,assays,search,references):
    by_assay={x['assay_id']:x for x in assays};by_search={x['tcell_id']:x for x in search};out=[];seen={};cache={}
    if len(by_assay)!=len(assays) or len(by_search)!=len(search):raise ValueError('duplicate assay source IDs')
    for receptor in exports:
        if receptor.get('assay__type')!='T cell':continue
        for aid in assay_ids(receptor.get('assay__iedb_ids')):
            a=by_assay[aid];s=by_search[aid];rid=receptor['receptor__iedb_receptor_id'];gid=receptor['receptor_group_id']
            reasons=[]
            if rid not in (s.get('receptor_ids') or []) or gid not in (s.get('tcr_receptor_group_ids') or []):reasons.append('assay_receptor_link_mismatch')
            if not receptor.get('reference__iedb_iri') or receptor['reference__iedb_iri'].rstrip('/').split('/')[-1] != str(s.get('reference_id')):reasons.append('reference_join_mismatch')
            if receptor.get('epitope__iedb_iri','').rstrip('/').split('/')[-1]!=str(s.get('structure_id')):reasons.append('epitope_join_mismatch')
            pep=a.get('epitope__name') or ''
            if pep!=s.get('linear_sequence') or pep!=receptor.get('epitope__name') or not set(pep)<=AA or not 5<=len(pep)<=20:reasons.append('noncanonical_or_mismatched_peptide')
            if a.get('epitope__modified_residues') or a.get('epitope__modifications'):reasons.append('modified_peptide')
            mhc=a.get('mhc_restriction__name') or ''
            if not re.fullmatch(r'HLA-[ABC]\*\d{2,3}:\d{2,3}',mhc) or mhc!=s.get('mhc_allele_name'):reasons.append('class_I_allele_not_unambiguous')
            try:alpha,beta=paired_sequences(receptor)
            except ValueError as exc:alpha=beta='';reasons.append(str(exc))
            fam=[]
            if alpha:
                k=(alpha,beta)
                if k not in cache:cache[k]=overlaps(alpha,beta,references)
                fam=cache[k]
                if fam:reasons.append('overlap_with_structural_release')
            method=a.get('assay__method') or '';response=a.get('assay__response_measured') or '';ep=endpoint(method,response)
            if ep=='other_assay':reasons.append('not_target_endpoint')
            outcome=observed_call(a.get('assay__qualitative_measurement'))
            if outcome=='unknown':reasons.append('unknown_observation')
            raw=a.get('assay__quantitative_measurement');unit=a.get('assay__units') or '';qual=a.get('assay__measurement_inequality') or ''
            if qual not in {'','='}:reasons.append('censored_measurement')
            source_variant=hashlib.sha256(json.dumps(receptor,sort_keys=True).encode()).hexdigest()
            key=(aid,rid,source_variant)
            row={'evidence_id':f'IEDB:{aid}:receptor:{rid}:export:{source_variant[:16]}','assay_id':aid,'receptor_id':rid,'receptor_group_id':gid,'receptor_name':receptor.get('receptor__reference_name') or '',
                 'alpha_variable':alpha,'beta_variable':beta,'peptide':pep,'mhc_allele':mhc,'endpoint':ep,'method':method,'response':response,'observed_outcome':outcome,
                 'raw_value':raw if raw is not None else '', 'raw_unit':unit,'qualifier':qual,'value_nm':finite_nm(raw,unit) if ep=='tcr_binding_spr_kd' else '',
                 'pmid':a.get('reference__pmid') or '', 'paper_location':a.get('assay__location_of_assay_data_in_reference') or '',
                 'assay_comments':a.get('assay__comments') or '', 'immunization_comments':a.get('immunization__comments') or '',
                 'effector_cell':a.get('effector_cell__name') or '', 'effector_culture':a.get('effector_cell__culture_condition') or '',
                 'apc':a.get('antigen_presenting_cell__name') or '', 'apc_culture':a.get('antigen_presenting_cell__culture_condition') or '',
                 'antigen_tested':a.get('assay_antigen__name') or '', 'antigen_relation':a.get('assay_antigen__epitope_relation') or '',
                 'assay_url':f'https://www.iedb.org/assay/{aid}','reference_url':f'https://pubmed.ncbi.nlm.nih.gov/{a.get("reference__pmid")}/',
                 'structural_overlap_families':';'.join(fam),'screening_reasons':';'.join(sorted(set(reasons))), 'evidence_status':'database_curated_primary_review_pending','label':'','use':'review_only'}
            if key in seen:
                if seen[key]!=row:raise ValueError('conflicting duplicate assay/receptor export')
                continue
            seen[key]=row;out.append(row)
    duplicates=Counter((x['assay_id'],x['receptor_id']) for x in out)
    for x in out:
        if duplicates[x['assay_id'],x['receptor_id']]>1:x['screening_reasons']=';'.join(sorted(set(filter(None,x['screening_reasons'].split(';')))|{'ambiguous_receptor_export'}))
    return sorted(out,key=lambda x:(x['assay_id'],x['receptor_id'],x['evidence_id']))


def condition_conflicts(rows):
    groups=defaultdict(list)
    for row in rows:
        key=(row['receptor_id'],row['peptide'],row['mhc_allele'],row['endpoint'],row['method'])
        groups[key].append(row)
    return [{'receptor_id':k[0],'peptide':k[1],'mhc_allele':k[2],'endpoint':k[3],'method':k[4],
             'assay_ids':sorted(x['assay_id'] for x in v),'resolution':'retain_conditions_and_primary_review; no positive-over-negative voting'}
            for k,v in sorted(groups.items()) if {'positive','conditional_negative'} <= {x['observed_outcome'] for x in v}]


def verify_pilot(rows,references):
    if len({x['evidence_id'] for x in rows})!=len(rows):raise ValueError('duplicate pilot evidence')
    groups=defaultdict(list)
    for x in rows:
        for key in ['primary_url','primary_locator','source_sha256','alpha_variable','beta_variable','conditions','negative_interpretation','receptor_sequence_source']:
            if not x.get(key):raise ValueError('pilot provenance missing '+key)
        if not re.fullmatch(r'HLA-[ABC]\*\d{2,3}:\d{2,3}',x.get('mhc_allele','')):raise ValueError('pilot MHC allele is ambiguous')
        if x['label'] not in [0,1] or x['endpoint']!='tcr_binding_spr_detection':raise ValueError('pilot endpoint/label mismatch')
        if not set(x['alpha_variable']+x['beta_variable']+x['peptide'])<=AA or min(len(x['alpha_variable']),len(x['beta_variable']))<70:raise ValueError('pilot sequence invalid')
        if x['label']==0 and (x['raw_value']!='not_detected' or x['kd_nm'] is not None):raise ValueError('non-detection is not a zero or infinite KD')
        if x['label']==1 and (x['kd_nm'] is None or not math.isfinite(x['kd_nm']) or x['kd_nm']<=0):raise ValueError('invalid measured KD')
        if overlaps(x['alpha_variable'],x['beta_variable'],references):raise ValueError('pilot/structural family overlap')
        groups[pair_key(x['alpha_variable'],x['beta_variable'])].append(x)
    for group in groups.values():
        if {x['label'] for x in group}!={0,1}:raise ValueError('pilot needs matched measured controls per receptor')
    return {'n_records':len(rows),'n_receptors':len(groups),'n_positive':sum(x['label'] for x in rows),'n_conditional_negative':sum(x['label']==0 for x in rows),
            'n_studies':len({x['pmid'] for x in rows}),'training_allowed':False,'benchmark_allowed':False,'use':'external_sanity_panel_only',
            'limitations':['single-study selected panel','wild-type/mutant peptide pair only','no broad specificity or clinical accuracy claim','negative only within stated SPR concentration and assay conditions']}


def write_tsv(path,rows):
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


def build(snapshot,reference_path,pilot_path,output):
    # Verify downloaded source fingerprints before using their derived all-row views.
    sources=json.loads((snapshot/'retrieval_manifest.json').read_text())
    for item in sources:
        p=(snapshot/item['file']).resolve()
        if not p.is_relative_to(snapshot.resolve()) or sha256(p)!=item['sha256']:raise ValueError('changed source snapshot: '+item['file'])
    references=json.loads(reference_path.read_text());arrays={}
    for name in ['tcr_export','tcell_export','tcell_search']:
        batches=sorted(snapshot.glob(name+'_batch_*.json'),key=lambda p:int(p.stem.rsplit('_',1)[1]))
        if not batches or any(p.name not in {x['file'] for x in sources} for p in batches):raise ValueError('untracked or missing API batch')
        arrays[name]=[x for p in batches for x in json.loads(p.read_text())]
    rows=make_ledger(arrays['tcr_export'],arrays['tcell_export'],arrays['tcell_search'],references)
    pilot=json.loads(pilot_path.read_text())
    known={x['file'] for x in sources}
    for row in pilot:
        if row['primary_file'] not in known:raise ValueError('untracked primary paper')
        if sha256(snapshot/row['primary_file'])!=row['source_sha256']:raise ValueError('changed primary paper')
        entities=json.loads(row['receptor_sequence_source'])
        if len(entities)!=2:raise ValueError('paired sequence sources required')
        for source,seq,cdr in zip(entities,[row['alpha_variable'],row['beta_variable']],[row['cdr3_alpha'],row['cdr3_beta']]):
            if source['file'] not in known:raise ValueError('untracked receptor source')
            p=snapshot/source['file']
            if sha256(p)!=source['sha256']:raise ValueError('changed receptor source')
            construct=json.loads(p.read_text())['entity_poly']['pdbx_seq_one_letter_code_can']
            if domain(construct)[0]!=seq or cdr not in seq:raise ValueError('curated receptor sequence does not match source')
    summary=verify_pilot(pilot,references);conflicts=condition_conflicts(rows)
    output.mkdir(parents=True,exist_ok=False)
    write_tsv(output/'assay_evidence.tsv',rows);write_tsv(output/'external_sanity_panel.tsv',pilot)
    (output/'condition_conflicts.json').write_text(json.dumps(conflicts,indent=2))
    registry={'purpose':'reserved_external_test_only','scope':'entire_receptor_family_and_future_predicted_copies','rule':'reject either corresponding alpha/beta variable-domain near match, or CDR3 core >=0.9 SequenceMatcher similarity with minimum 8 residues; conservative reservation, not specificity equivalence','receptors':[{k:x[k] for k in ['receptor_name','alpha_variable','beta_variable','pmid','receptor_sequence_source']} for x in pilot if x['label']==1]}
    (output/'holdout_registry.json').write_text(json.dumps(registry,indent=2))
    summary.update(n_database_evidence=len(rows),observed_outcomes=dict(Counter(x['observed_outcome'] for x in rows)),n_primary_review_pending=len(rows),n_condition_conflict_groups=len(conflicts),n_screened_candidates=sum(not x['screening_reasons'] for x in rows))
    (output/'summary.json').write_text(json.dumps(summary,indent=2))
    m={'schema_version':1,'status':'frozen_external_sanity_panel','training_allowed':False,'benchmark_allowed':False,'pilot_source_sha256':sha256(pilot_path),'structural_reference_sha256':sha256(reference_path),'retrieval_sources':sources,'code_files':{str(p.relative_to(Path(__file__).resolve().parents[3])):sha256(p) for p in [Path(__file__),Path(__file__).with_name('fetch_experimental.py'),Path(__file__).with_name('build_release.py'),Path(__file__).with_name('release.py')]},
       'artifacts':{p.name:sha256(p) for p in output.iterdir() if p.is_file()},'selection':'human paired alpha/beta receptor groups with at least one IEDB negative observation; targeted discovery, not representative sampling'}
    (output/'evidence_manifest.json').write_text(json.dumps(m,indent=2));print(json.dumps(summary));return summary


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--structural-records',type=Path,required=True);p.add_argument('--pilot',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();build(a.snapshot,a.structural_records,a.pilot,a.out)



def _cdr3_core(sequence):
    s=sequence or ''
    if s.startswith('C'):s=s[1:]
    if s.endswith(('F','W')):s=s[:-1]
    return s


from functools import lru_cache
from difflib import SequenceMatcher

@lru_cache(maxsize=300000)
def _near_cdr3(a,b):
    a,b=_cdr3_core(a),_cdr3_core(b)
    if min(len(a),len(b))<8 or 2*min(len(a),len(b))/(len(a)+len(b))<.9:return False
    return SequenceMatcher(None,a,b,autojunk=False).ratio()>=.9


def reserved_hits(record,registry):
    result=[]
    alpha=record.get('alpha') or record.get('alpha_variable') or ''
    beta=record.get('beta') or record.get('beta_variable') or ''
    ca=record.get('cdr3_alpha') or (domain(alpha)[1] if alpha else '')
    cb=record.get('cdr3_beta') or (domain(beta)[1] if beta else '')
    for r in registry:
        if ((alpha and similar(alpha,r['alpha_variable'])) or (beta and similar(beta,r['beta_variable']))
            or _near_cdr3(ca,domain(r['alpha_variable'])[1]) or _near_cdr3(cb,domain(r['beta_variable'])[1])):result.append(r['receptor_name'])
    return sorted(set(result))


def assert_reserved_absent(records,registry_path=None):
    if registry_path is None:registry_path=Path(__file__).resolve().parents[3]/'docs/evidence/holdout_registry.json'
    if not registry_path.exists():raise ValueError('external holdout reservation ledger is missing')
    registry=json.loads(registry_path.read_text())['receptors']
    checked=set()
    for row in records:
        key=tuple(row.get(x,'') for x in ['alpha','beta','alpha_variable','beta_variable','cdr3_alpha','cdr3_beta'])
        if key in checked:continue
        checked.add(key);hits=reserved_hits(row,registry)
        if hits:raise ValueError('reserved external receptor family cannot enter training release: '+','.join(hits))


def filter_associations(source,output,quarantine,registry_path):
    registry=json.loads(registry_path.read_text())['receptors'];kept=excluded=0
    if output.exists() or quarantine.exists():raise ValueError('refusing to overwrite association files')
    output.parent.mkdir(parents=True,exist_ok=True);quarantine.parent.mkdir(parents=True,exist_ok=True)
    tmp=output.with_suffix(output.suffix+'.tmp');qtmp=quarantine.with_suffix(quarantine.suffix+'.tmp')
    with source.open(newline='') as src,tmp.open('x',newline='') as out,qtmp.open('x',newline='') as reject:
        reader=csv.DictReader(src,delimiter='\t');fields=reader.fieldnames
        if not fields or not {'cdr3_alpha','cdr3_beta','label'}<=set(fields):raise ValueError('association schema missing')
        writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t');writer.writeheader()
        qw=csv.DictWriter(reject,fieldnames=fields+['source_line','reservation_reason'],delimiter='\t');qw.writeheader()
        for line,row in enumerate(reader,2):
            if row['label']:raise ValueError('association-only input must be unlabeled')
            hits=reserved_hits(row,registry)
            if hits:qw.writerow(dict(row,source_line=line,reservation_reason='reserved_external:'+','.join(hits)));excluded+=1
            else:writer.writerow(row);kept+=1
    tmp.replace(output);qtmp.replace(quarantine)
    return {'source_sha256':sha256(source),'output_sha256':sha256(output),'quarantine_sha256':sha256(quarantine),'registry_sha256':sha256(registry_path),'n_input':kept+excluded,'n_retained':kept,'n_quarantined':excluded,'rule':'either corresponding alpha/beta CDR3 core >=0.9 SequenceMatcher similarity, minimum 8 residues; conservative reservation, not specificity equivalence','training_labels_added':False}


if __name__=='__main__':main()
