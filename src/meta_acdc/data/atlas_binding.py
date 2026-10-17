"""ATLAS experimental affinity ledger; no synthetic or fabricated binary labels.

Only an exact human class-I complex, paired source sequences, a measured SPR
KD and unmutated assay MHC can enter the quantitative exploratory dataset.
Database curation is provenance, not a claim of independent primary-paper review.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re
from collections import Counter,defaultdict
from pathlib import Path
from difflib import SequenceMatcher
from .release import sha256
from .build_release import domain,similar
from .experimental_evidence import AA,reserved_hits

NUM=r'(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'

def kd_measurement(raw):
    m=re.fullmatch(r'\s*([<>~]?)\s*('+NUM+r')(?:\s*(?:\+/-|±)\s*('+NUM+r'))?\s*',str(raw))
    if not m:return {'qualifier':None,'kd_um':None,'error_um':None,'pkd':None,'status':'missing_or_unparsed'}
    qualifier=m[1] or '=';value=float(m[2]);error=float(m[3]) if m[3] else None
    if not math.isfinite(value) or value<=0:return {'qualifier':qualifier,'kd_um':None,'error_um':error,'pkd':None,'status':'invalid'}
    return {'qualifier':qualifier,'kd_um':value,'error_um':error,'pkd':6-math.log10(value) if qualifier=='=' else None,'status':'exact' if qualifier=='=' else 'censored_or_approximate'}


def fasta_records(text):
    rows=[]
    for block in text.split('>')[1:]:
        lines=block.splitlines();parts=lines[0].split('|')
        if len(parts)<4:raise ValueError('FASTA chain/source header missing')
        chains=[]
        for token in re.sub(r'^Chains?\s+','',parts[1]).split(','):
            auth=re.search(r'\[auth ([^]]+)\]',token)
            chains.append(auth[1] if auth else token.strip())
        seq=''.join(lines[1:]).strip()
        rows.append({'header':lines[0],'chains':chains,'sequence':seq,'human':'(9606)' in '|'.join(parts[3:])})
    if not rows:raise ValueError('empty FASTA')
    return rows


def paired_from_sources(fasta,annotations):
    pair=[];cdr=[]
    for kind in ['Alpha','Beta']:
        ann=[x for x in annotations if x['tcr_type']==kind]
        choices={}
        for x in ann:
            cores=set(x['cdr3_sequences'].split(','));chains=set(x['cdr3_chains'].split(','))
            for row in fasta:
                if chains.intersection(row['chains']) and row['human'] and len(row['sequence'])>=70 and set(row['sequence'])<=AA:
                    found=[c for c in cores if c and c in row['sequence']]
                    if len(found)==1:choices[row['sequence']]=found[0]
        if len(choices)!=1:raise ValueError('missing/ambiguous human '+kind+' construct')
        seq,core=next(iter(choices.items()));pair.append(domain(seq)[0]);cdr.append(core)
    return pair,cdr


def grouped_folds(rows,*,unseen_peptide=True):
    parent=list(range(len(rows)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    def union(i,j):parent[find(i)]=find(j)
    for i,a in enumerate(rows):
        for j,b in enumerate(rows[:i]):
            if (a['pmid']==b['pmid'] or a['pair_id']==b['pair_id'] or similar(a['alpha_variable'],b['alpha_variable']) or similar(a['beta_variable'],b['beta_variable']) or (unseen_peptide and SequenceMatcher(None,a['peptide'],b['peptide'],autojunk=False).ratio()>=.8)):union(i,j)
    groups=defaultdict(list)
    for i in range(len(rows)):groups[find(i)].append(i)
    blocks=sorted(groups.values(),key=lambda indices:(-len(indices),min(rows[i]['evidence_id'] for i in indices)))
    if len(blocks)<3:raise ValueError('fewer than three independent connected blocks')
    folds={x:[] for x in ['train','validation','test']};targets={'train':.6*len(rows),'validation':.2*len(rows),'test':.2*len(rows)}
    for indices in blocks:
        fold=min(folds,key=lambda f:len(folds[f])/targets[f]);folds[fold].extend(indices)
        key='B'+hashlib.sha256('|'.join(sorted(rows[i]['evidence_id'] for i in indices)).encode()).hexdigest()[:16]
        for i in indices:rows[i]['split_block']=key;rows[i]['fold']=fold
    for fold in folds:folds[fold].sort()
    return folds,{'n_blocks':len(blocks),'block_sizes':list(map(len,blocks)),'policy':'same study, exact receptor-ligand pair, near corresponding alpha OR beta variable chain; additionally peptide similarity >=0.8 when enabled','unseen_peptide':unseen_peptide}


def build(snapshot,registry_path,output):
    sources=json.loads((snapshot/'retrieval_manifest.json').read_text());tracked={x['file']:x for x in sources}
    for source in sources:
        path=(snapshot/source['file']).resolve()
        if not path.is_relative_to(snapshot.resolve()) or sha256(path)!=source['sha256']:raise ValueError('changed source '+source['file'])
    atlas=list(csv.DictReader((snapshot/'atlas_Mutants.tsv').open(),delimiter='\t'));complexes={x['PDB_ID']:x for x in csv.DictReader((snapshot/'tcr3d_complexes.tsv').open(),delimiter='\t')};chains=defaultdict(list)
    for x in csv.DictReader((snapshot/'tcr3d_chain.tsv').open(),delimiter='\t'):chains[x['pdb_id']].append(x)
    clinical_path=Path(__file__).resolve().parents[3]/'docs/clinical_primary_20261018/clinical_reservations.json'
    if not clinical_path.exists():raise ValueError('clinical holdout reservations missing')
    registry=json.loads(registry_path.read_text())['receptors']+json.loads(clinical_path.read_text())['receptors'];ledger=[];candidate=[]
    for raw in atlas:
        row={'evidence_id':'ATLAS:'+raw['Ind'],'raw_record':raw,'training_label':None,'measurement':kd_measurement(raw['Kd_microM']),'reasons':[],'primary_review':'database_curated; not independently reviewed per paper'};why=row['reasons'];pdb=raw['true_PDB'].lower();meta=complexes.get(pdb)
        if not meta or meta['TCR_complex']!='CLASSI' or meta['TCR_organism']!='Human':why.append('not_exact_human_class_I_complex')
        if raw['MHC_mut']!='WT':why.append('mutant_MHC_not_supported')
        if not re.fullmatch(r'HLA-[ABC]\*\d{2,3}:\d{2,3}',raw['MHCname']):why.append('allele_not_explicit_classical_HLA')
        if row['measurement']['status']!='exact':why.append(row['measurement']['status'])
        if not (raw['Exp_Method'].startswith('Biacore') or raw['Exp_Method']=='SPR'):why.append('not_SPR_endpoint')
        if not raw['PMID'].isdigit():why.append('missing_primary_reference')
        pep=raw['PEPseq']
        if not 5<=len(pep)<=20 or not set(pep)<=AA:why.append('peptide_not_canonical')
        file='rcsb_'+pdb+'.fasta'
        if not why:
            try:
                if file not in tracked:raise ValueError('missing tracked exact sequence source')
                f=fasta_records((snapshot/file).read_text());pair,cdr=paired_from_sources(f,chains[pdb])
                peptides={x['sequence'] for x in f if 5<=len(x['sequence'])<=20}
                if peptides!={pep}:raise ValueError('peptide sequence or complex ambiguity')
                identity={'alpha_variable':pair[0],'beta_variable':pair[1],'cdr3_alpha':cdr[0],'cdr3_beta':cdr[1]}
                hits=reserved_hits(identity,registry)
                if hits:raise ValueError('reserved_external_family:'+','.join(hits))
                row.update(identity);row.update(receptor_name=raw['TCRname'],peptide=pep,mhc_allele=raw['MHCname'],pmid=raw['PMID'],pdb=pdb,method=raw['Exp_Method'],temperature_K=raw['Temperature_K'],source_file=file,source_sha256=tracked[file]['sha256'],source_url=tracked[file]['url'],pair_id=hashlib.sha256('|'.join(pair+[pep,raw['MHCname']]).encode()).hexdigest(),pkd=row['measurement']['pkd'])
                candidate.append(row)
            except ValueError as exc:why.append(str(exc))
        ledger.append(row)
    # Same construct/ligand/study/conditions with discordant KDs is quarantined, not averaged.
    duplicates=defaultdict(list)
    for x in candidate:duplicates[(x['pair_id'],x['pmid'],x['method'],x['temperature_K'],x['raw_record']['SPR_SensorChip'],x['raw_record']['Coupling_Method'],x['raw_record']['Analyte'])].append(x)
    eligible=[]
    for values in duplicates.values():
        if len({x['pkd'] for x in values})>1:
            for x in values:x['reasons'].append('discordant_same_pair_study_conditions')
        else:
            eligible.append(values[0])
            for x in values[1:]:x['reasons'].append('duplicate_measurement_not_independent')
    eligible.sort(key=lambda x:x['evidence_id']);folds,grouping=grouped_folds(eligible)
    # Never recycle previously exposed validation/test points into training after filtering.
    parent_path=Path(__file__).resolve().parents[3]/'docs/database_binding/affinity_dataset.json'
    if parent_path.exists():
        previous={x['evidence_id']:x for x in json.loads(parent_path.read_text())}
        if any(x['evidence_id'] not in previous for x in eligible):raise ValueError('new assay rows require an explicit new split protocol')
        retained_ids={x['evidence_id'] for x in eligible}
        blocked={x['split_block'] for key,x in previous.items() if key not in retained_ids}
        for x in eligible:
            if previous[x['evidence_id']]['split_block'] in blocked:x['reasons'].append('reserved_parent_connected_block')
        eligible=[x for x in eligible if previous[x['evidence_id']]['split_block'] not in blocked]
        for x in eligible:
            old=previous[x['evidence_id']]
            if old['pair_id']!=x['pair_id'] or old['pkd']!=x['pkd']:raise ValueError('parent assay changed')
            x['fold']=old['fold'];x['split_block']=old['split_block']
        folds={f:[i for i,x in enumerate(eligible) if x['fold']==f] for f in ['train','validation','test']}
        grouping['parent_fold_inheritance']='old fold and connected block retained; exclusions never promote validation/test to train'
        grouping['parent_dataset_sha256']=sha256(parent_path)
        grouping['n_blocks']=len({x['split_block'] for x in eligible})
        grouping['block_sizes']=list(Counter(x['split_block'] for x in eligible).values())

    output.mkdir(parents=True,exist_ok=False)
    def save(name,value):(output/name).write_text(json.dumps(value,indent=2)+'\n')
    save('measurement_ledger.json',ledger);save('affinity_dataset.json',eligible);save('splits.json',folds)
    report={'n_source_records':len(atlas),'n_eligible_exact_SPR':len(eligible),'n_studies':len({x['pmid'] for x in eligible}),'n_receptor_names':len({x['receptor_name'] for x in eligible}),'fold_counts':{k:len(v) for k,v in folds.items()},'fold_blocks':{k:len({eligible[i]['split_block'] for i in v}) for k,v in folds.items()},'grouping':grouping,'exclusions':dict(Counter(reason for x in ledger for reason in x['reasons'])),'task':'experimental_SPR_pKD_regression','binary_binding_training_allowed':False,'clinical_training_allowed':False,'scope':'exploratory database-curated measured-affinity research only; censored and unobserved pairs are not negatives; no broad affinity or binary-binding accuracy accepted'}
    save('summary.json',report)
    manifest={'schema_version':1,'task':report['task'],'status':'frozen_exploratory_affinity_dataset','binary_training_allowed':False,'clinical_training_allowed':False,'parent_dataset_sha256':sha256(parent_path) if parent_path.exists() else None,'clinical_registry_sha256':sha256(clinical_path),'external_registry_sha256':sha256(registry_path),'sources':sources,'source_snapshot':str(snapshot.resolve()),'registry_path':str(registry_path.resolve()),'code_files':{name:sha256(Path(__file__).resolve().parents[3]/name) for name in ['src/meta_acdc/data/atlas_binding.py','src/meta_acdc/data/experimental_evidence.py','src/meta_acdc/data/build_release.py','src/meta_acdc/data/release.py','src/meta_acdc/models/affinity_baseline.py']},'files':{p.name:sha256(p) for p in output.iterdir() if p.is_file()},'summary':report};save('manifest.json',manifest);return report


def verify_release(release):
    manifest=json.loads((release/'manifest.json').read_text());root=Path(__file__).resolve().parents[3]
    primary_path=root/'docs/spr_primary_20261018/review_ledger.json'
    additions_path=root/'docs/spr_primary_20261018/clinical_sequence_reservations.json'
    if not primary_path.exists() or manifest.get('primary_review_sha256')!=sha256(primary_path):raise ValueError('primary measurement audit absent or changed; database candidates cannot train')
    if not additions_path.exists() or manifest.get('clinical_additions_sha256')!=sha256(additions_path):raise ValueError('clinical sequence additions absent or changed')
    clinical_path=root/'docs/clinical_primary_20261018/clinical_reservations.json'
    parent_path=Path(manifest.get('parent_dataset_path',str(root/'docs/database_binding/affinity_dataset.json')))
    if manifest.get('parent_dataset_sha256') is not None and (not parent_path.exists() or sha256(parent_path)!=manifest['parent_dataset_sha256']):raise ValueError('parent split history changed')
    if not clinical_path.exists() or manifest.get('clinical_registry_sha256')!=sha256(clinical_path):raise ValueError('clinical reservation absent or changed; rebuild a new release')
    if manifest['task']!='experimental_SPR_pKD_regression':raise ValueError('not quantitative affinity data')
    if not manifest.get('code_files') or not manifest.get('source_snapshot'):raise ValueError('affinity dependency provenance missing')
    for name,digest in manifest['code_files'].items():
        path=(root/name).resolve()
        if not path.is_relative_to(root) or sha256(path)!=digest:raise ValueError('affinity builder/consumer code changed')
    for name,digest in manifest['files'].items():
        path=(release/name).resolve()
        if not path.is_relative_to(release.resolve()) or sha256(path)!=digest:raise ValueError('changed release artifact')
    snapshot=Path(manifest['source_snapshot']).resolve()
    for source in manifest['sources']:
        path=(snapshot/source['file']).resolve()
        if not path.is_relative_to(snapshot) or sha256(path)!=source['sha256']:raise ValueError('changed affinity raw source')
    if sha256(manifest['registry_path'])!=manifest['external_registry_sha256']:raise ValueError('external registry changed')
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--registry',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();print(json.dumps(build(a.snapshot,a.registry,a.out),indent=2))
if __name__=='__main__':main()
