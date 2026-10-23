"""Curate complete deposited variable pairs; quarantine all unsupported correspondences."""
import argparse,hashlib,json,re
from pathlib import Path
from collections import Counter,defaultdict
import torch
from meta_acdc.structure.paired_variable_identity import germline_references,numbered_domains,validate_domain,pair_key,cif_declared_sequences,fasta_records
from meta_acdc.structure.receptor_identity_audit import is_subsequence
ROOT=Path(__file__).resolve().parents[1]
MAIN=Path('/home/administrator/docs/meta-acdc')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out-dir',type=Path,required=True);a=ap.parse_args();out=a.out_dir
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);base=out.parent;src=base/'sources'
    old=ROOT/'docs/fingerprint_retrain_20261023/input/fingerprint_dataset.pt';data=torch.load(old,map_location='cpu',weights_only=False)
    audit=json.loads((ROOT/'docs/fingerprint_hard_family_20261024/identity_audit/identity_audit.json').read_text())
    def source(name):
        retry=src/'retry1'/name
        manifest=json.loads((src/'retry1/manifest.json').read_text())
        return retry if any(x['file']==name and x['status']=='200' for x in manifest) else src/name
    refs=germline_references(list(src.glob('Homo_sapiens_TR*.fasta'))+list(src.glob('Mus_musculus_TR*.fasta')))
    crystal={x['index']:x for x in audit['rows']};rows=[];valid={}
    for index,oldrow in crystal.items():
        row=dict(index=index,complex_id=oldrow['complex_id'],domain='crystal',old_receptor_id=oldrow['observed_sequence_id'],old_family_id=oldrow['family_id'],status='quarantined',reasons=[])
        try:
            if not oldrow['eligible_construct_evidence']:raise ValueError('original chain role or SEQRES correspondence invalid')
            text=source(oldrow['pdb']+'.html').read_text();domains=numbered_domains(text)
            role_errors=[]
            for role,kind in [('alpha','VA'),('beta','VB')]:
                chain=oldrow[role+'_chain'];ann=domains.get(chain)
                path=source(oldrow['pdb']+'_'+chain+'.fasta')
                row[role]=dict(full_construct=oldrow['seqres_'+role],observed=data['records'][index][role],chain=chain,domain_candidate=ann)
                try:
                    if ann is None or ann['kind']!=kind:raise ValueError('external domain role missing')
                    item=validate_domain(oldrow['seqres_'+role],ann,path.read_text(),text,refs)
                    row[role].update(item,annotation_fasta=str(path.relative_to(ROOT)),annotation_sha256=sha(path),boundary_status='pass')
                except (ValueError,OSError) as error:
                    row[role]['boundary_status']='quarantined';role_errors.append(role+': '+str(error))
            if role_errors:raise ValueError('; '.join(role_errors))
            # Depositor chains were assigned by unique alpha/beta contacts in the frozen inventory.
            row['receptor_id']=pair_key(row['alpha']['variable'],row['beta']['variable'])
            row['construct_id']=pair_key(row['alpha']['full_construct'],row['beta']['full_construct'])
            row['status']='admitted';row['source_pdb_sha256']=oldrow['source_sha256'];valid[index]=row
        except (ValueError,KeyError,OSError) as e:row['reasons'].append(str(e))
        rows.append(row)
    # Explicit role FASTAs are inputs, not an assumption that PDB-name equality proves identity.
    role_refs=defaultdict(list)
    for path in sorted((ROOT/'af3_pending_upload').rglob('*.fasta')):
        fr=dict(fasta_records(path.read_text()))
        if set(fr)=={'alpha','beta','peptide','mhc','b2m'}:
            role_refs[(fr['alpha'],fr['beta'])].append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path),peptide=fr['peptide']))
    inv=torch.load(MAIN/'data/releases/20261017-structural-v7-atomic-repair/native_inventory_graphs.pt',map_location='cpu',weights_only=False)
    inventory={x['record']['complex_id']:x['record'] for x in inv['assets']}
    excerpt=out/'cif_declared_sequences';excerpt.mkdir()
    for index,record in enumerate(data['records']):
        if index in crystal:continue
        row=dict(index=index,complex_id=record['complex_id'],domain='predicted',old_receptor_id=record['receptor_id'],old_family_id=record['family_id'],status='quarantined',reasons=[])
        try:
            if record['source']=='v7-inventory':
                asset=inventory[record['complex_id']];path=Path(asset['source_path']);expected=asset['source_sha256']
            else:path=MAIN/'data/raw/af3_predictions'/(record['complex_id']+'_model_0.cif');expected=data['graph_hashes'][index]
            if sha(path)!=expected:raise ValueError('CIF source hash differs from frozen graph provenance')
            declared=cif_declared_sequences(path.read_text());pair=(record['alpha'],record['beta'])
            mappings=[(ca,cb) for ca,sa in declared.items() for cb,sb in declared.items() if ca!=cb and (sa,sb)==pair]
            if len(mappings)!=1:raise ValueError('declared input/observed paired chain mismatch or ambiguous')
            proof=role_refs.get(pair,[])
            if not proof:raise ValueError('no exact explicit alpha/beta role FASTA input')
            ca,cb=mappings[0]
            candidates={v['receptor_id']:v for v in valid.values() if pair[0].count(v['alpha']['variable'])==1 and pair[1].count(v['beta']['variable'])==1}
            if len(candidates)!=1:raise ValueError('complete paired variable sequence unavailable or ambiguous; no nearest match')
            parent=next(iter(candidates.values()))
            for role,seq,chain in [('alpha',pair[0],ca),('beta',pair[1],cb)]:
                variable=parent[role]['variable'];start=seq.index(variable)
                row[role]=dict(variable=variable,v_allele=parent[role]['v_allele'],full_construct=seq,observed=record[role],chain=chain,variable_start0=start,variable_end0=start+len(variable),n_terminal_extension=seq[:start],c_terminal_extension=seq[start+len(variable):],evidence_parent_crystal_index=parent['index'])
            row.update(receptor_id=parent['receptor_id'],construct_id=pair_key(*pair),source_cif_sha256=expected,exact_role_fasta_evidence=proof,status='admitted')
            (excerpt/(str(index)+'.json')).write_text(json.dumps(dict(index=index,source_path=str(path),source_sha256=expected,declared_sequences=declared,chain_mapping=dict(alpha=ca,beta=cb)),indent=2)+'\n')
        except (ValueError,KeyError,OSError) as e:row['reasons'].append(str(e))
        rows.append(row)
    rows.sort(key=lambda x:x['index']);assert [x['index'] for x in rows]==list(range(len(data['records'])))
    # Retain old family links and add exact complete single-chain links. Never split known relatives.
    parents=list(range(len(rows)))
    def find(x):
        while parents[x]!=x:parents[x]=parents[parents[x]];x=parents[x]
        return x
    def union(a,b):parents[find(b)]=find(a)
    groups={}
    for row in rows:
        if row['status']!='admitted':continue
        keys=[('old_family',row['old_family_id']),('paired',row['receptor_id']),('alpha',row['alpha']['variable']),('beta',row['beta']['variable'])]
        for k in keys:
            if k in groups:union(row['index'],groups[k])
            else:groups[k]=row['index']
    components=defaultdict(list)
    for row in rows:
        if row['status']=='admitted':components[find(row['index'])].append(row['index'])
    fold=json.loads((ROOT/'docs/fingerprint_retrain_20261023/fold_assignment.json').read_text());foldof={family:i for i,values in enumerate(fold['family_folds']) for family in values}
    for indices in components.values():
        cid='C'+hashlib.sha256('|'.join(sorted({rows[i]['receptor_id'] for i in indices})).encode()).hexdigest()[:16]
        folds={foldof[rows[i]['old_family_id']] for i in indices}
        crosses=(0 in folds) or (1 in folds and bool(folds&{2,3,4}))
        for i in indices:
            rows[i]['family_id']=cid;rows[i]['original_fold']=foldof[rows[i]['old_family_id']]
            rows[i]['development_partition']='quarantined_partition' if crosses else 'validation' if 1 in folds else 'fit'
            if crosses:rows[i]['partition_reason']='entire linked component touches old fold0 or crosses fixed fit/validation boundary; no reassignment'
    admitted=[x for x in rows if x['status']=='admitted'];development=[x for x in admitted if x['development_partition'] in ['fit','validation']]
    summary=dict(n_original=len(rows),n_identity_admitted=len(admitted),n_identity_quarantined=len(rows)-len(admitted),n_development=len(development),counts_by_domain_status=dict(Counter(x['domain']+'/'+x['status'] for x in rows)),n_complete_paired_variable_ids=len({x['receptor_id'] for x in admitted}),n_development_paired_variable_ids=len({x['receptor_id'] for x in development}),partitions=dict(Counter(x['development_partition'] for x in admitted)),quarantine_reasons=dict(Counter(reason for x in rows for reason in x['reasons'])))
    result=dict(summary=summary,identity_semantics='exact deposited complete paired variable sequences; engineered substitutions preserved; not clinical clone identity or binding truth',admission='both canonical domain boundaries + complete source construct + chain role proof; no missing-residue imputation',rows=rows,input_sha256=sha(old),reference_source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted(src.rglob('*')) if p.is_file()},independent_acceptance_accessed=False)
    (out/'identity_registry.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
