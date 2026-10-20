"""Freeze publisher-annotated MIRA activation benchmark, never SPR binding labels."""
from __future__ import annotations
import argparse,csv,hashlib,json,re
from collections import Counter,defaultdict
from pathlib import Path
from .experimental_evidence import AA
from .release import sha256


def validate_rows(raw):
    rows=[];seen=set();by_tcr=defaultdict(list)
    for line,x in enumerate(raw,2):
        alpha,beta=x['tcra_trimmed'],x['tcrb_trimmed'];pep=x['peptide'];hla=x['hla'];label=x['label']
        if not set(alpha+beta+pep)<=AA or min(len(alpha),len(beta))<70:raise ValueError('invalid paired canonical variable sequence')
        if not x['tcra_cdr3'] or not x['tcrb_cdr3'] or x['tcra_cdr3'] not in alpha or x['tcrb_cdr3'] not in beta:raise ValueError('CDR3 mismatch')
        if label not in {'0','1'} or not re.fullmatch(r'[AB]\*\d{2}:\d{2}',hla):raise ValueError('label/HLA scope')
        if len(pep)!=9:raise ValueError('not a benchmark 9-mer')
        pair=hashlib.sha256((alpha+'|'+beta).encode()).hexdigest();key=(pair,pep,hla)
        if key in seen:raise ValueError('duplicate receptor-peptide pair')
        seen.add(key);by_tcr[pair].append((pep,hla,int(label)))
        rows.append({'evidence_id':'IMMREP25:'+str(line),'source_line':line,'receptor_pair_id':pair,'alpha_variable':alpha,'beta_variable':beta,'cdr3_alpha':x['tcra_cdr3'],'cdr3_beta':x['tcrb_cdr3'],'peptide':pep,'mhc_allele':'HLA-'+hla,'endpoint':'MIRA_pool_activation','observed_activation':int(label),'binding_label':None,'clinical_label':None,'negative_basis':'conditional nonactivation inferred from multiplexed pool activation pattern; not direct SPR or universal nonbinding','hla_basis':'NetMHCpan prediction constrained by donor HLA genotype, per publisher methods','sequence_basis':'publisher-supplied complete variable domains and pairSEQ chain pairing; not independently resequenced','use':'external_functional_benchmark_reserved; exposure audit required before scoring'})
    for group in by_tcr.values():
        if len(group)!=10 or sum(x[2] for x in group)!=1 or len({x[1] for x in group})!=1:raise ValueError('publisher one-positive/nine-negative same-HLA matrix invalid')
    return rows


def build(snapshot,output):
    sources=json.loads((snapshot/'retrieval_manifest.json').read_text());metadata=json.loads((snapshot/'figshare_30491582.json').read_text());files={x['name']:x for x in metadata['files']};required=['immrep2025_for_release.tsv','immrep25_dataset_and_methods.pdf','figshare_30491582.json']
    for name in required:
        item=next(x for x in sources if x['file']==name);raw=(snapshot/name).read_bytes()
        if sha256(snapshot/name)!=item['sha256']:raise ValueError('changed functional source')
        if name in files and (len(raw)!=files[name]['size'] or hashlib.md5(raw).hexdigest()!=files[name]['computed_md5']):raise ValueError('publisher checksum mismatch')
    with (snapshot/'immrep2025_for_release.tsv').open(newline='') as f:rows=validate_rows(csv.DictReader(f,delimiter='\t'))
    output.mkdir(parents=True,exist_ok=False)
    (output/'functional_observations.json').write_text(json.dumps(rows,indent=2)+'\n');unique={x['receptor_pair_id']:x for x in rows};registry={'purpose':'reserved external functional test; do not merge into physical-binding/clinical training labels','receptors':[dict(receptor_name='IMMREP25:'+k,**{key:x[key] for key in ['alpha_variable','beta_variable','cdr3_alpha','cdr3_beta']}) for k,x in sorted(unique.items())]};(output/'reservation_registry.json').write_text(json.dumps(registry,indent=2)+'\n')
    report={'records':len(rows),'paired_receptors':len(unique),'peptides':len({x['peptide'] for x in rows}),'hla_alleles':sorted({x['mhc_allele'] for x in rows}),'activation_counts':dict(Counter(x['observed_activation'] for x in rows)),'endpoint':'MIRA_pool_activation','use':'external functional benchmark candidate only; not a training source in this run','physical_binding_validation_allowed':False,'clinical_validation_allowed':False,'license':metadata['license'],'raw_pool_measurements_available_in_export':False,'individual_donor_ids_available_in_export':False,'limitations':['One source/platform; 20 selected viral peptides and two HLA contexts','Conditional pool nonactivation is not proof of absent physical binding','HLA assignment uses prediction plus donor genotype','No dose/detection curves or patient toxicity outcomes','Related receptor families and peptide exposure must be excluded from any training sources before scoring']};(output/'summary.json').write_text(json.dumps(report,indent=2)+'\n');(output/'manifest.json').write_text(json.dumps({'task':'external_MIRA_activation_benchmark','sources':[x for x in sources if x['file'] in required],'files':{p.name:sha256(p) for p in output.iterdir()},'code_sha256':sha256(__file__),'training_used':False},indent=2)+'\n');return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();print(json.dumps(build(a.snapshot,a.out),indent=2))
if __name__=='__main__':main()
