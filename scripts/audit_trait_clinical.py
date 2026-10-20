"""Read-only clinical-source triage. CDR3 candidates never become clinical risk labels."""
from __future__ import annotations
import argparse, hashlib, io, json, re, zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def rows_from_zip(path):
    with zipfile.ZipFile(path) as outer:
        names=[n for n in outer.namelist() if n.endswith('.xlsx')]
        if len(names)!=1:raise ValueError('expected one clinical workbook')
        content=outer.read(names[0])
    with zipfile.ZipFile(io.BytesIO(content)) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            root=ET.fromstring(z.read('xl/sharedStrings.xml'))
            strings=[''.join(x.itertext()) for x in root.findall('s:si',NS)]
        sheets=[x for x in z.namelist() if re.fullmatch(r'xl/worksheets/sheet[0-9]+.xml',x)]
        if len(sheets)!=1:raise ValueError('review multi-sheet source explicitly')
        root=ET.fromstring(z.read(sheets[0]));rows=[]
        for row in root.findall('.//s:sheetData/s:row',NS):
            values=[None]*23
            for cell in row.findall('s:c',NS):
                if cell.find('s:f',NS) is not None:raise ValueError('formula source requires explicit evaluation')
                letters=re.match('[A-Z]+',cell.attrib['r']).group();index=0
                for c in letters:index=index*26+ord(c)-64
                value=cell.find('s:v',NS);text=value.text if value is not None else None
                if cell.attrib.get('t')=='s':text=strings[int(text)]
                elif cell.attrib.get('t')=='inlineStr':text=''.join(cell.find('s:is',NS).itertext())
                if index>23:raise ValueError('source layout changed')
                values[index-1]=text
            rows.append((int(row.attrib['r']),values))
    expected={0:'Clinical trials',2:'TCR name',3:'TRAV',6:'TRBV',11:'Target antigen',12:'Epitope',21:'PMID'}
    if any(rows[0][1][k]!=v for k,v in expected.items()):raise ValueError('unreviewed source header')
    return rows[1:]


def audit(path,exposure_records):
    canon=lambda s:re.sub(r'\s+','',str(s or '')).upper()
    present=lambda s:bool(s) and str(s).strip().lower() not in {'na','none','-','n/a'}
    named=[(n,r) for n,r in rows_from_zip(path) if present(r[2])]
    observations=[];pairs=set();paired_rows=0
    for n,r in named:
        a,b=canon(r[5]),canon(r[9]);strict=bool(re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]{8,35}',a) and re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]{8,35}',b))
        overlaps=[]
        if strict:
            pairs.add((a,b));paired_rows+=1
            for source,x in exposure_records:
                alpha=x.get('alpha_variable',x.get('alpha',''));beta=x.get('beta_variable',x.get('beta',''))
                # Conservative candidate exclusion: either known chain can flag exposure.
                if (a and a in alpha) or (b and b in beta):overlaps.append(source+':'+str(x.get('receptor_name',x.get('complex_id',x.get('evidence_id','unknown')))))
        observations.append(dict(excel_row=n,name=r[2],antigen=r[11],epitope=r[12],hla=r[13],
            cdr3_alpha=r[5],cdr3_beta=r[9],cdr_modification=r[10],trial=r[0],pmid=r[21],
            canonical_paired_cdr3=strict,known_chain_exposures=sorted(set(overlaps)),
            full_construct_identity_verified=False,clinical_label=None))
    return dict(source_url='https://pgx.zju.edu.cn/download.trait/Therapeutics.zip',
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),named_rows=len(named),
        strict_canonical_paired_rows=paired_rows,distinct_canonical_cdr3_pairs=len(pairs),
        strict_rows_with_known_chain_exposure=sum(bool(x['known_chain_exposures']) for x in observations),
        full_alpha_beta_sequence_columns_present=False,independent_receptor_count=None,
        clinical_training_allowed=False,scope='source_triage_not_model_validation',observations=observations,
        quarantine_examples=[dict(excel_row=6,reason='gp100 row labelled DMF4 with MART-1 CDR3 pair; PMID19451549 distinguishes DMF5 and murine gp100154')])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--zip',type=Path,required=True)
    ap.add_argument('--structures',type=Path,required=True);ap.add_argument('--affinity',type=Path,required=True)
    ap.add_argument('--registry',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():raise ValueError('refusing to overwrite source audit')
    exposure=[('structure-v4',r) for r in json.loads(args.structures.read_text())]
    exposure += [('affinity-v3',r) for r in json.loads(args.affinity.read_text())]
    exposure += [('reserved-panels',r) for r in json.loads(args.registry.read_text())['receptors']]
    report=audit(args.zip,exposure)
    report['exposure_file_sha256']={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [args.structures,args.affinity,args.registry]}
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in {'observations','exposure_file_sha256'}},indent=2))


if __name__=='__main__':main()
