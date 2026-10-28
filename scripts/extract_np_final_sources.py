"""Extract narrowly scoped final NP construct and paired-repertoire facts from publisher workbooks."""
from pathlib import Path
import argparse,hashlib,json,zipfile,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_boundaries_20261029'
def workbook_cells(path):
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        strings=[''.join(t.itertext()) for t in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)] if 'xl/sharedStrings.xml' in z.namelist() else []
        rel={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))};sheets={}
        for sheet in ET.fromstring(z.read('xl/workbook.xml')).findall('s:sheets/s:sheet',ns):
            target=rel[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']];member=target.lstrip('/') if target.startswith('/') else 'xl/'+target;cells={}
            for c in ET.fromstring(z.read(member)).findall('.//s:sheetData/s:row/s:c',ns):
                assert c.find('s:f',ns) is None,'formula cell requires separate treatment'
                v=c.find('s:v',ns)
                if v is not None:cells[c.attrib['r']]=strings[int(v.text)] if c.get('t')=='s' else v.text
                elif c.get('t')=='inlineStr':cells[c.attrib['r']]=''.join(c.find('s:is',ns).itertext())
            sheets[sheet.attrib['name']]=cells
    return sheets

def paired_repertoire_support(row,primary):
    return all(row[k]==primary[k] for k in ['CDR3_alpha','CDR3_beta','V_alpha','J_alpha','V_beta','J_beta'])
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    old=json.loads((ROOT/'docs/fingerprint_expanded_controls_20261029/natural_primary_cells.json').read_text());p=BASE/'sources/final/final_mmc7.xlsx';book=workbook_cells(p);new=json.loads(json.dumps(old))
    new.update(source_url='https://ars.els-cdn.com/content/image/1-s2.0-S2211124726009861-mmc7.xlsx',article_doi='10.1016/j.celrep.2026.117908',article_status='final publisher supplement; matched title and DOI metadata',workbook_member='final_mmc7.xlsx',workbook_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),source_license='publisher article metadata: CC BY-NC-ND 4.0; original workbook retained locally; factual sequence cells published')
    for row in new['rows']:
        n=str(row['source_row']);info=book['TCR constructs for Structures']
        for col,key in [('A','name'),('B','TRAV_TRAJ'),('C','TRBV_TRBJ'),('D','CDR3_alpha'),('E','CDR3_beta'),('H','constant_domain')]:row[key]=info[col+n]
        for c in row['chains'].values():c['primary_labeled_variable']=book[c['sheet']][c['cell']]
    oldbook=workbook_cells(ROOT/'docs/fingerprint_expanded_controls_20261029/sources/expanded/media-5.xlsx')
    keys=['name','TRAV_TRAJ','TRBV_TRBJ','CDR3_alpha','CDR3_beta','constant_domain']
    new['preprint_comparison']=dict(n_paired_receptors=6,n_recorded_variable_cells=12,variable_cells_unchanged=all(x['chains']==y['chains'] for x,y in zip(old['rows'],new['rows'])),CDR3_and_gene_and_constant_table_unchanged=all({k:x[k] for k in keys}=={k:y[k] for k in keys} for x,y in zip(old['rows'],new['rows'])),added_sheets=sorted(set(book)-set(oldbook)),preprint_workbook_sha256=old['workbook_sha256'])
    (a.out/'final_primary_cells.json').write_text(json.dumps(new,indent=2)+'\n')
    p=BASE/'sources/final/final_mmc3.xlsx';sheet='Data S2. Clonotype list';cells=workbook_cells(p)[sheet];rows=[]
    for n in sorted({int(''.join(c for c in k if c.isdigit())) for k in cells}):
        if n==1:continue
        n=str(n);rows.append(dict(clonotype=cells['A'+n],source_row=int(n),V_alpha=cells['C'+n],J_alpha=cells['D'+n],CDR3_alpha=cells['E'+n],V_beta=cells['G'+n],J_beta=cells['H'+n],CDR3_beta=cells['I'+n]))
    checks=[]
    for r in new['rows']:
        va,ja=r['TRAV_TRAJ'].rsplit('/',1);vb,jb=r['TRBV_TRBJ'].rsplit('/',1);target=dict(CDR3_alpha=r['CDR3_alpha'],CDR3_beta=r['CDR3_beta'],V_alpha=va,J_alpha=ja,V_beta=vb,J_beta=jb)
        hits=[x for x in rows if paired_repertoire_support(x,target)];beta_only=[x for x in rows if x['CDR3_beta']==target['CDR3_beta'] and not paired_repertoire_support(x,target)]
        checks.append(dict(name=r['name'],target=target,paired_CDR3_and_VJ_matches=hits,rejected_beta_only_or_context_mismatch=beta_only,n_matching_paired_rows=len(hits),interpretation='Paired CDR3/VJ supports primary lineage provenance only; multiple nucleotide clonotypes do not increase amino-acid receptor count; no independent donor identity inferred.'))
    report=dict(source_url='https://ars.els-cdn.com/content/image/1-s2.0-S2211124726009861-mmc3.xlsx',workbook_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),sheet=sheet,n_published_clonotypes=len(rows),rows=checks,independent_acceptance_accessed=False)
    (a.out/'final_repertoire_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print([(x['name'],x['n_matching_paired_rows']) for x in checks])
if __name__=='__main__':main()
