"""Offline replay and preservation checks for the expanded public development audit."""
from pathlib import Path
import argparse,gzip,hashlib,json,os,subprocess,sys,zipfile,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'docs/fingerprint_expanded_controls_20261029'
def sha(raw):return hashlib.sha256(raw).hexdigest()
def verify_local_workbook(facts):
    path=BASE/'sources/expanded/media-5.xlsx'
    if not path.exists():return 'not_packaged; obtain original from recorded primary URL for workbook re-extraction'
    assert sha(path.read_bytes())==facts['workbook_sha256']
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        strings=[''.join(t.itertext()) for t in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
        rel={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheets={}
        for sheet in ET.fromstring(z.read('xl/workbook.xml')).findall('s:sheets/s:sheet',ns):
            target=rel[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
            member=target.lstrip('/') if target.startswith('/') else 'xl/'+target
            sheets[sheet.attrib['name']]=ET.fromstring(z.read(member))
        def cell(sheet,reference):
            c=sheets[sheet].find('.//s:c[@r="'+reference+'"]',ns);assert c is not None and c.find('s:f',ns) is None
            v=c.find('s:v',ns);assert v is not None
            return strings[int(v.text)] if c.get('t')=='s' else v.text
        for row in facts['rows']:
            n=str(row['source_row'])
            for col,key in [('A','name'),('B','TRAV_TRAJ'),('C','TRBV_TRBJ'),('D','CDR3_alpha'),('E','CDR3_beta'),('H','constant_domain')]:assert cell('TCR constructs for Structures',col+n)==row[key]
            for chain in row['chains'].values():assert cell(chain['sheet'],chain['cell'])==chain['primary_labeled_variable']
    return 'local_original_hash_and_all_recorded_cells_exact'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--replay-dir',type=Path,required=True);args=ap.parse_args();out=args.replay_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,PYTHONPATH='src:.')
    def run(script,*extra):
        p=subprocess.run([sys.executable,'scripts/'+script,*map(str,extra)],cwd=ROOT,env=env,capture_output=True,text=True)
        with (out/'replay.log').open('a') as f:f.write(p.stdout+p.stderr)
        assert p.returncode==0,p.stderr
    archive=json.loads((BASE/'archive_layout.json').read_text());public=0;local=0
    for row in archive:
        if row['packaging']=='local_only_primary_NC':
            path=BASE/row['source']
            if path.exists():assert len(path.read_bytes())==row['bytes'] and sha(path.read_bytes())==row['sha256'];local+=1
            continue
        raw=(BASE/row['packaged_file']).read_bytes()
        if row['packaging']=='gzip_exact':raw=gzip.decompress(raw)
        assert len(raw)==row['bytes'] and sha(raw)==row['sha256'],row['source'];public+=1
    run('prepare_expanded_controls.py','--out',out/'plan')
    for name in ['search_plan.json','requests.json']:assert (out/'plan'/name).read_bytes()==(BASE/'complete_search'/name).read_bytes(),name
    for group in ['complete_search','natural_search']:
        run('join_native_search.py','--base',BASE/group,'--out',out/(group+'_results.json'))
        assert (out/(group+'_results.json')).read_bytes()==(BASE/group/'search_results.json').read_bytes()
    run('validate_expanded_controls.py','--out',out/'primary_validation_both_ends.json')
    assert (out/'primary_validation_both_ends.json').read_bytes()==(BASE/'complete_search/primary_validation_both_ends.json').read_bytes()
    run('audit_expanded_controls.py','--out',out/'audit_final.json')
    assert (out/'audit_final.json').read_bytes()==(BASE/'audit_final.json').read_bytes()
    facts=json.loads((BASE/'natural_primary_cells.json').read_text());workbook=verify_local_workbook(facts)
    plan=json.loads((BASE/'natural_search/search_plan.json').read_text())
    for target,row in zip(plan['targets'],facts['rows']):assert target['sequences']=={k:v['primary_labeled_variable'] for k,v in row['chains'].items()} and target['full_boundary_validated'] is False
    previous=[]
    folders=['fingerprint_retrain_20261023','fingerprint_hard_family_20261024','fingerprint_paired_identity_20261024','fingerprint_cross_source_20261024','fingerprint_cross_source_extension_20261024','fingerprint_audit_integration_20261027','fingerprint_positive_controls_20261027','fingerprint_construct_followup_20261027','fingerprint_native_search_20261029']
    for folder in folders:
        base=ROOT/'docs'/folder;manifest=json.loads((base/'DELIVERY_MANIFEST.json').read_text())
        for name,expected in manifest['files'].items():assert sha((base/name).read_bytes())==expected['sha256'],(folder,name)
        previous.append(dict(folder=folder,n_files=len(manifest['files']),unchanged=True))
    gates=json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/DELIVERY_MANIFEST.json').read_text())
    for name,expected in {**gates['unchanged_gates'],'src/meta_acdc/structure/primary_role_evidence.py':gates['code']['src/meta_acdc/structure/primary_role_evidence.py']}.items():assert sha((ROOT/name).read_bytes())==expected['sha256'],name
    modules=['tests.test_construct_role_adjudication','tests.test_primary_role_evidence','tests.test_cross_source_extension','tests.test_cross_source_retrieval','tests.test_paired_variable_identity','tests.test_native_search','tests.test_expanded_controls']
    p=subprocess.run([sys.executable,'-m','unittest',*modules,'-v'],cwd=ROOT,env=env,capture_output=True,text=True);(out/'tests.log').write_text(p.stdout+p.stderr);assert p.returncode==0 and 'Ran 51 tests' in p.stderr and 'skipped=' not in p.stderr,p.stderr
    result=dict(status='pass',n_tests=51,test_skips=0,n_public_raw_response_payloads_verified=public,n_local_only_primary_payloads_hash_checked=local,primary_workbook_check=workbook,plan_search_primary_and_summary_byte_exact=True,previous_deliveries=previous,n_frozen_files=sum(p['n_files'] for p in previous),role_boundary_atomic_gate_code_unchanged=True,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
