"""Recheck all-functionality reference claims without changing frozen gates."""
from pathlib import Path
import argparse,hashlib,json,re
from meta_acdc.structure.j_reference_evidence import read_j_evidence,suffix_evidence
from meta_acdc.structure.paired_variable_identity import numbered_domains,germline_references
from scripts.audit_np_boundaries import raw_j_references
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_tk3_j_evidence_20261101'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def audit():
    refs_dir=ROOT/'docs/fingerprint_paired_identity_20261024/sources'
    paths=sorted(refs_dir.glob('*_TRAJ.fasta'))+sorted(refs_dir.glob('*_TRBJ.fasta'));raw=read_j_evidence(paths);old=raw_j_references(paths);permitted=germline_references(paths)
    previous=ROOT/'docs/fingerprint_control_closure_20261101/audit_final_v3.json';case=json.loads(previous.read_text())['independent_4PRH_numbering']['rows'];cases=[]
    for row in case:
        prefix='TRAJ' if row['role']=='alpha' else 'TRBJ';seq=row['observed_variable'];matches=suffix_evidence(seq,raw,'Homo sapiens',prefix)
        cases.append(dict(role=row['role'],variable=seq,old_reported_all_status_matches=row['all_function_status_human_or_mouse_J_exact_matches'],corrected_all_status_suffix_evidence=matches,functional_gate_matches=row['functional_human_J_matches'],fixed_identity_admitted=False))
    data=json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/source_balanced_design.json').read_text());seen=set();rows=[];skips=[]
    for row in data['rows']:
        pdb=row['pdb'];p=ROOT/'docs/fingerprint_cross_source_extension_20261024/sources/stcrdab_all'/(pdb+'.html')
        if not p.exists():p=ROOT/'docs/fingerprint_cross_source_20261024/sources/stcrdab'/(pdb+'.html')
        if not p.exists():skips.append(dict(pdb=pdb,reason='no stored annotation page'));continue
        text=p.read_text();domains=numbered_domains(text);species=set(re.findall(r'<b>Species</b></td><td[^>]*>([^<]+)</td>',text))
        if len(species)!=1:skips.append(dict(pdb=pdb,reason='ambiguous or missing annotation species'));continue
        sp=next(iter(species)).strip();fullsp={'human':'Homo sapiens','mouse':'Mus musculus'}.get(sp)
        if not fullsp:skips.append(dict(pdb=pdb,reason='outside reference species'));continue
        for role,kind,prefix in [('alpha','VA','TRAJ'),('beta','VB','TRBJ')]:
            chain=row[role+'_chain'];key=(pdb,chain,role)
            if key in seen:continue
            seen.add(key);ann=domains.get(chain)
            if not ann or ann['kind']!=kind:continue
            seq=ann['sequence'];construct=(row.get(role) or {}).get('full_construct','')
            if construct.count(seq)!=1:continue
            new=suffix_evidence(seq,raw,fullsp,prefix);before=[x for x in old if x['gene'].startswith(prefix) and '|'+fullsp+'|' in x['header'] and seq.endswith(x['sequence'][-10:])];functional=[g for (s,g),vv in permitted.items() if s==sp and g.startswith(prefix) and any(seq.endswith(x['amino_acids'][-10:]) for x in vv)]
            rows.append(dict(pdb=pdb,chain=chain,role=role,observed_numbered_variable=seq,annotation_file=str(p.relative_to(ROOT)),annotation_sha256=sha(p),source_study=row.get('source_study_key'),functional_stopfree_J_matches=sorted(functional),old_stopfree_any_status_matches=sorted({x['gene'] for x in before}),unfiltered_suffix_evidence=new,newly_visible_stopped_or_ambiguous_references=[x for x in new if x['old_stopfree_filter_would_drop']],identity_admitted=False))
    affected=[r for r in rows if r['newly_visible_stopped_or_ambiguous_references']]
    result=dict(date='2026-11-01',n_reference_records=len(raw),n_stop_containing_references=sum(bool(r['stop_positions0']) for r in raw),n_ambiguous_references=sum(bool(r['ambiguous_positions0']) for r in raw),TK3_cases=cases,n_observed_external_domains=len(rows),n_domains_with_newly_visible_references=len(affected),affected_PDBs=sorted({r['pdb'] for r in affected}),affected_genes=sorted({x['gene'] for r in affected for x in r['newly_visible_stopped_or_ambiguous_references']}),rows=rows,unreviewed_annotation_skips=skips,reference_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},fresh_human_TRAJ_byte_exact=(BASE/'sources/IMGT_Homo_sapiens_TRAJ.fasta').read_bytes()==(refs_dir/'Homo_sapiens_TRAJ.fasta').read_bytes(),old_helper_remains_for_frozen_replay=True,admission_gates_changed=False,new_identity_admissions=0,new_QC_cross_study_pairs=0,protocol_frozen=False,training_started=False,independent_acceptance_accessed=False)
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    r=audit();a.out.write_text(json.dumps(r,indent=2)+'\n');print({k:v for k,v in r.items() if k.startswith('n_') or k.startswith('affected_')});print('TK3',[(x['role'],[(m['gene'],m['longest_exact_suffix_length']) for m in x['corrected_all_status_suffix_evidence']]) for x in r['TK3_cases']])
if __name__=='__main__':main()
