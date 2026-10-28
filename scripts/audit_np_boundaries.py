"""Reconcile NP domain boundaries, independent numbering, and fixed functional references."""
from pathlib import Path
import argparse,hashlib,json,re
from meta_acdc.structure.paired_variable_identity import germline_references,fasta_records,_CODONS
from scripts.curate_cross_source_receptors import chain_pairs
from meta_acdc.structure.paired_variable_identity import numbered_domains
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_boundaries_20261029'

def checked_numbering(full,numbering,start,end_inclusive,chain_type,expected):
    if chain_type!=expected:raise ValueError('numbered chain role differs from primary role')
    actual=[(tuple(p),aa) for p,aa in numbering if aa!='-']
    if not actual or len({p for p,aa in actual})!=len(actual):raise ValueError('missing or duplicate numbered residues')
    seq=''.join(aa for _,aa in actual)
    if start is None or end_inclusive is None or full[start:end_inclusive+1]!=seq or full.count(seq)!=1:raise ValueError('numbered domain is not the exact unique deposited interval')
    if min(p[0] for p,aa in actual)>3 or max(p[0] for p,aa in actual)<(128 if expected=='A' else 127):raise ValueError('terminal numbering incomplete')
    index=start;mapping=[]
    for position,aa in actual:
        mapping.append(dict(deposited_index0=index,IMGT_position=str(position[0])+str(position[1]).strip(),amino_acid=aa));index+=1
    return dict(start0=start,end0_exclusive=end_inclusive+1,length=len(seq),sequence=seq,first_IMGT=mapping[0]['IMGT_position'],last_IMGT=mapping[-1]['IMGT_position'],mapping=mapping)

def raw_j_references(paths):
    result=[]
    for path in paths:
        for header,seq in fasta_records(path.read_text()):
            p=header.split('|')
            if len(p)<8 or p[7] not in ['1','2','3']:continue
            nt=seq.replace('.','').upper()[int(p[7])-1:];aa=''.join(_CODONS.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))
            if '*' not in aa and 'X' not in aa:result.append(dict(gene=p[1],function=p[3],sequence=aa,header=header))
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    neural=json.loads((BASE/'numbering/numbering_raw.json').read_text());classic=json.loads((BASE/'numbering/classic_raw.json').read_text());facts=json.loads((BASE/'final_primary_cells.json').read_text())
    primary={r['name']+'_'+role:(r,c) for r in facts['rows'] for role,c in r['chains'].items()}
    paths=list((BASE/'sources/imgt').glob('*.fasta'));refs=germline_references(paths);rawj=raw_j_references([p for p in paths if p.stem.endswith(('TRAJ','TRBJ'))]);records=[]
    for i,item in enumerate(neural['inputs']):
        name=item['name'];expected={'alpha':'A','beta':'B'}[item['role']];row,cell=primary[name];full=item['sequence'];nr=neural['results'][name];cr=classic['numbered'][i];cd=classic['details'][i];methods={}
        if not nr.get('error'):
            try:methods['ANARCII']=dict(status='numbered',**checked_numbering(full,nr['numbering'],nr['query_start'],nr['query_end'],nr['chain_type'],expected))
            except ValueError as e:methods['ANARCII']=dict(status='unresolved',reason=str(e))
        else:methods['ANARCII']=dict(status='unresolved',reason=nr['error'],reported_chain_type=nr['chain_type'])
        if cr and len(cr)==1:
            numbered,start,end=cr[0]
            try:methods['ANARCI']=dict(status='numbered',**checked_numbering(full,numbered,start,end,cd[0]['chain_type'],expected),germline_nearest_match=cd[0].get('germlines'))
            except ValueError as e:methods['ANARCI']=dict(status='unresolved',reason=str(e))
        else:methods['ANARCI']=dict(status='unresolved',reason='no permitted alpha/beta numbering; best HMM classification may be delta',raw_hit_table=classic['hits'][i])
        valid=[v for v in methods.values() if v['status']=='numbered'];agreement=len(valid)==2 and all(valid[0][k]==valid[1][k] for k in ['start0','end0_exclusive','sequence','mapping'])
        boundary_agreement=len(valid)==2 and all(valid[0][k]==valid[1][k] for k in ['start0','end0_exclusive','sequence'])
        chosen=valid[0] if valid else None;vp='TRAV' if expected=='A' else 'TRBV';jp='TRAJ' if expected=='A' else 'TRBJ';sequence=chosen['sequence'] if chosen else None
        v_matches=sorted(g for (sp,g),rr in refs.items() if sp=='mouse' and g.startswith(vp) and sequence and any(sequence.startswith(x['amino_acids'][:10]) for x in rr))
        j_matches=sorted(g for (sp,g),rr in refs.items() if sp=='mouse' and g.startswith(jp) and sequence and any(sequence.endswith(x['amino_acids'][-10:]) for x in rr))
        raw_matches=[dict(r,exact_FR4_suffix_start0=len(sequence)-10,exact_FR4_suffix_end0_exclusive=len(sequence)) for r in rawj if sequence and r['gene'].startswith(jp) and sequence.endswith(r['sequence'][-10:])]
        text=(BASE/'sources/stcrdab'/(row['pdb']+'.html')).read_text();assert not chain_pairs(text) and not numbered_domains(text)
        primary_sequence=cell['primary_labeled_variable'];delta=primary_sequence[len(sequence):] if sequence and primary_sequence.startswith(sequence) else None
        record=dict(name=name,pdb=row['pdb'],role=item['role'],entity=item['entity'],author_chains=item['author_chains'],primary_cell=cell,methods=methods,two_method_position_exact_agreement=agreement,two_method_boundary_exact_agreement=boundary_agreement,
            primary_V_region_extra_after_numbered_domain=delta,numbered_boundary=({k:v for k,v in chosen.items() if k not in ['mapping','status','germline_nearest_match']} if chosen else None),
            functional_mouse_V_first10_matches=v_matches,functional_mouse_J_last10_matches=j_matches,all_function_status_J_last10_matches=raw_matches,
            fixed_terminal_reference_checks_pass=bool(chosen and v_matches and j_matches),
            V_J_rearrangement_junction='not uniquely inferred from construct amino acids; longest/suffix matches do not locate historical recombination cuts',
            external_STCRDab_pair_and_numbering_present=False,fixed_identity_admitted=False,
            reasons=['no valid external STCRDab full-domain paired annotation; independent tools are supplemental evidence']+([] if agreement else ['independent numbering does not agree at every numbered residue or one method failed'])+([] if j_matches else ['no permitted functional J/FR4 reference match; ORF matches do not pass the existing gate']))
        records.append(record)
    sourcehash={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result=dict(n_chains=12,n_ANARCII_numbered=sum(r['methods']['ANARCII']['status']=='numbered' for r in records),n_ANARCI_numbered=sum(r['methods']['ANARCI']['status']=='numbered' for r in records),n_two_method_exact_agreements=sum(r['two_method_position_exact_agreement'] for r in records),n_two_method_boundary_agreements=sum(r['two_method_boundary_exact_agreement'] for r in records),n_functional_terminal_checks_pass=sum(r['fixed_terminal_reference_checks_pass'] for r in records),n_admitted_chains=0,rows=records,reference_sha256=sourcehash,final_vs_preprint=facts['preprint_comparison'],protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k.startswith('n_')},indent=2))
if __name__=='__main__':main()
