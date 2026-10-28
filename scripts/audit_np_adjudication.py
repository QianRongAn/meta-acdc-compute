"""Resolve source-level NP chain roles and FR1 gaps while retaining fixed admission gates."""
from pathlib import Path
import argparse,collections,csv,gzip,hashlib,io,json
from meta_acdc.structure.paired_variable_identity import fasta_records,_CODONS
from scripts.audit_np_boundaries import raw_j_references
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_adjudication_20261029';PREV=ROOT/'docs/fingerprint_np_boundaries_20261029'
def csv_rows(path):return list(csv.DictReader(io.StringIO(gzip.decompress(path.read_bytes()).decode())))
def translate(nt):
    if len(nt)%3:return None
    return ''.join(_CODONS.get(nt[i:i+3],'X') for i in range(0,len(nt),3))
def high_quality(row):return all(row.get(k)=='True' for k in ['is_cell','high_confidence','full_length','productive'])
def unambiguous_pair(rows,alpha,beta,require_beta_constant=True):
    active=[r for r in rows if high_quality(r)]
    if len(active)!=2 or {r['chain'] for r in active}!={'TRA','TRB'}:return None
    if len({r['barcode'] for r in active})!=1:return None
    a=next(r for r in active if r['chain']=='TRA');b=next(r for r in active if r['chain']=='TRB')
    if a['cdr3']!=alpha or b['cdr3']!=beta:return None
    if not a.get('raw_clonotype_id') or a['raw_clonotype_id'] in ['None',''] or a['raw_clonotype_id']!=b['raw_clonotype_id']:return None
    if a.get('c_gene')!='TRAC' or not a.get('j_gene','').startswith('TRAJ') or not b.get('j_gene','').startswith('TRBJ'):return None
    if require_beta_constant and b.get('c_gene') not in ['TRBC1','TRBC2']:return None
    if any(translate(r['cdr3_nt'])!=r['cdr3'] for r in [a,b]):return None
    return a,b

def v_framework_references(paths):
    result=[]
    for path in paths:
        for header,nt in fasta_records(path.read_text()):
            h=header.split('|')
            if len(h)<8 or h[3] not in ['F','(F)'] or h[7] not in ['1','2','3']:continue
            nt=nt.upper()[int(h[7])-1:];triplets=[nt[i:i+3] for i in range(0,len(nt)-2,3)]
            aa=''.join('-' if t=='...' else _CODONS.get(t,'X') for t in triplets)
            if 'X' in aa or '*' in aa:continue
            # Alpha reference alignments can carry extra columns: do not equate every raw codon slot to an IMGT integer.
            anchors=[i for i,a in enumerate(aa) if a=='C' and 99<=i<=109]
            if not anchors:continue
            end=anchors[0]+1;prefix=aa[:end].replace('-','')
            result.append(dict(gene=h[1],header=header,framework_through_conserved_Cys=prefix,gapped_amino_acids=aa,raw_Cys_slot=end,
                raw_slots_directly_anchored_for_beta=h[1].startswith('TRBV') and len(aa)>103 and aa[22]=='C' and aa[103]=='C' and end==104))
    return result

def suffix_length(a,b):
    n=0
    for x,y in zip(reversed(a),reversed(b)):
        if x!=y:break
        n+=1
    return n

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    facts=json.loads((PREV/'final_primary_cells.json').read_text());numbered=json.loads((PREV/'boundary_audit_both_methods.json').read_text());byname={r['name']:r for r in numbered['rows']}
    refs=v_framework_references(list((PREV/'sources/imgt').glob('*TR?V.fasta')));jrefs=raw_j_references([PREV/'sources/imgt/Mus_musculus_TRAJ.fasta',PREV/'sources/imgt/Mus_musculus_TRBJ.fasta'])
    raw_inputs={r['name']:r for r in json.loads((PREV/'numbering/numbering_raw.json').read_text())['inputs']}
    role_checks=[];pair_records=[]
    samples=[('GSM9108719','NP1'),('GSM9108721','NP2'),('GSM9108725','D4')]
    for gsm,tag in samples:
        path=BASE/'sources/contigs'/f'{gsm}_{tag}_filtered_contig_annotations.csv.gz';rows=csv_rows(path);bybar=collections.defaultdict(list)
        for row in rows:bybar[row['barcode']].append(row)
        clones={r['clonotype_id']:r for r in csv_rows(BASE/'sources/contigs'/f'{gsm}_{tag}_clonotypes.csv.gz')}
        assert not {'sequence','sequence_aa','fwr1','fwr1_nt'}&set(rows[0]),'full native sequence availability changed; inspect source schema'
        for target in facts['rows']:
            candidate=[(bar,rr) for bar,rr in bybar.items() if target['CDR3_alpha'] in {r['cdr3'] for r in rr} and target['CDR3_beta'] in {r['cdr3'] for r in rr}]
            accepted=[];role_supported=[];rejected=[]
            for barcode,rr in candidate:
                pair=unambiguous_pair(rr,target['CDR3_alpha'],target['CDR3_beta'],require_beta_constant=False)
                if pair is None:rejected.append(barcode);continue
                alpha,beta=pair;clone=clones[alpha['raw_clonotype_id']];expected={'TRA:'+target['CDR3_alpha'],'TRB:'+target['CDR3_beta']}
                if set(clone['cdr3s_aa'].split(';'))!=expected:rejected.append(barcode);continue
                va,ja=target['TRAV_TRAJ'].rsplit('/',1);vb,jb=target['TRBV_TRBJ'].rsplit('/',1)
                if (alpha['v_gene'],alpha['j_gene'],beta['v_gene'],beta['j_gene'])!=(va,ja,vb,jb):rejected.append(barcode);continue
                keys=['contig_id','chain','v_gene','j_gene','c_gene','cdr3','cdr3_nt','raw_clonotype_id','raw_consensus_id']
                item=dict(barcode=barcode,alpha={k:alpha[k] for k in keys},beta={k:beta[k] for k in keys});role_supported.append(item)
                if beta.get('c_gene') in ['TRBC1','TRBC2']:accepted.append(item)
            pair_records.append(dict(name=target['name'],sample=gsm,sample_label=tag,n_both_loop_candidate_barcodes=len(candidate),n_unambiguous_quality_paired_barcodes=len(accepted),n_rejected_barcodes=len(candidate)-len(accepted),n_alpha_TRAJ_TRAC_and_paired_TRB_role_supported_barcodes=len(role_supported),n_role_pairs_with_missing_beta_constant=sum(x['beta']['c_gene'] in ['None',''] for x in role_supported),role_support_examples=role_supported[:2],accepted_barcodes=[x['barcode'] for x in accepted],paired_contig_examples=accepted if target['name']=='NP105' else accepted[:2],source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),interpretation='same-cell paired provenance only; cells/barcoding replicates are not independent receptors, individuals, or studies; full_length is a source flag, not access to the full native V(D)J sequence'))
    for target in facts['rows']:
        supports=[x for x in pair_records if x['name']==target['name'] and x['n_alpha_TRAJ_TRAC_and_paired_TRB_role_supported_barcodes']>0]
        for role,kind in [('alpha','TRAV'),('beta','TRBV')]:
            name=target['name']+'_'+role;r=byname[name];full=raw_inputs[name]['sequence'];matched=[v for v in refs if v['gene'].startswith(kind) and full.startswith(v['framework_through_conserved_Cys'])]
            domain=r['numbered_boundary']['sequence'];jprefix='TRAJ' if role=='alpha' else 'TRBJ';scores=[dict(gene=j['gene'],function=j['function'],header=j['header'],n_exact_suffix_amino_acids=suffix_length(domain,j['sequence'])) for j in jrefs if j['gene'].startswith(jprefix)];maximum=max(x['n_exact_suffix_amino_acids'] for x in scores);best=[x for x in scores if x['n_exact_suffix_amino_acids']==maximum]
            gap=None
            if target['name'] in ['NP202','NP101','NP207','NP203'] and role=='beta':
                anchored=[v for v in matched if v['raw_slots_directly_anchored_for_beta']];assert anchored
                maps=[]
                for v in anchored:
                    i=0;m=[]
                    for position,aa in enumerate(v['gapped_amino_acids'][:104],1):
                        if aa!='-':m.append(dict(deposited_index0=i,IMGT_position=str(position),amino_acid=aa));i+=1
                    maps.append(m)
                assert all(m==maps[0] for m in maps)
                official={m['deposited_index0']:m for m in maps[0]};neural={m['deposited_index0']:m for m in r['methods']['ANARCII']['mapping']};classic={m['deposited_index0']:m for m in r['methods']['ANARCI']['mapping']}
                disputes=[i for i in neural if neural[i]!=classic[i]];assert disputes==[13,14,15]
                assert all(official[i]==classic[i] for i in official)
                assert all(neural[i]==classic[i] for i in neural if i not in disputes)
                gap=dict(disputed_deposited_indices0=disputes,canonical_gap_IMGT=14,ANARCII_gap_IMGT=17,source_supported_method='ANARCI',is_CDR3_insertion_conflict=False,reference_genes=sorted(v['gene'] for v in anchored),official_V_prefix_mapping=maps[0],disputed_residues=[dict(deposited_index0=i,amino_acid=official[i]['amino_acid'],ANARCII_position=neural[i]['IMGT_position'],ANARCI_position=classic[i]['IMGT_position'],reference_position=official[i]['IMGT_position']) for i in disputes],raw_outputs_modified=False)
            role_checks.append(dict(name=name,biological_chain_role=role,biological_role_primary_and_natural_pair_supported=bool(supports),role_basis='final primary construct labels; natural TRA/TRAJ/TRAC with same-cell TRB/TRBJ. Both-constant-complete pairs are counted separately; missing beta constant calls remain missing. V/DV locus or HMM ranking alone is not expressed chain identity',n_source_supported_barcodes=sum(x['n_alpha_TRAJ_TRAC_and_paired_TRB_role_supported_barcodes'] for x in supports),n_both_constant_genes_recorded_barcodes=sum(x['n_unambiguous_quality_paired_barcodes'] for x in supports),missing_beta_C_gene_not_imputed=True,exact_complete_framework_V_alleles=sorted(v['gene'] for v in matched),exact_complete_framework_V_genes=sorted({v['gene'].split('*')[0] for v in matched}),V_gene_uniquely_determined_by_framework_amino_acids=len({v['gene'].split('*')[0] for v in matched})==1,longest_exact_J_suffix_matches=best,J_rearrangement_cut_uniquely_inferred=False,FR1_gap_adjudication=gap,fixed_identity_admitted=False))
    np105=[r for r in pair_records if r['name']=='NP105'];result=dict(n_NP_constructs=6,n_chain_role_sources_supported=sum(r['biological_role_primary_and_natural_pair_supported'] for r in role_checks),n_FR1_gap_conflicts_source_adjudicated=sum(r['FR1_gap_adjudication'] is not None for r in role_checks),NP105_sample_results=np105,same_cell_pair_provenance=pair_records,chain_reference_adjudications=role_checks,n_new_fixed_identity_admissions=0,remaining='STCRDab full paired numbering absent; full native variable/constant sequence contigs not published in these annotation CSVs; framework-equivalent loci/alleles and historical V(D)J cut positions remain unresolved; TRBJ1-6 ORF does not pass unchanged functional-reference gate',protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    a.out.write_text(json.dumps(result,indent=2)+'\n');print('SOURCE_SUPPORTED_ROLES',result['n_chain_role_sources_supported'],'FR1_ADJUDICATED',result['n_FR1_gap_conflicts_source_adjudicated']);print('NP105',[(r['sample_label'],r['n_both_loop_candidate_barcodes'],r['n_unambiguous_quality_paired_barcodes']) for r in np105])
if __name__=='__main__':main()
