"""Published engineered DNA is not a clinical lot; pending SPR is not a core test."""
from .clinical_primary import CODONS
from .evidence_diversity import _sources
from .spr_followup import paired_receptor_fasta,ROOT


def validate_1d3_reference(x,*,check_sources=True):
    dna=x['vector_DNA']
    if len(dna)!=6338 or set(dna)-set('ACGT') or dna!=''.join(t['nt'] for t in x['rows']):raise ValueError('published vector transcription')
    for t in x['rows']:
        if t['pdf_page'] not in {130,131,132} or ''.join(t['text'].split())!=t['nt']:raise ValueError('source row transcription')
    coordinates={'beta':(1103,2023),'linker':(2024,2089),'alpha':(2090,2899),'insert':(1103,2899)}
    for name,(start,end) in coordinates.items():
        p=x['parts'][name];nt=dna[start-1:end]
        protein=''.join(CODONS[nt[i:i+3]] for i in range(0,len(nt),3))
        if (p['start'],p['end'],p['nt'],p['protein'])!=(start,end,nt,protein):raise ValueError('published coordinates or translation')
        if p['length_nt']!=len(nt) or p['length_aa']!=len(protein):raise ValueError('part length')
    a=x['parts']['alpha']['protein'];b=x['parts']['beta']['protein'];link=x['parts']['linker']['protein'];full=x['parts']['insert']['protein']
    if not a.endswith('*') or '*' in a[:-1] or '*' in b or '*' in link:raise ValueError('ORF stop')
    if a[179]!='C' or b[190]!='C' or x['alpha_cysteine_position']!=180 or x['beta_cysteine_position']!=191:raise ValueError('engineered cysteines')
    if x['cdr3_alpha']!='CALGGGADGLTF' or x['cdr3_beta']!='CASSLSGTTADEQYF' or x['cdr3_alpha'] not in a or x['cdr3_beta'] not in b:raise ValueError('1D3 paired identity')
    if link!='GSGATNFSLLKQAGDVEENPGP' or x['chain_order']!='beta-P2A-alpha' or full!=b+link+a:raise ValueError('chain order or linker')
    if x['full_insert_DNA']!=x['parts']['insert']['nt'] or x['reference_protein']!=full[:-1] or len(x['reference_protein'])!=598:raise ValueError('coding insert')
    if (x['clinical_exact_insert_verified'] or x['clinical_lot_verified'] or x['clinical_vector_sequence_verified'] or x['training_allowed'] or x['model_scored'] or x['clinical_risk_label'] is not None or x['fold']!='clinical_reference_reserved'):raise ValueError('reference promoted to clinical proof')
    if check_sources:_sources(x)
    return True


def validate_pending_pap(rows,*,check_sources=True):
    # These are the equilibrium column of primary FigS6, not the neighbouring kinetic estimates.
    expected={'wt':30.,'S32Halpha':61.,'S32Malpha':5.3,'S32Qalpha':9.8,'S30E32Qalpha':11.5}
    if len(rows)!=5 or {x['variant'] for x in rows}!=set(expected):raise ValueError('PAP coverage')
    for x in rows:
        if (x['fold']!='external_pending' or x['scope']!='reserved_format_uncertain_SPR_not_classical_core' or x['classical_core_eligible'] is not False or x['training_allowed'] is not False or x['model_scored'] is not False or x['clinical_risk_label'] is not None or x['binary_binding_label'] is not None):raise ValueError('pending SPR promotion')
        if (x['lineage']!='TCR156' or x['study']!='10.1126/science.adx3162' or x['species']!='Homo sapiens' or x['peptide']!='TLMSAMTNL' or x['mhc_allele']!='HLA-A*02:01' or x['endpoint']!='tcr_binding_spr'):raise ValueError('PAP source identity')
        if x['kd_um']!=expected[x['variant']] or x['raw_unit']!='uM' or x['qualifier']!='=' or x['fit']!='equilibrium':raise ValueError('wrong primary KD column or units')
        if x['temperature_K']!=298.15 or x['instrument']!='Biacore T100' or 'tetramer' not in x['immobilized_ligand_reported']:raise ValueError('source method altered')
        if x['buffer'] is not None or x['flow_rate_ul_min'] is not None or x['replication_count'] is not None or x['primary_error_um'] is not None:raise ValueError('unreported condition imputed')
        for c in ('alpha','beta'):
            if x['cdr3_'+c] not in x[c+'_variable']:raise ValueError('paired identity')
        if check_sources:
            _sources(x);_sources({'sources':x['sequence_sources']})
            seq=paired_receptor_fasta((ROOT/x['sequence_sources'][0]['path']).read_text(),'TCR156')
            for c,delimiter in [('alpha','NIQNP'),('beta','EDLKN')]:
                if seq[c].split(delimiter)[0]!=x[c+'_variable']:raise ValueError('deposited variable chain changed')
    return dict(measurements=5,lineages=1,studies=1,eligible_classical_core=0)
