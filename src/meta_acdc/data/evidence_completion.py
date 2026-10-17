"""Versioned evidence completion: raw values, clinical references and external panels stay distinct."""
from pathlib import Path
import json,math,statistics
from .clinical_primary import CODONS
ROOT=Path(__file__).resolve().parents[3]
CURATION=ROOT/'docs/spr_followup_20261018'

def validate_c4(row):
    dna=row['dna'];protein=row['amino_acids']
    if len(dna)!=1809 or set(dna)-set('ACGT') or len(dna)%3:raise ValueError('C4 coding sequence length/alphabet')
    translated=''.join(CODONS[dna[i:i+3]] for i in range(0,len(dna),3))
    if translated!=protein+'*' or len(protein)!=602:raise ValueError('C4 translation mismatch')
    if row['beta_chain']+row['linker']+row['alpha_chain']!=protein or row['chain_order']!='beta-P2A-alpha':raise ValueError('C4 chain order')
    if row['beta_chain'][187]!='C' or row['alpha_chain'][176]!='C':raise ValueError('C4 cysteine positions')
    if row['cdr3_alpha'] not in row['alpha_chain'] or row['cdr3_beta'] not in row['beta_chain']:raise ValueError('C4 CDR identity')
    if row['clinical_lot_sequence_certificate'] or row['complete_vector_sequence'] or row['clinical_risk_label'] is not None or row['training_allowed']:raise ValueError('reference sequence promoted to clinical proof')
    return True

def validate_external(rows):
    ids=set()
    for x in rows:
        if x['evidence_id'] in ids:raise ValueError('duplicate external observation')
        ids.add(x['evidence_id'])
        if x['fold']!='external_test' or x['training_allowed'] or x['model_scored'] or x['clinical_risk_label'] is not None or x['binary_binding_label'] is not None:raise ValueError('external holdout/endpoint promotion')
        if x['scope']!='external_format_shift_challenge_only' or x['measurement_format']!='soluble TCR-antiCD3 bispecific ImmTAV':raise ValueError('external format hidden')
        if x['raw_unit']!='M':raise ValueError('external unit mismatch')
        values=x['raw_kd_replicates']
        if len(values)!=3 or any(not isinstance(v,(float,int)) or not math.isfinite(v) or v<=0 for v in values):raise ValueError('invalid external replicate')
        if not math.isclose(x['kd_um'],statistics.mean(values)*1e6,rel_tol=1e-12):raise ValueError('external KD conversion')
        if not x['alpha_variable'] or not x['beta_variable'] or not x['cdr3_alpha'] or not x['cdr3_beta']:raise ValueError('paired external identity missing')
    return dict(measurements=len(rows),lineages=len({x['lineage'] for x in rows}),studies=len({x['study'] for x in rows}),biological_receptors=len({x['lineage'] for x in rows}))

def structure_registry(records):
    """Historical structure records use alpha/beta, unlike reservation alpha_variable/beta_variable."""
    result=[]
    for x in records:
        alpha=x.get('alpha_variable') or x.get('alpha');beta=x.get('beta_variable') or x.get('beta')
        if not alpha or not beta:raise ValueError('historical paired identity missing')
        result.append(dict(receptor_name=x.get('complex_id') or x.get('receptor_name') or x.get('evidence_id'),alpha_variable=alpha,beta_variable=beta,cdr3_alpha=x.get('cdr3_alpha',''),cdr3_beta=x.get('cdr3_beta','')))
    return result
