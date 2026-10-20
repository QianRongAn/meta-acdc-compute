"""Distinct receptor evidence levels: complete protein reference versus paired-CDR3 SPR."""
import hashlib, math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
AA3=dict(zip('Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val'.split(),'ARNDCQEGHILKMFPSTWYV'))

def _sources(row):
    for source in row['sources']:
        p=(ROOT/source['path']).resolve()
        if not p.is_relative_to(ROOT) or hashlib.sha256(p.read_bytes()).hexdigest()!=source['sha256']:raise ValueError('source changed')

def validate_kras_reference(x,*,check_sources=True):
    for chain,n in [('alpha',269),('beta',303)]:
        try:sequence=''.join(AA3[t] for t in x[chain+'_rows'].split())
        except KeyError as e:raise ValueError('unrecognized amino acid transcription') from e
        if sequence!=x[chain] or len(sequence)!=n or x['cdr3_'+chain] not in sequence:raise ValueError('reference chain transcription')
    if x['linker']!='RAKRSGSGATNFSLLKQAGDVEENPGP' or x['chain_order']!='alpha-SEQ58-beta':raise ValueError('published linker/order')
    if x['reference_protein']!=x['alpha']+x['linker']+x['beta'] or len(x['reference_protein'])!=599:raise ValueError('reference assembly')
    if (x['full_DNA'] is not None or x['clinical_lot_verified'] or x['clinical_exact_insert_verified'] or x['complete_vector_sequence'] or x['training_allowed'] or x['clinical_risk_label'] is not None):raise ValueError('reference promoted to clinical DNA proof')
    if check_sources:_sources(x)
    return True

def validate_partial_spr(rows,*,check_sources=True):
    ids=set();reserved=exposed=0
    for x in rows:
        if x['evidence_id'] in ids:raise ValueError('duplicate measurement')
        ids.add(x['evidence_id'])
        if x['fold']!='external_pending' or x['training_allowed'] or x['model_scored'] or x['binary_binding_label'] is not None or x['clinical_risk_label'] is not None:raise ValueError('pending evidence promotion')
        if x['alpha_variable'] is not None or x['beta_variable'] is not None or x['full_variable_sequences_verified']:raise ValueError('full sequence imputation')
        if x['temperature_K'] is None or not 273<=x['temperature_K']<=320:raise ValueError('missing assay temperature')
        if not math.isfinite(x['kd_um']) or x['kd_um']<=0 or x['error_um']<0:raise ValueError('invalid SPR value')
        for chain in ('alpha','beta'):
            cdr=x['cdr3_'+chain]
            if len(cdr)<8 or set(cdr)-set('ACDEFGHIKLMNPQRSTVWY'):raise ValueError('paired CDR missing')
        if not x.get('trav') or not x.get('trbv'):raise ValueError('V genes missing')
        hits=any(x['overlap_checks'].values())
        if hits:
            if x['scope']!='historically_exposed_not_external_test':raise ValueError('historical exposure hidden')
            exposed+=1
        else:
            if x['scope']!='reserved_partial_identity_not_model_ready':raise ValueError('partial identity hidden')
            reserved+=1
        c=x['counterligand']
        if c['kd_um'] is not None or c['binary_binding_label'] is not None or c['max_tested_um'] is not None:raise ValueError('NB promoted to exact negative')
        if check_sources:_sources(x)
    return dict(reserved_partial_identity=reserved,historically_exposed=exposed,model_ready=0)
