"""Source-checked classical SPR holdouts; no clinical or binary endpoint promotion."""
import hashlib, math, re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]


def paired_receptor_fasta(text, clone):
    """Select named TCR entities explicitly; an MHC alpha chain is not a TCR."""
    result={}
    for record in text.split('>')[1:]:
        header,*lines=record.splitlines()
        h=header.lower()
        if not re.search(r'(?<![a-z0-9])'+re.escape(clone.lower())+r'(?![a-z0-9])',h):continue
        if not ('receptor' in h or re.search(r'\btcr(?:\d+)?\b',h)):continue
        for chain in ('alpha','beta'):
            if chain+' chain' not in h:continue
            if chain in result:raise ValueError('ambiguous receptor chain')
            sequence=''.join(lines).strip()
            if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWY'):raise ValueError('noncanonical chain')
            result[chain]=sequence
    if set(result)!={'alpha','beta'}:raise ValueError('paired named receptor missing')
    return result


def validate_classical_additions(rows, *, check_sources=True):
    ids=set()
    for x in rows:
        if x['evidence_id'] in ids:raise ValueError('duplicate external observation')
        ids.add(x['evidence_id'])
        if (x['fold']!='external_test' or x.get('training_allowed') is not False or x.get('model_scored') is not False
            or x.get('clinical_risk_label') is not None or x.get('binary_binding_label') is not None
            or x.get('training_label') is not None):raise ValueError('external endpoint promotion')
        if x['endpoint']!='tcr_binding_spr' or x['species']!='Homo sapiens' or x['receptor_type']!='alphabeta':raise ValueError('wrong biological endpoint')
        if not x['peptide'] or set(x['peptide'])-set('ACDEFGHIKLMNPQRSTVWY'):raise ValueError('noncanonical peptide')
        if x['qualifier']!='=' or x['raw_unit'] not in {'uM','nM'}:raise ValueError('censored or unknown units')
        value=float(x['raw_value'])*({'uM':1000.,'nM':1.}[x['raw_unit']])
        if not math.isfinite(value) or value<=0 or not math.isclose(value,x['kd_nm'],rel_tol=1e-12):raise ValueError('KD conversion')
        if x.get('temperature_K') is None or not 273<=x['temperature_K']<=320:raise ValueError('assay temperature missing')
        if x['fit'] not in {'equilibrium','kinetic'} or 'Biacore' not in x['instrument']:raise ValueError('unverified SPR method')
        for chain in ('alpha','beta'):
            cdr=x['cdr3_'+chain];v=x[chain+'_variable']
            if len(cdr)<8 or cdr not in v or set(v)-set('ACDEFGHIKLMNPQRSTVWY'):raise ValueError('paired CDR identity')
        for key in ('lineage','study','research_group','conditions','primary_locator','sources','sequence_sources'):
            if not x.get(key):raise ValueError('source identity missing '+key)
        if check_sources:
            for source in x['sources']+x['sequence_sources']:
                p=(ROOT/source['path']).resolve()
                if not p.is_relative_to(ROOT) or hashlib.sha256(p.read_bytes()).hexdigest()!=source['sha256']:raise ValueError('source changed')
            sequence_source=ROOT/x['sequence_sources'][0]['path']
            chains=paired_receptor_fasta(sequence_source.read_text(),x['receptor_name'])
            for chain in ('alpha','beta'):
                if x[chain+'_variable'] not in chains[chain]:raise ValueError('deposited chain identity mismatch')
    return dict(measurements=len(rows),lineages=len({x['lineage'] for x in rows}),studies=len({x['study'] for x in rows}),research_groups=len({x['research_group'] for x in rows}))
