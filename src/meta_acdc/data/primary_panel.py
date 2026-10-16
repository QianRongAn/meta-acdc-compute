"""Validate primary measurements without promoting assay calls to training labels."""
from __future__ import annotations
import json, math, re
from pathlib import Path
from .build_release import domain
from .experimental_evidence import AA, reserved_hits, finite_nm
from .release import sha256


def validate_measurements(rows, root):
    seen=set(); identities={}; counts={}
    for row in rows:
        if row['evidence_id'] in seen: raise ValueError('duplicate evidence ID')
        seen.add(row['evidence_id'])
        if row.get('training_label') is not None: raise ValueError('external observations cannot become training labels')
        if row['species']!='Homo sapiens' or row['receptor_type']!='alphabeta': raise ValueError('panel receptor scope')
        if not re.fullmatch(r'HLA-[ABC]\*\d{2,3}:\d{2,3}', row['mhc_allele']): raise ValueError('ambiguous allele')
        if not row['peptide'] or not set(row['peptide'])<=AA: raise ValueError('not an individual canonical peptide')
        for item in row['sources']:
            path=Path(item['path']); path=path if path.is_absolute() else root/path
            if sha256(path)!=item['sha256']: raise ValueError('primary source changed')
        if not row['sources'] or not row['primary_locator'] or not row['conditions']: raise ValueError('missing provenance/conditions')
        seqs=[]
        for source in row['sequence_sources']:
            path=root/source['path']
            if sha256(path)!=source['sha256']: raise ValueError('sequence source changed')
            seqs.append(domain(json.loads(path.read_text())['entity_poly']['pdbx_seq_one_letter_code_can'])[0])
        if seqs!=[row['alpha_variable'],row['beta_variable']] or len(seqs)!=2: raise ValueError('paired construct mismatch')
        key=(row['pmid'],row['receptor_name'])
        pair=tuple(seqs)
        if key in identities and identities[key]!=pair: raise ValueError('construct identity changed')
        identities[key]=pair
        ep=row['endpoint']; outcome=row['observed_outcome'];kd=row['kd_nm']
        if ep=='tcr_binding_spr':
            if outcome=='measured_binding':
                if row['qualifier']!='=' or not isinstance(kd,(int,float)) or not math.isfinite(kd) or kd<=0: raise ValueError('invalid exact affinity')
                expected=finite_nm(row['raw_value'],row['raw_unit'])
                if expected is None or not math.isclose(expected,kd,rel_tol=1e-10):raise ValueError('raw value/unit disagrees with normalized KD')
            elif outcome=='not_detected':
                if kd is not None or row['raw_value']!='NB': raise ValueError('non-detection cannot be exact or infinite KD')
            else: raise ValueError('unrecognized binding outcome; censored values need a separate endpoint')
        elif ep=='tcr_cellular_tnf_alpha':
            if kd is not None or outcome not in {'activation_detected','activation_not_detected'}: raise ValueError('cellular response is not intrinsic binding')
        else: raise ValueError('unsupported endpoint')
        counts[ep]=counts.get(ep,0)+1
    return {'passed':True,'records':len(rows),'receptors':len(identities),'studies':len({x['pmid'] for x in rows}),
            'endpoints':counts,'training_allowed':False,'population_benchmark_allowed':False}


def reservation_exclusions(structures, graph_records, registry):
    """Remove entire reserved families, including graphs using a reserved donor."""
    families={x['family_id'] for x in structures if reserved_hits(x,registry)}
    complexes={x['complex_id'] for x in structures if x['family_id'] in families}
    excluded=[];kept=[]
    for i,row in enumerate(graph_records):
        reason='reserved_receptor_family' if row['family_id'] in families else 'reserved_graft_donor' if row.get('donor') in complexes else None
        if reason: excluded.append({'original_index':i,'graph_id':row['graph_id'],'reason':reason})
        else: kept.append(i)
    return kept,excluded,complexes
