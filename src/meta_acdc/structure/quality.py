"""Per-model QC and identity attestation for predicted class-I complexes."""
from __future__ import annotations
import csv,json,math,re
from collections import defaultdict
from pathlib import Path
from .identity import chain_sequence,validate_peptide_atoms
from .cif import parse_cif
from ..data.build_release import domain
from ..data.release import sha256


def _same_variable(observed,predicted):
    a=domain(observed)[0];b=domain(predicted)[0]
    return min(len(a),len(b))>=70 and (a in b or b in a) and min(len(a),len(b))/max(len(a),len(b))>=.95


def annotate_prediction(residues, reference, expected_peptide):
    chains=defaultdict(list)
    for r in residues:chains[r.chain].append(r)
    seqs={c:chain_sequence(rs) for c,rs in chains.items()}
    used=set();roles={};kinds={}
    for kind,ref in [('tra',reference['alpha']),('trb',reference['beta']),('mhc_i',reference['mhc_sequence']),('b2m',reference['b2m_sequence']),('peptide',expected_peptide)]:
        hits=[]
        for c,s in seqs.items():
            if c in used:continue
            ok=_same_variable(ref,s) if kind in {'tra','trb'} else s==ref if kind=='peptide' else (ref in s or s in ref) and min(len(ref),len(s))/max(len(ref),len(s))>=.95
            if ok:hits.append(c)
        if len(hits)!=1:raise ValueError(f'identity {kind}: expected exactly one sequence match, got {len(hits)}')
        c=hits[0];used.add(c);kinds[c]=kind;roles[c]='tcr' if kind in {'tra','trb'} else 'mhc' if kind=='mhc_i' else kind
    if set(chains)!=used:raise ValueError('unexpected extra chains in prediction')
    for r in residues:r.role=roles[r.chain];r.chain_kind=kinds[r.chain]
    return residues


def confidence(path):
    stem=path.stem;job,idx=stem.rsplit('_model_',1)
    candidates=[path.with_name(f'{job}_conf_{idx}.json'),path.with_name(f'{job}_summary_confidences_{idx}.json')]
    if idx=='0':candidates.append(path.with_name(job+'_conf.json'))
    for p in candidates:
        if not p.exists():continue
        try:
            d=json.loads(p.read_text());v=float(d.get('iptm',d.get('tcr-pmhc_ipTM')))
            if math.isfinite(v) and 0<=v<=1:return v,str(p),sha256(p)
        except (ValueError,TypeError,json.JSONDecodeError):pass
    return None,'',''


def qc_file(path, references):
    job=path.stem.rsplit('_model_',1)[0];pdb=job.split('_',1)[0].lower();expected=job.split('_',1)[1] if '_' in job else ''
    row={'instance':path.parent.name,'job_id':job,'model_id':path.stem,'source_path':str(path.resolve()),'file_sha256':sha256(path),'iptm':'','confidence_path':'','confidence_sha256':'','reference_id':'','reference_sha256':'','family_id':'','native':'no','flag':''}
    flags=[];v,cp,ch=confidence(path);row.update(iptm='' if v is None else v,confidence_path=cp,confidence_sha256=ch)
    if v is None:flags.append('no_model_confidence')
    elif v<.75:flags.append('low_iptm')
    try:
        rs=parse_cif(path)
        if not rs:raise ValueError('empty structure')
        accepted=[]
        for ref in references.get(pdb,[]):
            try:annotate_prediction(rs,ref,expected);accepted.append(ref)
            except ValueError:continue
        identities={(x['alpha'],x['beta'],x['mhc_sequence'],x['b2m_sequence']) for x in accepted}
        if len(identities)!=1:raise ValueError('no unique reference receptor/pMHC sequence match')
        ref=sorted(accepted,key=lambda x:x['complex_id'])[0];annotate_prediction(rs,ref,expected)
        validate_peptide_atoms(rs)
        row.update(reference_id=ref['complex_id'],reference_sha256=ref['source_sha256'],family_id=ref['family_id'],native='yes' if expected==ref['peptide'] else 'no')
        from .graph import build_interface_graph_from_residues
        g=build_interface_graph_from_residues(rs)
        if g.n_nodes<20:raise ValueError('small interface')
    except Exception as exc:
        flags.append('structure_or_identity_failed');row['error']=f'{type(exc).__name__}: {exc}'
    row['flag']=';'.join(flags);return row


def build_qc(roots,references,output):
    by=defaultdict(list)
    for record in references:
        if record.get('domain','crystal')=='crystal':by[record['pdb']].append(record)
    rows=[]
    for root in roots:
        for p in sorted(root.rglob('*_model_*.cif')):
            rows.append(qc_file(p,by))
            if len(rows)%250==0:print(f'QC files={len(rows)} passed={sum(not x["flag"] for x in rows)}',flush=True)
    fields=['instance','job_id','model_id','source_path','file_sha256','iptm','confidence_path','confidence_sha256','reference_id','reference_sha256','family_id','native','flag','error']
    with output.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    return rows
