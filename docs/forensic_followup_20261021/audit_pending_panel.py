"""Inventory the original tracked pending panel without modifying sequences."""
from pathlib import Path
import csv,difflib,hashlib,io,json,subprocess,tarfile,sys
ROOT=Path(__file__).resolve().parents[2]
ROLES={'alpha','beta','peptide','mhc','b2m'}
def sha(data): return hashlib.sha256(data).hexdigest()
def parse(text):
    records={}; key=None
    for line in text.splitlines():
        line=line.strip()
        if line.startswith('>'):
            key=line[1:]
            if key in records: raise ValueError('Duplicate role '+key)
            records[key]=''
        elif line:
            if key is None: raise ValueError('Sequence before header')
            records[key]+=line
    if set(records)!=ROLES: raise ValueError('Missing or unknown roles')
    if any(not s or set(s)-set('ACDEFGHIKLMNPQRSTVWY') for s in records.values()):
        raise ValueError('Empty or noncanonical sequence')
    return records

def components(pairs,method):
    parent=list(range(len(pairs)));edges=[]
    def find(i):
        while parent[i]!=i: i=parent[i]
        return i
    for i,a in enumerate(pairs):
        for j,b in enumerate(pairs[:i]):
            if method=='equal_length_hamming':
                good=all(len(x)==len(y) and sum(c==d for c,d in zip(x,y))/len(x)>=.98 for x,y in zip(a,b))
            else:
                good=all(sum(m.size for m in difflib.SequenceMatcher(None,x,y,autojunk=False).get_matching_blocks())/max(len(x),len(y))>=.98 for x,y in zip(a,b))
            if good:
                edges.append([i,j]);parent[find(i)]=find(j)
    groups={}
    for i in range(len(pairs)):groups.setdefault(find(i),[]).append(i)
    return dict(n_components=len(groups),members=list(groups.values()),edges=edges)

def audit():
    files=sorted((ROOT/'af3_pending_upload').rglob('*.fasta'));rows=[]
    origin=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','4f730d2'],text=True).strip()
    raw=subprocess.check_output(['git','-C',str(ROOT),'archive',origin,'af3_pending_upload'])
    with tarfile.open(fileobj=io.BytesIO(raw)) as tf:
        originals={m.name:sha(tf.extractfile(m).read()) for m in tf if m.isfile() and m.name.endswith('.fasta')}
    for path in files:
        data=path.read_bytes();seq=parse(data.decode());rel=path.relative_to(ROOT).as_posix()
        rows.append(dict(path=rel,job_id=path.stem,pdb=path.stem.split('_')[0],sha256=sha(data),
            matches_first_commit=originals.get(rel)==sha(data),sequences=seq,
            lengths={k:len(v) for k,v in seq.items()}))
    pairs=sorted({(r['sequences']['alpha'],r['sequences']['beta']) for r in rows})
    with (ROOT/'af3_pending_upload/pending_upload_manifest.tsv').open() as f:
        manifest=list(csv.DictReader(f,delimiter='\t'))
    manifest_jobs=[r['job_id'] for r in manifest];jobs=[r['job_id'] for r in rows]
    comparisons={}
    for base in [Path('/home/administrator/docs/meta-acdc'),Path('/home/administrator/autodl-tmp/meta-acdc-compute')]:
        found={p.relative_to(base).as_posix():sha(p.read_bytes()) for p in (base/'af3_pending_upload').rglob('*.fasta')}
        comparisons[str(base)]=dict(n_files=len(found),exact_inventory_match=found=={r['path']:r['sha256'] for r in rows})
    return dict(record_status='FORENSIC_WITHDRAWN_INPUT_PANEL',clinical_interpretation_allowed=False,
        first_commit=origin,first_commit_fasta_count=len(originals),rows=rows,
        counts=dict(files=len(rows),pdb_prefixes=len({r['pdb'] for r in rows}),
            exact_receptor_pairs=len(pairs),unique_sequences={k:len({r['sequences'][k] for r in rows}) for k in sorted(ROLES)},
            beta_starts_NA=sum(r['sequences']['beta'].startswith('NA') for r in rows),
            mhc_starts_GS=sum(r['sequences']['mhc'].startswith('GS') for r in rows)),
        manifest=dict(n_rows=len(manifest),unique_jobs=len(set(manifest_jobs)),
            matches_fasta_job_set=set(manifest_jobs)==set(jobs),duplicate_fasta_jobs=len(jobs)-len(set(jobs))),
        all_bytes_match_first_commit=all(r['matches_first_commit'] for r in rows),copy_comparisons=comparisons,
        pair_index=[dict(alpha=a,beta=b) for a,b in pairs],
        clustering_sensitivity={m:components(pairs,m) for m in ['equal_length_hamming','matching_blocks_max_length']},
        clustering_limit='Both chains separately >=98%, connected components; full construct sequences, not validated V-gene families. Original report algorithm unavailable; 73 not assumed.',
        biological_independence_established=False)
if __name__=='__main__':
    output=Path(sys.argv[1])
    if output.exists():raise FileExistsError(output)
    output.write_text(json.dumps(audit(),indent=2)+'\n')
