"""Check all severe-overlap counts directly from the public coordinate archive.

Uses a Python standard-library cell grid instead of the original KD-tree. No raw PDB,
AF3 download, torch checkpoint, or training dependency is required.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,math
from itertools import product
from collections import defaultdict
HERE=Path(__file__).resolve().parent

def spatial_grid(points,cutoff=1.2):
    if not math.isfinite(cutoff) or cutoff<=0:raise ValueError('invalid cutoff')
    grid=defaultdict(list)
    for point in points:
        if len(point)!=3 or not all(math.isfinite(v) for v in point):raise ValueError('invalid coordinates')
        grid[tuple(math.floor(v/cutoff) for v in point)].append(point)
    return grid

def count_close_pairs(peptide,grid,cutoff=1.2):
    count=0
    for p in peptide:
        if len(p)!=3 or not all(math.isfinite(v) for v in p):raise ValueError('invalid coordinates')
        cell=tuple(math.floor(v/cutoff) for v in p)
        for offset in product((-1,0,1),repeat=3):
            for q in grid.get(tuple(a+b for a,b in zip(cell,offset)),()):
                count+=sum((a-b)**2 for a,b in zip(p,q))<=cutoff**2
    return count

def verify(folder):
    folder=Path(folder);receipt=json.loads((folder/'export_receipt.json').read_text())
    data=(folder/'atomic_inputs.json.gz').read_bytes()
    if hashlib.sha256(data).hexdigest()!=receipt['archive_sha256']:raise ValueError('coordinate archive changed')
    a=json.loads(gzip.decompress(data));expected=json.loads((HERE.parents[1]/'docs/n_review/structural_atom_audit.json').read_text());old={r['graph_id']:r for r in expected['rows']}
    seen=set();counts={};bad=[];grids={}
    for row in a['graphs']:
        if row['graph_id'] in seen:raise ValueError('duplicate graph')
        seen.add(row['graph_id']);env=[at[1:] for r in a['native_residues'][row['complex_id']] if r['role']!='peptide' for at in r['atoms'] if not at[0].upper().startswith('H')]
        pep=[at[1:] for r in row['peptide_residues'] for at in r['atoms'] if not at[0].upper().startswith('H')]
        if not env or not pep:raise ValueError('empty atomic partition')
        key=row['complex_id']
        if key not in grids:grids[key]=spatial_grid(env)
        n=count_close_pairs(pep,grids[key])
        if n!=old[row['graph_id']]['construction']['cross_chain_pairs_below_cutoff']:bad.append(row['graph_id'])
        stats=counts.setdefault(row['kind'],{'n':0,'n_overlap':0});stats['n']+=1;stats['n_overlap']+=n>0
    if seen!=set(old):raise ValueError('audit coverage differs')
    if bad:raise ValueError('independent distance mismatch '+str(bad[:5]))
    return dict(all_counts_match=True,n_graphs=len(seen),counts=counts,method='Python standard-library spatial cell grid; <=1.2 A; heavy-atom peptide versus environment',coordinates_sha256=receipt['archive_sha256'],source_audit_sha256=hashlib.sha256((HERE.parents[1]/'docs/n_review/structural_atom_audit.json').read_bytes()).hexdigest(),training_or_scoring=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise ValueError('refusing verification overwrite')
    r=verify(a.bundle);a.out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
