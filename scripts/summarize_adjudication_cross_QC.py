"""Assess current coordinates for the closest existing complete cross-study receptor identities."""
from pathlib import Path
import argparse,json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_adjudication_20261029'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    qc=json.loads((BASE/'fresh_cross_QC/native_atomic_qc.json').read_text());old=json.loads((ROOT/'docs/fingerprint_construct_followup_20261027/joint_readiness_final.json').read_text())['QC_pass_units'];rows=[]
    for pdb in ['4jfe','4jfh','6eqb','7n1f','7rtr']:
        entries=[r for r in qc['rows'] if r['pdb']==pdb];passed=[r for r in entries if r.get('status')=='native_graph_QC_pass'];hashes=[]
        for r in passed:
            previous=[o['native_graph_sha256'] for o in old if (o['pdb'],o['alpha_chain'],o['beta_chain'])==(r['pdb'],r['alpha_chain'],r['beta_chain'])];assert r['graph_sha256'] in previous;hashes.append(r['graph_sha256'])
        rows.append(dict(pdb=pdb,n_fixed_QC_pass=len(passed),accepted_graph_hashes_unchanged=bool(passed),graph_sha256=hashes,records=entries))
    by={r['pdb']:r for r in rows};pairs=[dict(anchor=x,partner=y,both_ends_fixed_QC_pass=bool(by[x]['n_fixed_QC_pass'] and by[y]['n_fixed_QC_pass'])) for x,y in [('7n1f','7rtr'),('6eqb','4jfe'),('6eqb','4jfh')]]
    out=dict(n_fresh_coordinate_PDBs=5,n_existing_complete_identity_cross_study_pairs_checked=len(pairs),n_current_QC_units_pass=sum(r['n_fixed_QC_pass'] for r in rows),n_new_qualified_cross_study_pairs=sum(r['both_ends_fixed_QC_pass'] for r in pairs),pairs=pairs,rows=rows,scope='focused current-coordinate revalidation of existing complete identities; two passed graphs reproduce existing development inventory and are not new data',protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    assert out['n_new_qualified_cross_study_pairs']==0 and out['n_current_QC_units_pass']==2
    a.out.write_text(json.dumps(out,indent=2)+'\n');print(out['pairs'])
if __name__=='__main__':main()
