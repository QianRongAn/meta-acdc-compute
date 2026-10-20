"""Bounded forensic replay of withdrawn outputs; never a scientific release.

Run in the existing trusted environment. Inputs stay read-only. Uses exactly the
archived 5c70164 source, not corrected current features. No training/test labels.
"""
from pathlib import Path
import csv, hashlib, io, json, statistics, subprocess, sys, tarfile, tempfile
ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = Path('/home/administrator/docs/meta-acdc')
COMMIT = '5c70164'
JOBS = ['5brz_EVDPIGHLY', '5brz_ESDPIVAQY', '5brz_ILAKFLHWL']
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    output = Path(sys.argv[1])
    if output.exists(): raise FileExistsError(output)
    commit = subprocess.check_output(['git','-C',str(ROOT),'rev-parse',COMMIT], text=True).strip()
    all_recovered = '--all-recovered' in sys.argv[2:]
    if all_recovered:
        with (ROOT/'docs/delivery_audit_20261018/recovered/da_scores_seed0_ensemble.tsv').open() as f:
            jobs = sorted(r['job_id'] for r in csv.DictReader(f,delimiter='\t'))
    else: jobs = JOBS
    inventory = json.loads((ROOT/'docs/delivery_audit_20261018/recovered_inventory.json').read_text())
    checkpoints = [ORIGINAL/'data/processed'/f'da_seed{s}.model.pt' for s in range(5)]
    for p in checkpoints:
        if sha(p) != inventory[p.name]['sha256']: raise ValueError('Historical checkpoint changed: '+p.name)
    with tempfile.TemporaryDirectory(prefix='meta-acdc-forensic-') as td:
        archive = subprocess.check_output(['git','-C',str(ROOT),'archive',commit,'src'])
        with tarfile.open(fileobj=io.BytesIO(archive)) as tf:
            for m in tf:
                if not m.isfile(): continue
                parts = Path(m.name).parts
                if parts[0] != 'src' or '..' in parts: raise ValueError(m.name)
                dest = Path(td)/m.name; dest.parent.mkdir(parents=True,exist_ok=True)
                dest.write_bytes(tf.extractfile(m).read())
        sys.path.insert(0,str(Path(td)/'src'))
        import torch
        from meta_acdc.models.egnn import EGNN
        from meta_acdc.structure.cif import parse_cif
        from meta_acdc.structure.graph import build_interface_graph_from_residues
        from meta_acdc.structure.score_predictions import score_graph
        torch.set_num_threads(2)
        nets=[]
        for p in checkpoints:
            c=torch.load(p,map_location='cpu',weights_only=False)
            net=EGNN(node_dim=c['n_node_dim'],edge_dim=c['n_edge_dim'],depth=6,hidden=128)
            net.load_state_dict(c['model']);net.eval();nets.append(net)
        results=[]; inputs={}
        for job in jobs:
            scores=[[] for _ in nets]
            paths = sorted((ORIGINAL/'data/raw/af3_predictions').glob(f'{job}_model_*.cif'))
            if not paths: raise ValueError('Missing inputs for '+job)
            for path in paths:
                inputs[str(path.relative_to(ORIGINAL))]=sha(path)
                graph=build_interface_graph_from_residues(parse_cif(path))
                for seed,net in enumerate(nets):
                    value=score_graph(net,graph,'cpu',mask_plddt=True,temperature=1.0)
                    scores[seed].append(value)
                pass
            print(job, len(paths), flush=True)
            for seed,vals in enumerate(scores):
                table=ROOT/'docs/delivery_audit_20261018/recovered'/f'da_scores_seed{seed}_ensemble.tsv'
                with table.open() as f: old={r['job_id']:r for r in csv.DictReader(f,delimiter='\t')}
                # Legacy scorer rounded each AF3 score before aggregating.
                mean=statistics.mean(float(f'{v:.4f}') for v in vals)
                results.append(dict(job_id=job,seed=seed,per_af3_model_scores=vals,
                    replay_mean=f'{mean:.4f}',historical_mean=old[job]['score'],
                    historical_n_models=int(old[job]['n_models']),replay_n_models=len(vals),
                    matches_legacy_four_decimals=f'{mean:.4f}'==old[job]['score'],
                    model_coverage_matches=len(vals)==int(old[job]['n_models'])))
        receipt=dict(record_status='FORENSIC_WITHDRAWN',clinical_interpretation_allowed=False,
            historical_causal_attestation=False,source_commit=commit,torch_version=torch.__version__,
            device='cpu',mask_plddt=True,temperature=1.0,selection=('all 188 recovered-table jobs; every retained AF3 model; all 5 retained seeds' if all_recovered else 'three prespecified disputed jobs; all retained AF3 models; all 5 retained seeds'),
            checkpoints={p.name:sha(p) for p in checkpoints},input_sha256=inputs,
            source_archive_sha256=hashlib.sha256(archive).hexdigest(),results=results,
            matches=sum(r['matches_legacy_four_decimals'] for r in results),comparisons=len(results))
        output.write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__': main()
