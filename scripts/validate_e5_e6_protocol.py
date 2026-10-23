"""Validate prospective synthetic plans without generating study outcomes."""
from pathlib import Path
import argparse,hashlib,json,math

ROOT=Path(__file__).resolve().parents[1]

def stream_seed(experiment,version,scenario_id,replicate,stream):
    key='|'.join(map(str,[experiment,version,scenario_id,replicate,stream]))
    return int(hashlib.sha256(key.encode()).hexdigest()[:16],16)

def validate(root):
    d=root/'docs/e5_e6_protocol_20261023'
    a=json.loads((d/'E5_plan.json').read_text());b=json.loads((d/'E6_plan.json').read_text())
    for plan in [a,b]:
        assert plan['data_mode']=='synthetic_only'
        assert not any(plan[k] for k in ['biological_data_read','heldout_labels_read','heldout_scores_read','study_runs_executed'])
        assert plan['repetitions']==5000 and plan['version']=='1.0'
        assert plan['planned_git_tag']=='protocol/e5-e6-v1-20261023'
    primary=a['primary'];grid=a['grid']
    assert primary['n_selection'] in grid['n_selection'] and primary['k'] in grid['k'] and primary['rho'] in grid['rho']
    assert all(n>=4 and n%2==0 for n in grid['n_selection'])
    assert grid['k']==sorted(set(grid['k'])) and grid['k'][0]==1 and grid['k'][-1]==50
    assert all(0<=v<=1 for v in grid['rho']) and 1.0 in grid['rho']
    assert primary['n_shadow_positive']*2==primary['n_shadow']
    assert b['application_gate']['actual_minimum_benefit_delta'] is None
    assert b['application_gate']['status']=='not_authorized_by_this_simulation'
    assert b['bootstrap_draws']==1999
    assert len({s['id'] for s in b['scenarios']})==len(b['scenarios'])==13
    assert b['primary']['scenario'] in {s['id'] for s in b['scenarios']}
    count=0
    for s in b['scenarios']:
        assert 2<=s['blocks']<=s['n']
        rho=s['prediction_error_source_fraction'];technical=s['technical_SD'];gamma=s['paired_latent_error_correlation']
        assert 0<=rho<1 and 0<=technical<1 and 0<=gamma<=1
        assert math.isclose(rho+(1-rho)/2+(1-rho)/2,1)
        for delta in b['effect_grid_pKD_RMSE']:
            c2=((1-delta)**2-technical**2)/(1-technical**2)
            assert 0<c2<=1
            assert math.isclose(math.sqrt((1-technical**2)+technical**2),1)
            assert math.isclose(math.sqrt(c2*(1-technical**2)+technical**2),1-delta)
            count+=1
    assert 1/(2**3)>b['primary_false_pass_alpha']
    assert stream_seed('E5','1.0','n14_rho0.5',0,'labels')!=stream_seed('E5','1.0','n14_rho0.5',0,'shadow')
    return dict(status='protocol_parameter_and_algebra_checks_passed',E6_scenarios=13,E6_effect_cells_checked=count,
                E5_grid_cells=len(grid['n_selection'])*len(grid['k'])*len(grid['rho']),
                study_runs_executed=False,biological_data_read=False,
                note='No simulated performance, power, CI coverage or model outcome was calculated.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);a=parser.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    result=validate(ROOT);a.out.write_bytes((json.dumps(result,indent=2)+'\n').encode());print(json.dumps(result))
