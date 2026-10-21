"""Regenerate the manuscript's descriptive tables and figures from pinned public reports.

No training/scoring function is imported. The consumed test is represented only
by its already-public aggregate report; fresh reserve files are never opened.
"""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,csv,hashlib,json,math,platform,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def load_sources(root=ROOT):
    manifest=json.loads((HERE/'source_manifest.json').read_text())
    result={}
    for name,meta in manifest['files'].items():
        p=(root/name).resolve()
        if not p.is_relative_to(root.resolve()):raise ValueError('source escapes root')
        data=p.read_bytes()
        if hashlib.sha256(data).hexdigest()!=meta['sha256']:raise ValueError('source changed: '+name)
        result[name]=json.loads(data)
    return result

def load_context(root=ROOT):
    manifest=json.loads((HERE/'context_manifest.json').read_text())
    result={}
    for name,digest in manifest['files'].items():
        p=(root/name).resolve()
        if not p.is_relative_to(root.resolve()):raise ValueError('context escapes root')
        raw=p.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=digest:raise ValueError('context changed: '+name)
        result[name]=json.loads(raw)
    return result

def derive(s):
    rows=s['docs/n_review/structural_atom_audit.json']['rows']
    if len(rows)!=888:raise ValueError('historical audit coverage changed')
    parents={r['graph_id']:r for r in rows if r['kind']=='native'}
    if len(parents)!=300:raise ValueError('nonunique native IDs')
    atom=[];families=[];paired=[]
    for kind in ['native','graft','displaced']:
        part=[r for r in rows if r['kind']==kind]
        by_family=defaultdict(list)
        for r in part:
            c=r['construction']
            if c['severe_distance_A']!=1.2:raise ValueError('threshold changed')
            flag=int(c['cross_chain_pairs_below_cutoff']>0)
            by_family[r['family_id']].append(flag)
            if kind!='native':
                parent=r['graph_id'].split('<-',1)[0] if kind=='graft' else r['graph_id'].removesuffix('-disp')
                p=parents[parent]
                if p['family_id']!=r['family_id']:raise ValueError('parent family mismatch')
                paired.append(dict(control_id=r['graph_id'],native_id=parent,kind=kind,family_id=r['family_id'],native_overlap=int(p['construction']['cross_chain_pairs_below_cutoff']>0),control_overlap=flag))
        n=len(part);nflag=sum(sum(v) for v in by_family.values())
        family_rate=sum(sum(v)/len(v) for v in by_family.values())/len(by_family)
        atom.append(dict(kind=kind,n_records=n,n_audit_families=len(by_family),n_severe_overlap=nflag,severe_overlap_fraction=nflag/n,equal_family_overlap_fraction=family_rate,n_incomplete_interface=sum(bool(r['missing_nodes']) for r in part)))
        families.extend(dict(kind=kind,family_id=k,n_records=len(v),n_overlap=sum(v),overlap_fraction=sum(v)/len(v)) for k,v in sorted(by_family.items()))
    discriminator=[]
    for a in atom[1:]:
        p=[x for x in paired if x['kind']==a['kind']]
        specificity=1-sum(x['native_overlap'] for x in p)/len(p)
        discriminator.append(dict(control=a['kind'],n_pairs=len(p),sensitivity=a['severe_overlap_fraction'],specificity=specificity,balanced_accuracy=(a['severe_overlap_fraction']+specificity)/2,interpretation='in-sample retrospective finite-collection description; no biological or held-out performance'))
    report=s['docs/affinity_baseline_20261018/acceptance/report.json']
    if report['test_labels_used_for_selection'] or report['post_test_fit_allowed']:raise ValueError('unexpected test policy')
    affinity=[]
    for role,key in [('primary','primary'),('sensitivity','all15_sensitivity')]:
        v=report[key]
        for method,field in [('Ridge','model'),('Training mean','train_mean')]:
            affinity.append(dict(analysis=role,method=method,n=v['n'],n_sources=v['n_sources'],n_lineages=v['n_lineages'],**v[field]))
    for name,v in report['by_source'].items():
        for method,field in [('Ridge','model'),('Training mean','train_mean')]:
            affinity.append(dict(analysis=name,method=method,n=v['n'],n_sources=v['n_sources'],n_lineages=v['n_lineages'],**v[field]))
    compare=s['docs/affinity_development_20261018/bounded_run/comparison.json']
    if compare['test_used_for_selection'] or compare['consumed_test_evaluated'] or compare['new_test_evaluated']:raise ValueError('bounded development policy changed')
    validation=[dict(method=c['name'],alpha=c['alpha'],**c['validation'],prediction_sd=c['validation_prediction_sd']) for c in compare['candidates']]
    validation.append(dict(method='Training mean',alpha=None,**compare['train_mean_validation'],prediction_sd=0.))
    panel=s['docs/forensic_followup_20261021/pending_panel.json']
    actual_pairs={(r['sequences']['alpha'],r['sequences']['beta']) for r in panel['rows']}
    if len(actual_pairs)!=panel['counts']['exact_receptor_pairs']:raise ValueError('panel count mismatch')
    old=s['docs/forensic_followup_20261021/legacy_full_replay.json']
    per=s['docs/forensic_followup_20261021/per_structure_reconciliation.json']
    groups=[dict(label='FASTA files',count=len(panel['rows'])),dict(label='PDB prefixes',count=panel['counts']['pdb_prefixes']),dict(label='Exact receptor pairs',count=len(actual_pairs)),dict(label='Hamming components',count=panel['clustering_sensitivity']['equal_length_hamming']['n_components']),dict(label='Matching-block components',count=panel['clustering_sensitivity']['matching_blocks_max_length']['n_components'])]
    summary=dict(atom=atom,qc_only_discriminator=discriminator,panel=groups,affinity=affinity,validation=validation,
      rmse_relative_improvement=report['primary']['rmse_improvement']/report['primary']['train_mean']['rmse_pkd'],
      forensic=dict(job_seed_matches=old['matches'],job_seed_comparisons=old['comparisons'],structure_seed_matches=per['matches'],structure_seed_comparisons=per['comparisons']),
      scope='retrospective single-project case study; no new biological labels, model fits, or reserved-test scoring')
    return summary,families,paired

def write_csv(path,rows):
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def build(out,*,figures=True):
    out=Path(out)
    if out.exists():raise ValueError('refusing output overwrite; choose a fresh output directory')
    sources=load_sources();summary,families,paired=derive(sources)
    out.mkdir(parents=True)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    for name,rows in [('structural_qc',summary['atom']),('family_qc',families),('paired_qc',paired),('qc_only_discriminator',summary['qc_only_discriminator']),('input_units',summary['panel']),('affinity_metrics',summary['affinity']),('development_comparison',summary['validation'])]:write_csv(out/(name+'.csv'),rows)
    context=load_context()
    report=context['docs/six_blockers_20261021/ergo_comparison/report.json']
    if report['new_independent_acceptance'] or report['model_fit_or_tuning_performed']:raise ValueError('context endpoint/selection changed')
    rows=[dict(method=label,**report['metrics'][key],endpoint='MIRA activation association',new_independent_acceptance=False) for label,key in [('Official fixed ERGO-I','ergo_score'),('Historical SPR reference score','spr_reference_pkd_score')]]
    write_csv(out/'functional_context.csv',rows)
    if figures:plot(summary,sources,out)
    outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file()}
    receipt=dict(analysis='descriptive report assembly',source_manifest_sha256=hashlib.sha256((HERE/'source_manifest.json').read_bytes()).hexdigest(),context_manifest_sha256=hashlib.sha256((HERE/'context_manifest.json').read_bytes()).hexdigest(),builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=sys.version,platform=platform.platform(),outputs=outputs,training_performed=False,new_test_scoring=False,individual_test_predictions_loaded=False,fresh_reserve_loaded=False)
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return summary

def plot(d,s,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    blue='#26547c';orange='#c86532';grey='#697786'
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=190,bbox_inches='tight')
        fig.savefig(out/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    labels=['Native','Grafted','Displaced'];a=d['atom']
    for j,(key,title) in enumerate([('severe_overlap_fraction','Severe cross-chain overlap'),('equal_family_overlap_fraction','Equal weight per audit family')]):
        vals=[r[key]*100 for r in a];axs[j].bar(labels,vals,color=[blue,orange,grey]);axs[j].set_ylim(0,113);axs[j].set_ylabel('Flagged records (%)' if j==0 else 'Mean within-family flag rate (%)');axs[j].set_title(title)
        for i,v in enumerate(vals):axs[j].text(i,v+2,f'{v:.1f}%',ha='center')
    fig.suptitle('Figure 1. Fixed 1.2-A screen on historical structures',fontsize=12)
    save(fig,'figure1_structural_qc')
    fig,ax=plt.subplots(figsize=(8.8,4.2),layout='constrained');rows=d['panel'];labels=['FASTA files','PDB prefixes','Exact receptor\npairs','98% Hamming\ncomponents','98% matching-block\ncomponents'];vals=[r['count'] for r in rows]
    ax.bar(labels,vals,color=[grey,grey,blue,orange,orange]);ax.set_ylabel('Count');ax.set_ylim(0,390)
    for i,v in enumerate(vals):ax.text(i,v+6,str(v),ha='center')
    ax.set_title('Figure 2. Input counts depend on the unit and similarity rule');save(fig,'figure2_input_units')
    fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
    pri=[r for r in d['affinity'] if r['analysis']=='primary']
    for j,key in enumerate(['rmse_pkd','mae_pkd']):
        vals=[r[key] for r in pri];axs[j].bar(['Ridge','Training mean'],vals,color=[blue,grey]);axs[j].set_ylim(0,1.2);axs[j].set_ylabel(('RMSE' if j==0 else 'MAE')+' (pKD)')
        for i,v in enumerate(vals):axs[j].text(i,v+.025,f'{v:.4f}',ha='center')
        axs[j].set_title('Already-consumed frozen test; 14 observations')
    fig.suptitle('Figure 3. Minimal improvement over a constant comparator',fontsize=12);save(fig,'figure3_affinity')
    fig,axs=plt.subplots(1,2,figsize=(10,4.3),layout='constrained')
    rows=d['validation'];labels=['Original\nreference','Sequence\nbalanced','CDR3\nchemistry','CDR3\nwords','Training\nmean']
    for j,(key,ylabel) in enumerate([('rmse_pkd','Validation RMSE (pKD)'),('prediction_sd','Validation prediction SD (pKD)')]):
        vals=[r[key] for r in rows];axs[j].bar(range(len(rows)),vals,color=[blue,orange,orange,orange,grey]);axs[j].set_xticks(range(len(rows)),labels);axs[j].set_ylabel(ylabel);axs[j].set_ylim(0,max(vals)*1.25)
        for i,v in enumerate(vals):axs[j].text(i,v+max(vals)*.025,f'{v:.3f}',ha='center',fontsize=9)
    fig.suptitle('Figure 4. Bounded development: output spread does not ensure lower error',fontsize=12);save(fig,'figure4_development')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--no-figures',action='store_true');a=p.parse_args();x=build(a.out,figures=not a.no_figures);print(json.dumps({'output':str(a.out),'graphs':sum(r['n_records'] for r in x['atom']),'new_test_scoring':False}))
