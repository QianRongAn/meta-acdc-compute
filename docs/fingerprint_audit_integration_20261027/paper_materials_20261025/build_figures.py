"""Build audit figures and family-level source tables from frozen artifacts."""
import csv
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import torch
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
WORK = HERE.parent
ROOT = Path('/home/administrator/docs/meta-acdc-fingerprint-review')
sys.path.insert(0, str(WORK / 'audit_followup_20261025'))
from verify_metrics import summarize

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'svg.fonttype': 'none', 'savefig.facecolor': 'white'})
BLUE, ORANGE, GRAY, GREEN = '#245B83', '#B65D2D', '#65717C', '#48745C'
INPUTS = {}
ROWS = []


def read(path):
    INPUTS[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return json.loads(path.read_text(encoding='utf-8-sig'))


def point(ax, estimate, y, color, name):
    mean, (lo, hi) = estimate['mean'], estimate['ci95']
    assert np.isfinite([mean, lo, hi]).all() and lo <= mean <= hi
    ax.errorbar(mean, y, xerr=[[mean-lo], [hi-mean]], fmt='o', color=color,
                markersize=6, capsize=3, linewidth=1.5)
    ROWS.append(dict(figure=FIGURE, panel=name, mean=mean, ci_low=lo, ci_high=hi))


def save(fig, name, footnote):
    fig.text(.02, .015, footnote, fontsize=8, color='#444444', va='bottom')
    fig.savefig(HERE / f'{name}.png', dpi=200)
    fig.savefig(HERE / f'{name}.svg')
    plt.close(fig)


if __name__ == '__main__':
    assert roc_auc_score([0, 1], [.5, .5]) == .5
    assert not (HERE / 'figure_data.csv').exists(), 'Use a fresh directory; results are not overwritten.'
    audit = read(WORK / 'audit_followup_20261025/independent_metric_replay.json')
    domain = read(WORK / 'audit_followup_20261025/cross_domain_diagnostic.json')
    original = read(WORK / 'audit_handoff_20261024/retrain_evaluation.json')
    index = read(WORK / 'audit_handoff_20261024/records_index_exact.json')
    data_path = ROOT / 'docs/fingerprint_retrain_20261023/input/fingerprint_dataset.pt'
    INPUTS[str(data_path)] = hashlib.sha256(data_path.read_bytes()).hexdigest()
    assert INPUTS[str(data_path)] == index['summary']['dataset_sha256']
    data = torch.load(data_path, map_location='cpu', weights_only=False)
    records = data['records']
    assert len(records) == len(index['records']) == 331
    for q, (record, exported) in enumerate(zip(records, index['records'])):
        assert exported['graph_index'] == q
        assert all(record[k] == exported[k] for k in ['complex_id', 'family_id', 'receptor_id'])
    rid = np.array([r['receptor_id'] for r in records])
    fam = np.array([r['family_id'] for r in records])
    comp = np.stack([np.asarray(g[0])[:, 4:24].mean(0) for g in data['graphs']])
    comp /= np.maximum(np.linalg.norm(comp, axis=1, keepdims=True), 1e-12)
    size = np.log(np.array([len(g[0]) for g in data['graphs']], dtype=float))
    bases = {'composition': comp @ comp.T, 'size_log_distance': -abs(size[:, None]-size[None, :])}
    galleries, labels = {}, {}
    for st in ['all', 'same_family']:
        galleries[st], labels[st] = [], []
        for q in range(331):
            mask = np.arange(331) != q
            if st == 'same_family':
                mask &= (rid == rid[q]) | (fam == fam[q])
            gallery = np.flatnonzero(mask)
            galleries[st].append(gallery)
            labels[st].append((rid[gallery] == rid[q]).astype(int))

    def aucs(matrix, st):
        return np.array([roc_auc_score(y, matrix[q, g]) if len(set(y)) == 2 else np.nan
                         for q, (g, y) in enumerate(zip(galleries[st], labels[st]))])

    baseline_rows = {}
    for name, matrix in bases.items():
        baseline_rows[name] = {}
        for st in galleries:
            baseline_rows[name][st] = aucs(matrix, st)
            estimate = summarize(baseline_rows[name][st], fam)
            expected = original['baselines'][name][st]
            assert np.isclose(estimate['mean'], expected['mean'], atol=1e-12, rtol=0)
            assert np.allclose(estimate['ci95'], expected['ci95'], atol=1e-12, rtol=0)

    FIGURE = 'fig1_performance'
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8))
    fig.subplots_adjust(left=.20, right=.98, bottom=.18, top=.82, wspace=.55)
    methods = [('Composition', 'composition', GREEN), ('Graph size', 'size_log_distance', GRAY),
               ('Model: true labels', 'true', BLUE), ('Model: shuffled labels', 'shuffled', ORANGE)]
    for ax, st, title in zip(axes, ['all', 'same_family'], ['A  All negatives', 'B  Same-family negatives']):
        for y, (label, key, color) in enumerate(methods):
            est = summarize(baseline_rows[key][st], fam) if key in bases else audit['estimates'][key][st]['auroc']
            point(ax, est, y, color, f'{st}/{key}')
        ax.set(yticks=range(4), yticklabels=[m[0] for m in methods], xlim=(0, 1.025), ylim=(3.5, -.5), xlabel='Family-equal AUROC')
        ax.axvline(.5, color='#999999', ls='--', lw=1)
        ax.set_title(title, loc='left', fontweight='bold', pad=24)
        count = '226 queries / 48 families' if st == 'all' else '62 queries / 11 families'
        ax.text(0, 1.02, count, transform=ax.transAxes, fontsize=9, color=GRAY)
    save(fig, FIGURE, 'Dots: family-equal means. Bars: conditional 95% family-bootstrap intervals.\nPreviously exposed 331-graph development panel; legacy observed-sequence labels.')

    # Per-family detail: same eligibility and seed averaging as the primary estimate.
    models = []
    for seed in original['protocol']['seeds']:
        assembled = np.zeros((331, 331))
        for fold in range(5):
            path = ROOT / 'docs/fingerprint_retrain_20261023/run/true/model_run_v2' / f'sim_fold{fold}_seed{seed}.npy'
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            rel = f'true/model_run_v2/{path.name}'
            assert digest == original['artifact_sha256'][rel]
            INPUTS[str(path)] = digest
            matrix = np.load(path, allow_pickle=False)
            idx = [q for q, r in enumerate(index['records']) if r['fold'] == fold]
            assembled[idx] = matrix[idx]
        models.append(aucs(assembled, 'same_family'))
    values = np.stack(models, axis=1)
    valid = np.isfinite(values).all(axis=1)
    model_mean = np.full(331, np.nan)
    model_mean[valid] = values[valid].mean(axis=1)
    details = []
    for family in sorted(set(fam[valid])):
        selected = valid & (fam == family)
        a, b = model_mean[selected].mean(), baseline_rows['composition']['same_family'][selected].mean()
        details.append(dict(family_id=str(family), n_queries=int(selected.sum()), model=float(a), composition=float(b), delta=float(a-b)))
    assert len(details) == 11 and sum(r['n_queries'] for r in details) == 62
    assert np.isclose(np.mean([r['delta'] for r in details]), original['primary']['mean'], atol=1e-12, rtol=0)
    with (HERE / 'family_detail.csv').open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(details[0])); writer.writeheader(); writer.writerows(details)
    FIGURE = 'fig2_family_detail'
    fig, ax = plt.subplots(figsize=(9, 6.4))
    fig.subplots_adjust(left=.33, right=.96, bottom=.17, top=.86)
    for y, row in enumerate(details):
        ax.plot([row['composition'], row['model']], [y, y], color='#BFC7CD', lw=1.5)
        ax.scatter(row['composition'], y, color=GREEN, s=35, label='Composition' if y == 0 else None, zorder=3)
        ax.scatter(row['model'], y, color=BLUE, s=35, label='True-label model' if y == 0 else None, zorder=3)
    ax.axvline(.5, color='#999999', ls='--', lw=1)
    ax.set(yticks=range(11), yticklabels=[f"{r['family_id']}  (q={r['n_queries']})" for r in details], xlim=(-.025, 1.025), ylim=(10.7, -.7), xlabel='Within-family mean query AUROC')
    ax.set_title('Same-family retrieval: all 11 eligible families', loc='left', fontweight='bold', pad=30)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.01), ncol=2, frameon=False, fontsize=9)
    save(fig, FIGURE, 'One row per legacy family, in family-ID order. Model values average all three seeds within each query.\nLines connect paired estimates; they are not confidence intervals. Families are not assumed source-independent.')

    FIGURE = 'fig3_stability'
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3))
    fig.subplots_adjust(left=.18, right=.97, bottom=.21, top=.83, wspace=.52)
    for ax, st, title in zip(axes, ['all', 'same_family'], ['A  Full-gallery rankings', 'B  Same-family rankings']):
        for y, arm, color in [(0, 'true', BLUE), (1, 'shuffled', ORANGE)]:
            point(ax, audit['estimates'][arm][st]['spearman'], y, color, f'{st}/{arm}')
        ax.axvline(0, color='#999999', ls='--', lw=1)
        ax.set(yticks=[0, 1], yticklabels=['True labels', 'Shuffled labels'], xlim=(-.12, .67), ylim=(1.6, -.6), xlabel='Family-equal ranking Spearman')
        ax.set_title(title, loc='left', fontweight='bold', pad=26)
        ax.text(0, 1.02, '331 queries / 110 families' if st == 'all' else '62 queries / 11 families', transform=ax.transAxes, fontsize=9, color=GRAY)
    save(fig, FIGURE, 'Conditional 95% family-bootstrap intervals; three seed pairs averaged within query.\nRanking reproducibility does not establish retrieval gain over composition.')

    FIGURE = 'fig4_cross_domain'
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    fig.subplots_adjust(left=.18, right=.97, bottom=.23, top=.80, wspace=.52)
    for ax, key, title in zip(axes, ['all_negatives', 'same_family_negatives'], ['A  All negatives', 'B  Same-family negatives']):
        result = domain['results'][key]
        for y, name, color in [(0, 'composition', GREEN), (1, 'corrected_model', BLUE)]:
            point(ax, result[name]['cross_minus_same'], y, color, f'{key}/{name}')
        ax.axvline(0, color='#999999', ls='--', lw=1)
        ax.set(yticks=[0, 1], yticklabels=['Composition', 'True-label model'], xlim=(-.21, .16), ylim=(1.6, -.6), xlabel='AUROC: cross-domain minus same-domain')
        ax.set_title(title, loc='left', fontweight='bold', pad=26)
        ax.text(0, 1.02, f"{result['n_common_queries']} common queries / {result['n_common_groups']} families", transform=ax.transAxes, fontsize=9, color=GRAY)
    save(fig, FIGURE, 'Exploratory paired comparison: identical negatives, common eligible queries, legacy labels.\nConditional 95% family-bootstrap intervals. Domain is observational and confounded with source/construct.')
    with (HERE / 'figure_data.csv').open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ROWS[0])); writer.writeheader(); writer.writerows(ROWS)
    for path, expected in INPUTS.items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected, f'Input changed: {path}'
    report = dict(status='PASS', n_figures=4, n_families=11, n_queries=62, size_baseline_independently_replayed=True,
                  primary_delta=np.mean([r['delta'] for r in details]), source_sha256=INPUTS)
    (HERE / 'figure_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'source_sha256'}, indent=2))
