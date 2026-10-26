"""Independent read-only replay of frozen matrix metrics; no model training."""
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import torch
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path('/home/administrator/docs/meta-acdc-fingerprint-review')
OUT = Path(__file__).resolve().parent


def summarize(values, families):
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    counts = np.isfinite(values).sum(axis=1)
    query = np.nansum(values, axis=1) / np.maximum(counts, 1)
    query[counts == 0] = np.nan
    means = []
    for family in sorted(set(families)):
        subset = query[np.asarray(families) == family]
        subset = subset[np.isfinite(subset)]
        if len(subset):
            means.append(subset.mean())
    assert len(means) >= 2
    means = np.asarray(means)
    rng = np.random.default_rng(20261023)
    bootstrap = means[rng.integers(0, len(means), (2000, len(means)))].mean(axis=1)
    return dict(mean=float(means.mean()), ci95=np.quantile(bootstrap, [.025, .975]).tolist(),
                n_queries=int((counts > 0).sum()), n_groups=len(means))


def auc_rows(matrix, galleries, labels):
    return np.array([roc_auc_score(y, matrix[q, gallery]) if len(set(y)) == 2 else np.nan
                     for q, (gallery, y) in enumerate(zip(galleries, labels))])


def main():
    # Tie handling, seed averaging before group averaging, and missing queries.
    assert roc_auc_score([0, 1], [1, 1]) == .5
    test = summarize([[0, 1], [1, 1], [np.nan, np.nan], [0, 0]], ['a', 'a', 'a', 'b'])
    assert test['mean'] == .375 and test['n_queries'] == 3
    source_paths = set()

    def read_bytes(path):
        source_paths.add(path)
        return path.read_bytes()

    run = ROOT / 'docs/fingerprint_retrain_20261023/run'
    original = json.loads(read_bytes(run / 'evaluation.json'))
    dataset_path = ROOT / 'docs/fingerprint_retrain_20261023/input/fingerprint_dataset.pt'
    dataset_hash = hashlib.sha256(read_bytes(dataset_path)).hexdigest()
    assert dataset_hash == original['protocol']['dataset_sha256']
    data = torch.load(dataset_path, weights_only=False, map_location='cpu')
    records = data['records']
    folds = json.loads(read_bytes(ROOT / original['protocol']['fold_assignment']))
    assert folds['graph_records'] == [r['complex_id'] for r in records]
    n = len(records)
    assert n == 331
    families = [r['family_id'] for r in records]
    seeds = original['protocol']['seeds']
    rid = np.array([r['receptor_id'] for r in records])
    fam = np.array(families)
    population = {}
    for stratum in ['all', 'same_family']:
        galleries, labels = [], []
        for q in range(n):
            mask = np.arange(n) != q
            if stratum == 'same_family':
                mask &= (fam == fam[q]) | (rid == rid[q])
            gallery = np.flatnonzero(mask)
            galleries.append(gallery)
            labels.append((rid[gallery] == rid[q]).astype(int))
        population[stratum] = (galleries, labels)

    composition = np.stack([np.asarray(g[0])[:, 4:24].mean(0) for g in data['graphs']])
    composition /= np.maximum(np.linalg.norm(composition, axis=1, keepdims=True), 1e-12)
    base = composition @ composition.T
    values, ranking, results = {}, {}, {}
    for arm in ['true', 'shuffled']:
        matrices = []
        for seed in seeds:
            assembled = np.zeros((n, n))
            for fold in range(5):
                tag = f'fold{fold}_seed{seed}'
                directory = run / arm / 'model_run_v2'
                simpath, embpath = directory / f'sim_{tag}.npy', directory / f'emb_{tag}.npy'
                for path in [simpath, embpath]:
                    digest = hashlib.sha256(read_bytes(path)).hexdigest()
                    assert digest == original['artifact_sha256'][str(path.relative_to(run))]
                sim = np.load(simpath, allow_pickle=False)
                emb = np.load(embpath, allow_pickle=False).astype(float)
                assert sim.shape == (n, n) and np.isfinite(sim).all()
                assert np.allclose(sim, emb @ emb.T, atol=2e-6, rtol=0)
                idx = np.flatnonzero(np.asarray(folds['per_graph_fold']) == fold)
                assembled[idx] = sim[idx]
            matrices.append(assembled)
        results[arm], values[arm], ranking[arm] = {}, {}, {}
        for st, (galleries, labels) in population.items():
            scores = np.stack([auc_rows(matrix, galleries, labels) for matrix in matrices], axis=1)
            stable = np.full((n, 3), np.nan)
            for q, gallery in enumerate(galleries):
                if st == 'all' or np.isfinite(scores[q, 0]):
                    for col, (a, b) in enumerate(combinations(range(3), 2)):
                        stable[q, col] = spearmanr(matrices[a][q, gallery], matrices[b][q, gallery]).statistic
            values[arm][st], ranking[arm][st] = scores, stable
            results[arm][st] = dict(auroc=summarize(scores, families), spearman=summarize(stable, families))
            for metric, expected in [('auroc', original['family_equal_seed_averaged'][arm][st]),
                                     ('spearman', original['ranking_stability'][arm][st])]:
                actual = results[arm][st][metric]
                assert np.isclose(actual['mean'], expected['mean'], atol=1e-12, rtol=0)
                assert np.allclose(actual['ci95'], expected['ci95'], atol=1e-12, rtol=0)
    contrasts = {}
    for st, (galleries, labels) in population.items():
        baseline = auc_rows(base, galleries, labels)
        contrasts[st] = dict(composition=summarize(baseline, families),
            true_minus_composition=summarize(values['true'][st] - baseline[:, None], families),
            true_minus_shuffled=summarize(values['true'][st] - values['shuffled'][st], families),
            exploratory_stability_difference=summarize(ranking['true'][st] - ranking['shuffled'][st], families))
    for key, expected in [('true_minus_composition', original['primary']),
                          ('true_minus_shuffled', original['true_minus_shuffled']['same_family'])]:
        assert np.isclose(contrasts['same_family'][key]['mean'], expected['mean'], atol=1e-12, rtol=0)
        assert np.allclose(contrasts['same_family'][key]['ci95'], expected['ci95'], atol=1e-12, rtol=0)

    alignment = json.loads((OUT / 'procrustes_train_only/rashomon_results.json').read_text())
    align_checks = []
    for fold in alignment['fold_results']:
        # Locate the supplement without depending on its historical parent key.
        supplement = next(v for v in fold.values() if isinstance(v, dict) and 'fit_graph_indices' in v)
        train, test = set(supplement['fit_graph_indices']), set(supplement['evaluation_graph_indices'])
        assert not train & test and train | test == set(range(n))
        assert supplement['maximum_cosine_change_float64'] < 1e-10
        align_checks.append(dict(n_train=len(train), n_test=len(test), disjoint=True,
                                 max_cosine_change=supplement['maximum_cosine_change_float64']))
    for path, expected in alignment['source_sha256'].items():
        assert hashlib.sha256(read_bytes(Path(path))).hexdigest() == expected
    report = dict(status='PASS', n_instances=30, n_graphs=n, estimates=results, contrasts=contrasts,
                  alignment_checks=align_checks,
                  caveat='Fixed models/gallery, family conditional bootstrap; exploratory stability contrast is not a validation claim.',
                  source_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source_paths)})
    target = OUT / 'independent_metric_replay.json'
    with target.open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
    print(json.dumps({'status': report['status'], 'contrasts': contrasts, 'alignment': align_checks}, indent=2))


if __name__ == '__main__':
    main()
