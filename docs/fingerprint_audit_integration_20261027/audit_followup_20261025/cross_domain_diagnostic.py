"""Exploratory paired domain diagnostic under legacy observed-sequence labels."""
import hashlib
import json
from pathlib import Path

import numpy as np
import torch
from sklearn.metrics import roc_auc_score

from verify_metrics import ROOT, OUT, summarize


def main():
    target = OUT / 'cross_domain_diagnostic.json'
    assert not target.exists()
    path = ROOT / 'docs/fingerprint_retrain_20261023/input/fingerprint_dataset.pt'
    index_path = OUT.parent / 'audit_handoff_20261024/records_index_exact.json'
    index = json.loads(index_path.read_text(encoding='utf-8-sig'))
    assert hashlib.sha256(path.read_bytes()).hexdigest() == index['summary']['dataset_sha256']
    data = torch.load(path, map_location='cpu', weights_only=False)
    records = data['records']
    for i, (record, exported) in enumerate(zip(records, index['records'])):
        assert exported['graph_index'] == exported['matrix_row'] == exported['matrix_column'] == i
        for key in ['receptor_id', 'family_id', 'domain', 'complex_id']:
            assert record[key] == exported[key]
    assert len(records) == len(index['records']) == 331
    n = len(records)
    rid = np.array([r['receptor_id'] for r in records])
    fam = np.array([r['family_id'] for r in records])
    domain = np.array([r['domain'] for r in records])
    comp = np.stack([np.asarray(g[0])[:, 4:24].mean(0) for g in data['graphs']])
    comp /= np.maximum(np.linalg.norm(comp, axis=1, keepdims=True), 1e-12)
    run = ROOT / 'docs/fingerprint_retrain_20261023/run'
    evaluation = json.loads((run / 'evaluation.json').read_text())
    matrices = [comp @ comp.T]
    sources = {str(index_path): hashlib.sha256(index_path.read_bytes()).hexdigest(),
               str(path): hashlib.sha256(path.read_bytes()).hexdigest()}
    for seed in evaluation['protocol']['seeds']:
        assembled = np.zeros((n, n))
        for fold in range(5):
            simpath = run / 'true/model_run_v2' / f'sim_fold{fold}_seed{seed}.npy'
            digest = hashlib.sha256(simpath.read_bytes()).hexdigest()
            assert digest == evaluation['artifact_sha256'][str(simpath.relative_to(run))]
            sources[str(simpath)] = digest
            matrix = np.load(simpath, allow_pickle=False)
            idx = [i for i, r in enumerate(index['records']) if r['fold'] == fold]
            assembled[idx] = matrix[idx]
        matrices.append(assembled)
    result = {}
    for restriction in ['all_negatives', 'same_family_negatives']:
        scores = np.full((n, len(matrices), 2), np.nan)
        for q in range(n):
            negatives = np.flatnonzero(rid != rid[q])
            if restriction == 'same_family_negatives':
                negatives = negatives[fam[negatives] == fam[q]]
            positives = np.flatnonzero((rid == rid[q]) & (np.arange(n) != q))
            for col, is_same in enumerate([True, False]):
                selected = positives[(domain[positives] == domain[q]) == is_same]
                if len(selected) and len(negatives):
                    gallery = np.concatenate([selected, negatives])
                    y = np.concatenate([np.ones(len(selected)), np.zeros(len(negatives))])
                    for m, matrix in enumerate(matrices):
                        scores[q, m, col] = roc_auc_score(y, matrix[q, gallery])
        common = np.isfinite(scores).all(axis=(1, 2))
        n_groups = len(set(fam[common]))
        row = dict(n_common_queries=int(common.sum()), n_common_groups=n_groups,
                   same_only_eligible=int(np.isfinite(scores[:, 0, 0]).sum()),
                   cross_only_eligible=int(np.isfinite(scores[:, 0, 1]).sum()),
                   query_indices=np.flatnonzero(common).tolist(),
                   query_domains={str(d): int(((domain == d) & common).sum()) for d in set(domain)})
        if n_groups >= 2:
            for name, selection in [('composition', [0]), ('corrected_model', [1, 2, 3])]:
                a = scores[common][:, selection, :]
                row[name] = dict(same_domain=summarize(a[:, :, 0], fam[common]),
                                 cross_domain=summarize(a[:, :, 1], fam[common]),
                                 cross_minus_same=summarize(a[:, :, 1] - a[:, :, 0], fam[common]))
        else:
            row['status'] = 'insufficient_groups_for_summary'
        result[restriction] = row
    output = dict(status='PASS', protocol='Exploratory after seeing prior results; fixed identical negatives for both positive-domain conditions; paired common-query population; seed means then equal families.',
                  limitations='Legacy observed-chain labels, not validated complete receptor truth; domain is observational and confounded with source/construct. Not a causal domain or leakage experiment. No cross-fold pooled scores.',
                  results=result, source_sha256=sources)
    with target.open('x') as handle:
        json.dump(output, handle, indent=2, allow_nan=False)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
