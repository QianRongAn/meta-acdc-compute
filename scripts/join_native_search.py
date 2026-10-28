"""Join paired-chain discovery hits without turning failed requests into negatives."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def parse_hits(status, payload):
    if status == 204 and not payload:
        return []
    if status != 200:
        raise ValueError('search request did not succeed')
    data = json.loads(payload)
    if data.get('result_type') != 'polymer_entity' or not isinstance(data.get('total_count'), int):
        raise ValueError('response is not a polymer-entity search result')
    hits = data.get('result_set', [])
    if len(hits) != data['total_count']:
        raise ValueError('search response incomplete despite return_all_hits')
    return [hit['identifier'] for hit in hits]


def pair_candidates(alpha, beta):
    if alpha is None or beta is None:
        return None
    return sorted({x.rsplit('_',1)[0].lower() for x in alpha} &
                  {x.rsplit('_',1)[0].lower() for x in beta})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.base = args.base.resolve()
    if args.out.exists():
        raise FileExistsError(args.out)
    plan = json.loads((args.base/'search_plan.json').read_text())
    records = {}
    for directory in ['search','search_retry']:
        manifest = args.base/'sources'/directory/'manifest.json'
        if manifest.exists():
            for row in json.loads(manifest.read_text()):
                records[row['file']] = (manifest.parent,row)
    queries,lookup,failures = [],{},[]
    for query in plan['queries']:
        directory,row = records[query['file']]
        path = directory/query['file']
        payload = path.read_bytes() if path.exists() else b''
        if 'sha256' in row:
            assert hashlib.sha256(payload).hexdigest() == row['sha256']
        result = dict(file=query['file'],role=query['role'],source=str(path.relative_to(ROOT)),http_status=row['status'])
        try:
            hits = parse_hits(row['status'],payload)
            result.update(status='search_completed',hits=hits,n_entities=len(hits))
        except (ValueError,KeyError,TypeError) as error:
            hits = None
            result.update(status='unavailable_not_negative',reason=str(error))
            failures.append(query['file'])
        queries.append(result)
        lookup[(query['role'],query['sequence'])] = hits
    targets = []
    for target in plan['targets']:
        hits = pair_candidates(lookup[('alpha',target['sequences']['alpha'])],lookup[('beta',target['sequences']['beta'])])
        targets.append(dict(target_id=target['target_id'],anchor=target['anchor'],
                            full_boundary_validated=target['full_boundary_validated'],
                            anchor_study=target['anchor_study'],candidate_pdbs=hits,
                            status='requires_primary_pair_boundary_source_and_atomic_validation' if hits is not None else 'incomplete_search_not_negative'))
    result = dict(n_queries=len(queries),n_successful_queries=len(queries)-len(failures),failed_query_files=failures,
                  n_targets=len(targets),targets=targets,queries=queries,
                  discovery_pdbs=sorted({p for target in targets for p in target['candidate_pdbs'] or []}),
                  interpretation='Per-PDB intersection only discovers leads; neither 100% local match nor different PDB establishes complete paired receptor or source independence.',
                  protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['n_queries','n_successful_queries','failed_query_files','n_targets','discovery_pdbs']},indent=2))


if __name__=='__main__':
    main()
