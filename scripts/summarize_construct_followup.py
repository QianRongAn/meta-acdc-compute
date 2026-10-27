"""Join adjudicated identities with fixed native QC; never infer experimental negatives."""
from collections import defaultdict
from pathlib import Path
import argparse
import copy
import hashlib
import itertools
import json

from scripts.summarize_cross_source_extension import qualified_rows, control_groups
from meta_acdc.structure.source_balanced_retrieval import source_balanced_blocks, joint_components
from meta_acdc.structure.receptor_identity_audit import seqres_sequences

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'docs/fingerprint_construct_followup_20261027'
PREVIOUS = ROOT / 'docs/fingerprint_positive_controls_20261027'


def key(row):
    return row['pdb'], row['alpha_chain'], row['beta_chain']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--identity', type=Path, required=True)
    ap.add_argument('--qc', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    old = json.loads((PREVIOUS/'source_balanced_design.json').read_text())
    adjudication = json.loads(args.identity.read_text())
    rows = {key(r): copy.deepcopy(r) for r in old['rows']}
    for row in adjudication['rows']:
        assert rows[key(row)]['status'] == 'quarantined'
        rows[key(row)] = row
    identity_rows = [r for r in rows.values() if r['status'] == 'identity_admitted_structure_QC_pending']
    qc_old = json.loads((PREVIOUS/'native_atomic_qc.json').read_text())
    qc_new = json.loads(args.qc.read_text())
    passed = qualified_rows(identity_rows, qc_old['rows'] + qc_new['rows'])
    assert len({key(r) for r in passed}) == len(passed)
    controls = control_groups(passed)
    cross_pairs = []
    for a, b in itertools.combinations(passed, 2):
        if a['receptor_id'] == b['receptor_id'] and a['source_study_key'] and b['source_study_key'] and a['source_study_key'] != b['source_study_key']:
            cross_pairs.append(dict(pdb_a=a['pdb'], pdb_b=b['pdb']))
    blocks = source_balanced_blocks(passed)
    strict = source_balanced_blocks(passed, True)
    eligible_contexts = {b['V_context'] for b in blocks}
    eligible_components = [c for c in joint_components(passed) if any(r['V_context'] in eligible_contexts for r in c)]
    contexts = defaultdict(list)
    for row in identity_rows:
        if row['source_study_key']:
            contexts[(row['source_study_key'],row['V_context'])].append(row)
    identity_controls = []
    for (study, context), members in sorted(contexts.items()):
        if len({r['receptor_id'] for r in members}) < 2:
            continue
        pdbs = sorted({r['pdb'] for r in members})
        identity_controls.append(dict(study=study,V_context=context,pdbs=pdbs,
            independent_natural_clonotypes=False,
            reason='6218 parent and S110C derivative' if '7n5c' in pdbs else 'JDI parent and affinity-engineered derivative'))
    queries = json.loads((BASE/'exact_search_plan.json').read_text())
    searches = defaultdict(dict)
    for query in queries:
        retry = BASE/'sources/additional'/query['file']
        path = retry if retry.exists() else BASE/'sources/exact_search'/query['file']
        response = json.loads(path.read_text())
        if 'result_set' not in response or response['result_type'] != 'polymer_entity':
            raise ValueError('sequence search incomplete: '+query['file'])
        searches[query['pdb']][query['role']] = dict(
            identifiers=[hit['identifier'] for hit in response['result_set']],
            path=str(path.relative_to(ROOT)), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            interpretation='discovery at 100% local alignment identity; complete paired sequences still required')
    query_report = []
    for pdb, roles in sorted(searches.items()):
        a = {x.split('_')[0].lower() for x in roles['alpha']['identifiers']}
        b = {x.split('_')[0].lower() for x in roles['beta']['identifiers']}
        query_report.append(dict(query_pdb=pdb, roles=roles, same_pdb_chain_hit_intersection=sorted(a & b)))
    # Inspect actual deposited termini of the most promising historical 6218 lead.
    parent = next(r for r in adjudication['rows'] if r['pdb']=='7n5p')
    comparisons = []
    for pdb, relative in [('3pqy','sources/retry/3pqy.pdb'), ('7n4k','sources/tertiary/7n4k.pdb')]:
        path = BASE/relative
        sequences = seqres_sequences(path.read_text())
        previous = next(r for r in old['rows'] if r['pdb']==pdb)
        chain = previous['alpha_chain']; full = sequences[chain]
        assert full == previous['alpha']['full_construct']
        beta = sequences[previous['beta_chain']]
        comparisons.append(dict(pdb=pdb, alpha_chain=chain, beta_chain=previous['beta_chain'],
            actual_alpha_start=full[:18], complete_parent_alpha_start=parent['alpha']['variable'][:18],
            complete_parent_variable_is_exact_substring=parent['alpha']['variable'] in full,
            complete_parent_beta_is_exact_substring=parent['beta']['variable'] in beta,
            original_quarantine_reasons=previous['reasons'],
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            conclusion='same biological lineage is not complete deposited paired-variable equality; no N-terminal residue is prepended'))
    summary = dict(n_identity_admitted_pdbs=len({r['pdb'] for r in identity_rows}),
                   n_complete_paired_ids=len({r['receptor_id'] for r in identity_rows}),
                   n_QC_graphs=len(passed), n_QC_pdbs=len({r['pdb'] for r in passed}),
                   n_QC_complete_paired_ids=len({r['receptor_id'] for r in passed}),
                   n_QC_cross_study_pairs=len(cross_pairs),
                   n_identity_same_study_same_V_groups=len(identity_controls),
                   n_QC_same_study_same_V_control_groups=len(controls),
                   n_QC_independent_natural_control_groups=0,
                   n_QC_source_balanced_blocks=len(blocks), n_QC_author_disjoint_balanced_blocks=len(strict),
                   n_eligible_joint_components=len(eligible_components))
    result = dict(summary=summary, status='blocked_before_protocol_freeze',
        primary_role_receipts=str((BASE/('construct_adjudication_final.json' if any(r['pdb']=='8cx4' for r in adjudication['rows']) else 'construct_adjudication.json')).relative_to(ROOT)),
        admitted_pdbs=sorted({r['pdb'] for r in identity_rows}),
        QC_pass_units=[{k:r[k] for k in ['pdb','alpha_chain','beta_chain','receptor_id','V_context','source_study_key','native_graph_sha256']} for r in passed],
        same_source_same_V_identity_groups=identity_controls, QC_same_source_same_V_groups=controls,
        sequence_search_refresh=query_report, deposited_6218_comparisons=comparisons,
        cross_study_pairs=cross_pairs, source_balanced_blocks=blocks, author_disjoint_blocks=strict,
        blockers=['no fixed-QC passing cross-study complete alpha/beta positive pair',
                  'no qualified source-balanced query block or independent natural same-source same-V structural control',
                  'fewer than three eligible disjoint source/family components'],
        protocol_frozen=False, training_release_created=False, model_training_started=False,
        independent_acceptance_accessed=False,
        scope='sequence-supported structural-development candidates only; not experimental nonbinding labels')
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
