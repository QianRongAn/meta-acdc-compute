"""Validate fresh sequence/source metadata against all paired discovery leads."""
from pathlib import Path
import argparse
import hashlib
import itertools
import json

from scripts.curate_cross_source_receptors import study_key
from meta_acdc.structure.source_balanced_retrieval import author_keys

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'docs/fingerprint_native_search_20261029'


def exact_entity_pairs(entry, alpha, beta):
    entities = entry.get('polymer_entities') or []
    def matches(sequence):
        return [e for e in entities if (e.get('entity_poly') or {}).get('pdbx_seq_one_letter_code_can','').count(sequence) == 1]
    return [(a,b) for a,b in itertools.product(matches(alpha),matches(beta)) if a['rcsb_id'] != b['rcsb_id']]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    plan = json.loads((BASE/'search_plan.json').read_text())
    discovery = json.loads((BASE/'search_results.json').read_text())
    assert discovery['n_successful_queries'] == len(plan['queries'])
    metadata = {}
    manifest = json.loads((BASE/'sources/metadata/manifest.json').read_text())
    source_hashes = {}
    for source in manifest:
        path = BASE/'sources/metadata'/source['file']
        assert source['status'] == 200 and hashlib.sha256(path.read_bytes()).hexdigest() == source['sha256']
        data = json.loads(path.read_text())
        if data.get('errors'):
            raise ValueError('fresh metadata query incomplete')
        for entry in data['data']['entries']:
            assert entry is not None
            metadata[entry['rcsb_id'].lower()] = entry
        source_hashes[str(path.relative_to(ROOT))] = source['sha256']
    assert set(metadata) == set(discovery['discovery_pdbs'])
    old = json.loads((ROOT/'docs/fingerprint_positive_controls_20261027/source_balanced_design.json').read_text())
    identities = {(r['pdb'],r['alpha_chain'],r['beta_chain']):r for r in old['rows']}
    correction = ROOT/'docs/fingerprint_construct_followup_20261027'
    identities.update({(r['pdb'],r['alpha_chain'],r['beta_chain']):r
                       for r in json.loads((correction/'construct_adjudication_final.json').read_text())['rows']})
    qc_rows = json.loads((correction/'joint_readiness_final.json').read_text())['QC_pass_units']
    qc_keys = {(r['pdb'],r['alpha_chain'],r['beta_chain']) for r in qc_rows}
    atomic_failures = {}
    scope_failures = {}
    for report in [ROOT/'docs/fingerprint_positive_controls_20261027/native_atomic_qc.json',
                   correction/'native_qc_final/native_atomic_qc.json']:
        for row in json.loads(report.read_text())['rows']:
            if 'alpha_chain' in row and 'beta_chain' in row and row.get('status') != 'native_graph_QC_pass':
                atomic_failures[(row['pdb'],row['alpha_chain'],row['beta_chain'])] = row.get('reasons',[])
            if row.get('stage') == 'unit_extraction' and row.get('n_units') == 0:
                scope_failures[row['pdb']] = row.get('excluded',[])
    records = []
    for target,lead in zip(plan['targets'],discovery['targets']):
        assert target['target_id'] == lead['target_id']
        anchor_citation = metadata[target['anchor']].get('rcsb_primary_citation') or {}
        fresh_anchor_study = study_key(anchor_citation)
        for pdb in lead['candidate_pdbs']:
            entry = metadata[pdb]
            pairs = exact_entity_pairs(entry,target['sequences']['alpha'],target['sequences']['beta'])
            citation = entry.get('rcsb_primary_citation') or {}
            fresh_study = study_key(citation)
            cross_study = bool(pdb != target['anchor'] and fresh_study and fresh_anchor_study and fresh_study != fresh_anchor_study)
            pair_reports, qualified = [],False
            for a,b in pairs:
                for ac,bc in itertools.product(a['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'],
                                                b['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']):
                    if ac == bc:
                        continue
                    key = (pdb,ac,bc)
                    stored = identities.get(key)
                    sequences_unchanged = bool(stored and all(stored[role]['full_construct'] == entity['entity_poly']['pdbx_seq_one_letter_code_can']
                                                              for role,entity in [('alpha',a),('beta',b)]))
                    if stored and sequences_unchanged:
                        reasons = stored.get('reasons',[]) if stored['status'] == 'quarantined' else []
                        identity_valid = stored['status'] == 'identity_admitted_structure_QC_pending'
                        atomic_valid = key in qc_keys and identity_valid and target['full_boundary_validated']
                        if identity_valid and not atomic_valid:
                            reasons = atomic_failures.get(key,[]) or [x.get('reason','scope exclusion') for x in scope_failures.get(pdb,[])] or ['no fixed-QC-passing contact-assigned unit in the frozen candidate inventory']
                    else:
                        reasons = ['new chain pair or fresh construct differs; current complete pairing/coordinates must be curated']
                        identity_valid,atomic_valid = False,False
                    qualified |= atomic_valid
                    pair_reports.append(dict(alpha_entity=a['rcsb_id'],beta_entity=b['rcsb_id'],alpha_chain=ac,beta_chain=bc,
                        constructs_unchanged=sequences_unchanged,identity_gate_pass=identity_valid,existing_fixed_QC_pass=atomic_valid,
                        reasons=reasons,unreviewed_pair=not stored))
            if not pairs:
                decision = 'local_alignment_only_no_exact_paired_query_substrings'
            elif not target['full_boundary_validated']:
                decision = 'fragment_discovery_not_complete_identity'
            elif qualified and cross_study:
                decision = 'requires_new_protocol_readiness_review'
            elif qualified:
                decision = 'existing_fixed_QC_candidate_no_new_study'
            else:
                decision = 'complete_substrings_do_not_pass_existing_identity_or_atomic_scope'
            ak,bk = author_keys({'primary_citation':anchor_citation}),author_keys({'primary_citation':citation})
            records.append(dict(target_id=target['target_id'],anchor=target['anchor'],pdb=pdb,
                query_is_complete_identity=target['full_boundary_validated'],fresh_title=(entry.get('struct') or {}).get('title'),
                fresh_citation=citation,fresh_study=fresh_study,cross_study=cross_study,
                frozen_anchor_study=target['anchor_study'],fresh_anchor_study=fresh_anchor_study,
                n_declared_polymer_entities=len(entry.get('polymer_entities') or []),
                author_disjoint=bool(ak and bk and not ak&bk),
                full_query_pair_on_distinct_entities=bool(pairs),chain_pair_checks=pair_reports,
                decision=decision))
    cross = [r for r in records if r['query_is_complete_identity'] and r['cross_study'] and r['full_query_pair_on_distinct_entities']]
    eligible = [r for r in records if r['decision'] == 'requires_new_protocol_readiness_review']
    report = dict(n_fresh_metadata_pdbs=len(metadata),n_target_pdb_leads=len(records),
                  n_complete_cross_study_leads=len(cross),n_new_qualified_cross_study_leads=len(eligible),
                  complete_cross_study_leads=cross,rows=records,source_sha256=source_hashes,
                  interpretation='Exact deposited substrings still require current chain pairing and immutable atomic QC. Cartesian entity/chain combinations are leads, not contact-proven pairs.',
                  training_release_created=False,evaluation_protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    args.out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(n_fresh_metadata_pdbs=len(metadata),n_target_pdb_leads=len(records),
                         n_complete_cross_study_leads=len(cross),n_new_qualified_cross_study_leads=len(eligible),
                         cross_study_pdbs=[r['pdb'] for r in cross]),indent=2))


if __name__=='__main__':
    main()
