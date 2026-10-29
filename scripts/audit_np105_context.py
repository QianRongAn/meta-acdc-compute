"""Audit unchanged NP105 ORF contexts against archived neural outcomes; no renumbering."""
from pathlib import Path
import collections, gzip, hashlib, json
ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "docs/fingerprint_np105_recovery_20261029"
OUT = ROOT / "docs/fingerprint_np105_context_20261030"
def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()
def domain_context(row):
    boundary = row["boundary"]
    start, end = boundary["start0"], boundary["end0_exclusive"]
    full, domain = row["ORF"], boundary["sequence"]
    if not domain or full[start:end] != domain or full.count(domain) != 1:
        raise ValueError("domain must be the exact unique original observed interval")
    return full[:start], domain, full[end:]
def outcome(row):
    n = row["neural"]
    mapping = tuple((tuple(p), aa) for p, aa in (n.get("numbering") or []) if aa != "-")
    return n.get("chain_type"), n.get("error"), mapping

def main():
    classic = json.loads(gzip.decompress((OLD / "evidence/whole_classic_numbering.json.gz").read_bytes()))
    neural = json.loads(gzip.decompress((OLD / "evidence/whole_neural_numbering.json.gz").read_bytes()))
    lookup = {x["name"]: x for x in neural["rows"]}
    groups = collections.defaultdict(list)
    for cell in classic["cells"]:
        for role in ("TRA", "TRB"):
            for row in cell[role]:
                if row["complete"]:
                    left, domain, right = domain_context(row)
                    nn = lookup[row["name"]]
                    assert sha(row["ORF"]) == nn["ORF_sha256"]
                    groups[(role, domain)].append((row, nn, left, right))
    records = []
    for (role, domain), members in groups.items():
        classic_maps = {tuple((tuple(p), aa) for p, aa in row["numbered"][0][0] if aa != "-") for row, _, _, _ in members}
        assert len(classic_maps) == 1, "same-domain classic maps also vary; keep this unresolved"
        distinct = collections.defaultdict(list)
        by_full = collections.defaultdict(set)
        rows = []
        for row, nn, left, right in members:
            sign = outcome(nn)
            distinct[sign].append(row["name"])
            by_full[row["ORF"]].add(sign)
            rows.append(dict(name=row["name"], barcode=row["barcode"], ORF_sha256=sha(row["ORF"]), N_flank_sha256=sha(left), C_flank_sha256=sha(right), N_flank_length=len(left), C_flank_length=len(right), neural_chain=nn["neural"].get("chain_type"), neural_error=nn["neural"].get("error"), neural_score=nn["neural"].get("score"), identical_non_gap_map=nn["non_gap_position_map_identical"], role_agrees=nn["neural_role_matches_source"]))
        assert all(len(v) == 1 for v in by_full.values())
        matching = sum(x[1]["non_gap_position_map_identical"] and x[1]["neural_role_matches_source"] for x in members)
        records.append(dict(role=role, observed_domain=domain, domain_sha256=sha(domain), ORF_occurrences=len(members), distinct_classic_non_gap_mappings=len(classic_maps), unique_full_ORFs=len(by_full), distinct_neural_role_error_or_non_gap_map_outcomes=len(distinct), same_domain_different_context_different_outcomes=len(distinct)>1, agrees_with_classic=matching, disagreements=len(members)-matching, members=rows))
    assert sum(g["ORF_occurrences"] for g in records)==394
    assert sum(g["disagreements"] for g in records)==55
    report=dict(scope="all 394 original complete ORF occurrences; original unmodified neural outputs", groups=records, n_observed_domain_proteins=len(records), n_complete_ORF_occurrences=394, n_disagreements_retained=55, groups_with_context_dependent_outcomes=sum(g["same_domain_different_context_different_outcomes"] for g in records), new_model_inference=False, cropped_inputs=0, synthetic_sequences=0, reference_bases_added=0, duplicate_full_ORF_predictions_reused_in_original_driver=True, duplicate_consistency_not_independent_repeat_evidence=True, interpretation="Within exact same-role observed variable-domain groups, neural outcomes vary with whole-ORF context. This does not identify which flank causes the effect, establish biological correctness of flanks, or validate a revised numbering protocol.", original_Q30_pairs=54, original_joint_numbered_Q30_pairs=26, new_identity_admissions=0, protocol_frozen=False, independent_acceptance_accessed=False)
    (OUT / "context_audit.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps([{k:v for k,v in g.items() if k not in ("members","observed_domain")} for g in records],indent=2))
if __name__=="__main__":
    main()
