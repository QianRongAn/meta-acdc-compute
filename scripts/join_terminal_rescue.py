"""Replay fragment searches and verify exact substrings in deposited entities."""
import argparse
import hashlib
import json
from urllib.parse import parse_qs, urlencode, urlsplit

from scripts.join_native_search import ROOT, pair_candidates, parse_hits
from scripts.prepare_terminal_rescue import BASE
from scripts.validate_native_search import exact_entity_pairs


def load_receipt(path, receipt):
    payload = path.read_bytes() if path.is_file() else b""
    assert hashlib.sha256(payload).hexdigest() == receipt["sha256"], path
    return payload


def search_hits(plan):
    attempts = {}
    for folder in sorted((BASE / "sources").glob("search*")):
        manifest = folder / "manifest.json"
        if manifest.is_file():
            for receipt in json.loads(manifest.read_text()):
                attempts.setdefault(receipt["file"], []).append((folder / receipt["file"], receipt))
    results = {}
    for query in plan["queries"]:
        cached = query["cached_response"]
        receipts = [(ROOT / cached["file"], cached)] if cached else attempts.get(query["file"], [])
        result = {"query": query["file"], "hits": None, "attempts": []}
        for path, receipt in receipts:
            url = receipt.get("url") or receipt["URL"]
            assert json.loads(parse_qs(urlsplit(url).query)["json"][0]) == query["query"]
            payload = load_receipt(path, receipt)
            attempt = {"file": str(path.relative_to(ROOT)), "status": receipt["status"]}
            try:
                hits = parse_hits(int(receipt["status"]), payload)
            except (ValueError, KeyError, TypeError) as error:
                attempt["reason"] = str(error)
            else:
                assert result["hits"] is None or result["hits"] == hits
                result["hits"] = hits
            result["attempts"].append(attempt)
        results[query["sequence"]] = result
    return results


def metadata_requests(pdbs):
    requests = []
    for start in range(0, len(pdbs), 25):
        ids = json.dumps([p.upper() for p in pdbs[start:start + 25]])
        query = """{ entries(entry_ids: %s) {
          rcsb_id struct { title }
          rcsb_primary_citation { title pdbx_database_id_DOI pdbx_database_id_PubMed rcsb_authors year }
          polymer_entities { rcsb_id rcsb_polymer_entity { pdbx_description }
          entity_poly { pdbx_seq_one_letter_code_can }
          rcsb_polymer_entity_container_identifiers { auth_asym_ids } }
        } }""" % ids
        requests.append([f"batch_{start // 25:02}.json",
                         "https://data.rcsb.org/graphql?" + urlencode({"query": query})])
    return requests


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--requests")
    args = parser.parse_args()
    plan = json.loads((BASE / "search_plan.json").read_text())
    search = search_hits(plan)
    metadata = {}
    folder = BASE / "sources/metadata"
    if (folder / "manifest.json").is_file():
        for receipt in json.loads((folder / "manifest.json").read_text()):
            payload = load_receipt(folder / receipt["file"], receipt)
            if receipt["status"] == "200":
                response = json.loads(payload)
                assert not response.get("errors"), response
                for entry in response["data"]["entries"]:
                    if entry:
                        metadata[entry["rcsb_id"].lower()] = entry
    rows = []
    for target in plan["targets"]:
        alpha, beta = (target["sequences"][role] for role in ("alpha", "beta"))
        pdbs = pair_candidates(search[alpha]["hits"], search[beta]["hits"])
        exact = {}
        for pdb in pdbs or []:
            if pdb in metadata:
                matches = exact_entity_pairs(metadata[pdb], alpha, beta)
                if matches:
                    exact[pdb] = [[a["rcsb_id"], b["rcsb_id"]] for a, b in matches]
        rows.append({"target_id": target["target_id"], "candidate_pdbs": pdbs, "exact_entity_matches": exact})
    pdbs = sorted({pdb for row in rows for pdb in row["candidate_pdbs"] or []})
    counts = {
        "queries": len(search),
        "successful_queries": sum(r["hits"] is not None for r in search.values()),
        "complete_pair_searches": sum(r["candidate_pdbs"] is not None for r in rows),
        "fragment_pairs": len(rows),
        "candidate_pdbs": len(pdbs),
        "metadata_covered_pdbs": sum(p in metadata for p in pdbs),
        "exact_fragment_pdbs": len({p for r in rows for p in r["exact_entity_matches"]}),
    }
    output = dict(counts=counts, queries=list(search.values()), rows=rows, candidate_pdbs=pdbs,
                  missing_metadata=sorted(set(pdbs) - metadata.keys()),
                  interpretation="Distinct deposited-entity substring hits are leads only. Incomplete searches remain null.")
    with open(args.out, "x") as stream:
        stream.write(json.dumps(output, indent=2) + "\n")
    if args.requests:
        with open(args.requests, "x") as stream:
            stream.write(json.dumps(metadata_requests(pdbs), indent=2) + "\n")
    print(json.dumps(counts))


if __name__ == "__main__":
    main()
