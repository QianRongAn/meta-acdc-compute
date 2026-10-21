# Meta-ACDC preprint author-review package

This package turns the retained evidence into one bounded research manuscript. It is not a release of a validated binding or clinical model. The author requested a preprint and confirmed single authorship; personal and submission metadata will be supplied later.

- manuscript.md: complete English draft with abstract, methods, results, discussion, availability and disclosures.
- supplement.md: evidence tables, chronology, limits and reproduction instructions.
- output/pdf_author_review/: current rendered manuscript and supplementary PDFs for author review. Earlier PDF drafts remain in output/pdf/ and output/pdf_v2/ and later numbered draft directories.
- results_v2/: four data-derived figures, eight CSV tables, a summary and provenance receipt.
- atomic_bundle_v1/: coordinate extracts and an independent distance-count verification for all888 historical records.
- claim_evidence.json: supported and excluded claims.
- references.json: thirteen DOI/title-verified primary references.
- source_manifest.json: immutable core report-input hashes, pinned to e4d68de; context_manifest.json separately pins15eb5c6.
- analysis_plan.json: retrospective descriptive aggregation plan; no new model fitting or test scoring.
- submission_readiness.json: precise current state, unresolved author metadata and current eligibility requirement.

## Reproduce

From the repository root, use a fresh output directory:

    python docs/preprint_20261021/build_evidence.py --out /tmp/meta-acdc-paper-report
    python docs/preprint_20261021/verify_atomic_inputs.py --bundle docs/preprint_20261021/atomic_bundle_v1 --out /tmp/meta-acdc-atom-check.json

No private training data are needed to regenerate the paper's tables and figures. The Python standard library allows verification directly from the coordinate archive. The original Torch fit environment is recorded in the historical evidence. A complete historical neural replay additionally requires inputs not all included here; that limit is explicit.

## Submission boundary

The scientific deliverable is an author-review draft. The human author must critically review and shape its interpretations; funding, conflicts, correspondence and real affiliation remain to be supplied. No author name or institutional relationship is invented.

The official openRxiv statement dated20 October2026 requires verifiable human authors and an organizational affiliation for research-integrity accountability: https://openrxiv.org/approach-to-ai . Whether an unaffiliated researcher can satisfy submission eligibility must be clarified before submission. This requirement was discovered during preparation and is not hidden by entering an invented affiliation.

bioRxiv screening and posting are not guaranteed. No account was accessed, no message sent and no manuscript submitted in preparing this package. Future experimental collaborations and broader biological claims remain separate research work.

Latest user-requested pushes6d79c0e/dff57e3 were checked in latest_push_review.json. Report inputs match; the historical displacement result is not generalized to the new direction-aware helper.
