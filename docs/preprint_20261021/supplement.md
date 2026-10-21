# Supplementary information

Companion to: Decoy artifacts and limited source-disjoint affinity prediction in a TCR-pMHC computational workflow.

## S1. Scope and chronology

This is a retrospective, single-project analysis. Evidence was known before the manuscript analysis plan. The plan adds fixed descriptive aggregation and public coordinate replication; it is not described as a prospectively registered hypothesis test. The underlying affinity experiment was planned and model-locked before its test was opened. The bounded development experiment used a separate plan and did not reevaluate the consumed test.

Historical chronology is retained in Git: 0b8e6acf3fc56a628e6434f6985022e959117603 fixed the affinity protocol; c0c2095762950b212001dffd84600c6ec59d5083 locked the model before test access; b59b5da6240f0583cb5832407e5c2a8e1b804bf0 retained the completed acceptance. Later development plans were recorded at 7f8e9b1 and d41af3d before fitting the bounded alternatives. The present report uses the e4d68de core evidence snapshot, with 15eb5c6 contextual reports. User-identified updates 6d79c0e and dff57e3 were checked: all nineteen report inputs are unchanged. The main branch also contains a direction-aware displacement repair; the old fixed-axis overlap rate is not applied to that revised generator. These are version-controlled project records, not an external clinical-trial registration.

The core manuscript builder reads fifteen enumerated report/protocol JSON files, with four separately hashed later contextual reports. It never imports a training or scoring entry point, decodes individual test predictions, or accesses the fresh reserve. No threshold or model is selected for a better manuscript outcome. The previous test remains consumed and the fresh reserve remains unavailable for tuning.

## S2. Structural construction evidence

Table S1. Fixed severe-overlap screen. Percentages describe records in this collection. Audit families are inherited groups and are not biological replicates.

| Kind | Records | Audit families | Overlap flags | Record rate | Family mean | Incomplete interface |
| --- | --- | --- | --- | --- | --- | --- |
| native | 300 | 86 | 0 | 0.00% | 0.00% | 9 |
| graft | 288 | 83 | 146 | 50.69% | 57.99% | 9 |
| displaced | 300 | 86 | 300 | 100.00% | 100.00% | 9 |

The native/displacement sets span 86 inherited families and the graft set 83. Repetition within a family is retained. Equal-family weighting means averaging the within-family flag fractions; it does not make families representative of a population. Rows and all 588 parent-control mappings are supplied in structural_qc.csv, family_qc.csv and paired_qc.csv.

Table S2. Fixed QC-only rule on paired native/control records. A record is called synthetic if it contains a flagged severe cross-chain distance. Balanced accuracy is (sensitivity + specificity)/2. There is no fitted threshold, held-out construction test, biological class assignment, or inference of a neural model's learned mechanism.

| Control | Parent-control pairs | Sensitivity | Specificity | Balanced accuracy |
| --- | --- | --- | --- | --- |
| graft | 288 | 0.506944 | 1.000000 | 0.753472 |
| displaced | 300 | 1.000000 | 1.000000 | 1.000000 |

The public atomic_inputs.json.gz archive contains the 300 reconstructed source complexes and peptide coordinates for every historical native/control record. Each source is linked to a source-file SHA256, basename, domain and PDB code. Every graph retains the source complex ID, construction kind, inherited family, node keys and donor ID where applicable. Coordinates are serialized without additional decimal rounding. Source material is derivative PDB/prediction input, not new wet-laboratory evidence.

During export, the existing construction helper reproduced all 888 retained screen reports and their interface-missing-atom lists. The independent distance checker reads only this archive and the original audit JSON, uses a Python standard-library spatial cell-grid implementation, and confirms the count of every peptide-environment heavy-atom pair at <=1.2 A. It does not use torch models or download original structures. Agreement between two implementations is a reproducibility check; it does not validate the cutoff's biological meaning. The exporter reads only the project's own trusted pickle cache and is not suitable for arbitrary untrusted torch files.

## S3. Structural repair dispositions

The historical 888 records were not repurposed as binding labels. The 588 synthetic controls were retired; 168 native assets were retained without response labels; 132 native records were not admitted under the new policy. Their reasons included 97 contact-screen exclusions, nine incomplete interfaces, 25 predicted identity/reference-QC failures, and one changed historical-reference attribution. These categories describe policy decisions and do not prove that every excluded source structure is experimentally wrong. The new five construction pairs have no calibration/test pairs. No repaired structural-model accuracy is claimed.

## S4. Counting units and grouping sensitivity

Table S3. Actual panel counts. Full-construct identity and similarity are distinct from independently annotated variable-region or V-gene families.

| Unit | Count |
| --- | --- |
| FASTA files | 347 |
| PDB prefixes | 111 |
| Exact receptor pairs | 98 |
| Hamming components | 89 |
| Matching-block components | 57 |

All 347 FASTA hashes are retained in the pinned pending_panel.json, alongside sequences and the first committed version. There are 91 peptide strings, 85 alpha strings, 76 beta strings, 45 MHC strings and 11 B2M strings. Component counts can differ because changing the pairwise comparison creates different graph edges and transitive links. Neither 89 nor 57 is promoted as the correct biological family count, and the unavailable 73-cluster algorithm was not reconstructed by tuning these rules.

## S5. Affinity evaluation and bounded development

Table S4. Already-completed locked affinity experiment. All errors are in pKD. Test metrics use equal total source weight, equal lineage weight within each source, and equal observation weight within each lineage. The source groups are coauthor-connected components. No new per-observation predictions were made to construct this table.

| Subset | Method | Observations | Lineages | Sources | RMSE | MAE |
| --- | --- | --- | --- | --- | --- | --- |
| primary | Ridge | 14 | 7 | 3 | 0.996215 | 0.848597 |
| primary | Training mean | 14 | 7 | 3 | 1.003632 | 0.853082 |
| sensitivity | Ridge | 15 | 7 | 3 | 1.024121 | 0.876009 |
| sensitivity | Training mean | 15 | 7 | 3 | 1.031228 | 0.880495 |
| Baker-linked | Ridge | 5 | 4 | 1 | 1.116805 | 1.014953 |
| Baker-linked | Training mean | 5 | 4 | 1 | 1.127836 | 1.012459 |
| Immunocore | Ridge | 7 | 2 | 1 | 1.178506 | 0.992569 |
| Immunocore | Training mean | 7 | 2 | 1 | 1.179479 | 0.992258 |
| Wu et al. | Ridge | 2 | 1 | 1 | 0.584124 | 0.538268 |
| Wu et al. | Training mean | 2 | 1 | 1 | 0.598871 | 0.554530 |

Table S5. Primary studies underlying the external affinity panel. Variant names are not counted as independent receptor lineages. Full measurement-specific peptide, sequence, orientation, fit and source locators are preserved in the existing frozen records; this manuscript does not amend those records.

| PMID | First author | Year | DOI |
| --- | --- | --- | --- |
| 33468649 | Smith AR | 2021 | 10.1073/pnas.2018125118 |
| 35484264 | Chandran SS | 2022 | 10.1038/s41591-022-01786-3 |
| 39042676 | Gaglione SA | 2024 | 10.1073/pnas.2406615121 |
| 39500929 | Karuppiah V | 2024 | 10.1038/s41598-024-75140-7 |
| 39576860 | Wu D | 2024 | 10.1126/sciadv.adq6150 |

The coauthor graph joins the Smith, Chandran and Gaglione studies into one source component. Karuppiah and Wu form the other two components. Same-study, receptor-chain similarity and peptide similarity also define four connected test blocks. The use of source components is conservative; author-name ambiguity and incomplete provenance remain possible. Sharing a V gene alone is not an exclusion under the retained protocol, so this is not V-gene-disjoint testing.

The original training and validation sets each had ten observations and five connected blocks. The later bounded development experiment excluded ATLAS:485 from validation before candidate fitting because its source overlapped a new reserve. The remaining nine validation observations occupy four blocks. All original remaining folds and blocks were preserved, with no reshuffling to increase independent counts.

Table S6. All fixed development candidates on the nine-observation validation set. Alphas of new candidates were selected inside the five inherited training blocks only. The reference retains its previous configuration-selection history.

| Candidate | Alpha | Validation RMSE | Validation MAE | Prediction SD |
| --- | --- | --- | --- | --- |
| anchor_raw | 1.0 | 0.772873 | 0.599472 | 0.260809 |
| sequence_balanced | 100.0 | 0.793244 | 0.663391 | 0.026379 |
| cdr3_chemical | 0.1 | 0.917451 | 0.873417 | 0.278756 |
| cdr3_words | 0.1 | 0.907292 | 0.862141 | 0.272811 |
| Training mean | N/A | 0.802913 | 0.666416 | 0.000000 |

The 1,014-feature reference includes one/two-residue frequencies, length terms, composition-product interactions, MHC indicators, instrument indicators and temperature terms. The chemical variants and group normalization are fully defined in the versioned plan_v2.json and development_affinity.py. Their results do not establish a population ranking of representations. Unseen vocabulary entries are ignored at application; the small sample and its assay composition leave substantial uncertainty about transferability.

## S6. Historical score provenance

The five-seed reconstruction covered 188 jobs, 884 AF3 structures, 940 job-by-seed means and 4,420 structure-by-seed scores. The retained four-decimal matches were 938/940 and 4,417/4,420. All discrepancies are 0.0001 at the reported precision, with causes not yet isolated. Fourteen jobs have only one retained AF3 input; they do not provide five-model variability estimates. The separate 21-job scoring table cannot yet be bound to its original checkpoint.

The audit's historical numeric results remain marked FORENSIC_WITHDRAWN. Present-day reconstruction is not evidence of the exact historical command or file overwrite order. No historical clinical ranking is reinstated. MAG-IC3 and a3a are distinct engineered constructs, and near scores are not a statistical equivalence test.

## S7. Evidence tiers and availability

Tier A: all manuscript tables and figures can be regenerated from tracked, hash-pinned reports. Tier B: all severe-distance counts can be recomputed from the tracked 28.3-MB coordinate archive, using a second algorithm. Tier C: full historical neural-score replay additionally needs old checkpoints and 884 input files; the hashes and local replay results are tracked, but the complete Tier C inputs are not included in this manuscript package. No claim of complete public Tier C replication is made.

The old learned model and clinical interpretation are not required to substantiate the narrow Tier A/B results. Full coordinate-to-source identity can be checked against the recorded original-file hashes when those inputs are available. That chain does not remove sequence/role ambiguities already disclosed in the historical audit.

The project metadata declares MIT for code. Original database entries, structures and papers retain their original terms. Metadata access is not a blanket redistribution license. No third-party article PDF is bundled. The final preprint license is a separate author decision.

## S8. Reproduction commands

From the repository root, choose fresh output paths. Existing result directories are deliberately not overwritten.

    python docs/preprint_20261021/build_evidence.py --out /tmp/meta-acdc-paper-report

    python docs/preprint_20261021/verify_atomic_inputs.py --bundle docs/preprint_20261021/atomic_bundle_v1 --out /tmp/meta-acdc-atom-check.json

    PYTHONPATH=src python -m unittest discover -s tests -p test_preprint_evidence.py -v

The report builder uses Python's standard library and Matplotlib. Coordinate verification uses only the Python standard library. The existing historical fit environment is recorded in its fit_report.json; no training should be restarted with the consumed test profile. The PDF renderer additionally uses ReportLab and records its dependencies. Packages were reused from existing environments, with no installation in this round.

## S9. Post-snapshot contextual evidence from a separate endpoint

The four hashes in context_manifest.json bind the subsequent 15eb5c6 project update. These reports were completed before manuscript inclusion. No additional model was run or selected to improve the manuscript, and the core fifteen-file report builder is unchanged.

Table S7. Association with the same previously consumed MIRA functional observations. This is not SPR regression, clinical validation, or a new independent holdout. ERGO's target and beta-CDR3/peptide inputs differ from those of the SPR-reference score.

| Score | Per-peptide macro AUROC | Per-peptide macro AP |
| --- | --- | --- |
| Official fixed ERGO-I McPAS/LSTM | 0.511247 | 0.116626 |
| Historical SPR-reference score | 0.500078 | 0.119691 |

The official code commit is 5c6fc3719df08c84be83a7f365009cf49b2f23b7 and the checkpoint hash is 238a6f44b15eb499013cd76a439cf82853bd601fc784cd375907c8cc65eb915b. A technical amendment, recorded before metric calculation, addressed upstream last-batch truncation by appending one declared repeat of the first real query. All 10,001 returned pair identities and their order were checked; the repeated score agreed exactly. Only the declared final control was excluded, leaving all 10,000 scientific observations. The official code and weights were not modified. Exact beta-CDR3 overlaps affected 120 observations, with zero exact beta-CDR3/peptide pair matches; full alpha/beta family and source independence remains unestablished.

The conditional TCR4 evidence is traceable to Ma et al., Figure 1c/d, Methods, Source Data and 7L1D. Conditions were Biacore T200, 4 C, TCR immobilized and 1:1 steady-state fitting. Neoantigen KD was 62 +/- 6 microM across six reported affinity replicates. Wild-type non-detection at 208.24 microM was retained with a sample-integrity control; no quantified detection limit supports a numeric KD bound. The exposure audit identifies the previously evaluated TCR4 lineage and shared Baker provenance. New independent-lineage count is zero. Another paper's different receptor named TCR4 is not merged by name.

The newly indexed BATCAVE and P-TEAM candidate collections have not completed admission under this manuscript's affinity protocol. They are not added to training counts, are not treated as independent tests, and are not converted from functional silence to physical nonbinding. Their curation remains a separate project task.
