> 当前证据状态：[SPR第三轮续查、不同V基因候选与KRAS工程参考](DIVERSITY_FOLLOWUP_20261018.md)。严格池31条（原分区16/12/3，测试仅2块）；外部核心7条保持未评分。新增2个不同V基因的配对CDR3候选，完整可变区待补；KRAS10mer补到599aa公开工程蛋白参考，未冒充临床DNA。仍有10条原文条件不足、2条精确表值和2条温度缺口，训练继续阻断。

> 历史版本记录（已由上方续查更新）：[临床原文、工程构型和结局整理](CLINICAL_PRIMARY_CURATION_20261018.md)。临床资料仍按患者/组合保留；新SPR开发集51条，临床候选队列独立预留，旧版本数据量和结果仅作历史。

> **后续原子数据处置已落地，当前判决见 [N2／N8 修复报告](STRUCTURE_ATOMIC_REPAIR_20261017.md)：168 份无标签资料核验，5 对新控制不足以恢复重训。**

> **2026-10-17 update:** [Scoped structural retraining is now released](RETRAINING_DECISION_20261017.md); binding and clinical claims remain blocked. The release/guard state described below is historical.

> **最新状态：** [证据扩展与 structural-v3](EXPERIMENTAL_EVIDENCE_EXPANSION.md) 已完成机器验证；hy4 按用户要求暂缓。旧报告中的配额阻断和 v2 计数为历史状态；训练尚未放行。

> Historical repair record. Current experimental-data status and the user’s hy4 waiver are documented in [EXPERIMENTAL_EVIDENCE.md](EXPERIMENTAL_EVIDENCE.md). Existing manifest hashes are preserved; new guard code requires a new manifest before training.

# Pretraining readiness repair — 2026-10-16

Decision: machine validation passed for a structural-decoy prototype; **training remains blocked pending independent hy4 acceptance**. The provider returned HTTP 402, `Insufficient account funds`, before any review. No training, scheduling, replacement of historical artifacts, or restored performance claims occurred.

## Repaired blockers

| Blocker | Repair and evidence |
|---|---|
| Old VDJdb labels still loadable | Rebuilt 199,385 association rows with empty label and explicit confidence/unlabeled scope. Sequence baseline and active-learning entry points reject these as binding labels. ProtoNet's unsupported association-to-negative training path is explicitly disabled. |
| Actual IEDB CSV incompatible with the TSV importer | Detect and validate the two-row header before output; stream CSV/TSV, preserve assay IDs, raw values/units/qualifiers and conflicting calls at assay level. Exact endpoint and unit requirements prevent unspecified binding/elution evidence becoming IC50 or TCR labels. Full 9,207,473,993-byte input rebuilt to 4,845,721 assay records. |
| Wrong chain identities / incomplete or ambiguous complexes | Depositor COMPND identities plus unique heavy-atom contacts assign class-I alpha/beta complexes; ambiguous, incomplete and other-scope inputs are recorded in exclusions. Peptide canonical heavy atoms must be complete. No chain-letter guessing in the release builder. |
| Defective graphs, stale masks, identical opposite-label examples | Regenerated current heavy-atom-centroid graphs; explicit peptide mask, canonical amino-acid encoding, same-length complete-backbone grafts. Deduplicate model inputs across every class. Donors belong to the same frozen fold. |
| Crystal/AF3 and receptor leakage | Group same PDB and near paired variable-domain sequences before generating decoys; predicted copies inherit reference family/fold. Four frozen folds keep model fitting, validation, temperature fitting and final test separate. Family rule is paired-chain 90% sequence similarity and 90% CDR3 similarity when detected, using connected components. This operational definition is not proof that every biological homolog has been excluded. |
| Job-level QC borrowed model-0 confidence | Rebuilt QC for 2,109 individual files; 257 meet identity/parse/integrity/confidence filters. Missing confidence and failed parsing are failures, regardless of another model's confidence. Shared confidence files apply only to model 0. QC retains source/reference/confidence hashes. ipTM>=0.75 is a project filter, not binding evidence. |
| Clinical table stale despite source changes | Regenerated the index; removed structural links for unverified receptor identities. The table is an evidence index and cannot validate clinical risk. |
| Training / calibration / scoring differ | Per-graph normalization, separate peptide mask, graph-specific coordinate normalization. Checkpoints carry release hash and confidence-masking policy; calibration/scoring consume that policy. Calibration and test folds are separate from model fitting. Family bootstrap uses the release family ledger rather than PDB names. |
| Legacy data/checkpoints remain selectable | Structural trainers require an accepted manifest and frozen folds, hash-check sources/artifacts/code, reject legacy/mismatched inputs, ad hoc positive directories and overwriting checkpoints. A manifest with pending independent acceptance cannot start training. |

Good ideas adopted from the sibling project: separate evidence status from labels, explicit chain roles and peptide masks, provenance and failure ledgers, and matched structural controls. Its simulated labels and performance claims were not imported.

## Rebuilt local release

`data/releases/20261016-structural-v2/` (large data are local and git-ignored):

- 910 graphs: 306 native, 306 displaced and 298 graft controls.
- 237 crystal native units + 69 QC-accepted predicted natives; 88 receptor families.
- Train 559, validation 157, calibration 67, test 127 graphs, both structural classes in each fold.
- VDJdb 199,385 unlabeled associations. IEDB 306,616 binding-method-unspecified and 4,539,105 eluted-ligand records; these are MHC evidence, not experimental TCR binding labels.
- Independent array/source/fold validator passed all 910 graphs: finite tensors, feature identities, peptide masks, edge/RBF consistency, input uniqueness, donor-fold isolation, source hashes, model QC, family/PDB split isolation and evidence-table scopes.
- Unit suite: 40 tests, 37 passed and 3 skipped for missing legacy fixtures. Skipped: attribution with an old checkpoint, old CIF fixture, and legacy electrostatics fixture. New tests include single-versus-mixed-batch equality, translation invariance and training-versus-scoring equality using random weights; no optimizer run.
- The previous v1 draft is retained as history. v2 supersedes it after peptide atom completeness checks excluded additional sources. No release overwrote historical caches, model files or raw data.

Machine validation means reproducible input integrity for this operational structural task. Grafts/displacements are synthetic controls without experimental nonbinding evidence; graft relaxation and artifact sensitivity remain methodological limitations. These data cannot authorize binding AUROC, recall, cross-reactivity accuracy or clinical-risk probability. No historical performance number is restored.

## Reproduction and gate

Use the existing environment; no reinstall:

```bash
cd /home/administrator/docs/meta-acdc-readiness
export PYTHONPATH="$PWD/src" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
PY=/home/administrator/docs/meta-acdc/.venv/bin/python
$PY -m unittest discover -s tests -v
$PY -m meta_acdc.data.validate_release --release-dir data/releases/20261016-structural-v2
$PY scripts/pretrain_audit.py --project-root "$PWD" \
  --release data/releases/20261016-structural-v2/release_manifest.json \
  --report /tmp/meta-acdc-readiness.json
```

The first validator checks data only. The final gate reports `WAIT_INDEPENDENT_ACCEPTANCE` and exits 2 until an independent verdict is obtained; this is intentional. Trainers must refuse the current `training_allowed:false` manifest. Never edit it in place to bypass acceptance. After a successful independent review of the exact code and artifact hashes, issue a new approved manifest version referencing the same immutable data; retain the pending manifest and failure log. A later source, graph, QC or fold change requires revalidation and another review.

To rebuild, choose a new empty release directory and use `meta_acdc.data.build_release --structures /home/administrator/docs/meta-acdc/data/raw/structures --predicted-root /home/administrator/docs/meta-acdc/data/raw/af3_native_predictions --out NEW_RELEASE_DIR`; regenerate evidence tables with `data.ingest_vdjdb --in ... --out ...`, `data.ingest_iedb --in ... --out ...` and `data.clinical_cases --out ...`. Source paths and hashes are in input_sources.json and structure_records.json. The validation/publishing command requires every table and rejects an existing graph or published manifest rather than replacing them.

## Remaining release conditions

1. Independent hy4 review and acceptance of this exact revision. Model ID reused: `opencode-go/hy4-preview`. Session: `ses_ebbb117e5ffefxzj5cM2zKL1dj`; provider quota HTTP 402. No verdict exists. Do not retry automatically or duplicate work while a task is running.
2. For a future binding/activation task: traceable assay-specific TCR–peptide–MHC positives/negatives, conflict handling and independent family-isolated evaluation. Current association/MHC/clinical-index data do not provide these labels. `binding_training_allowed` stays false.
3. Structure scoring and biological conclusions require new compatible checkpoints and evaluation after approval; old models remain historical. The present changes contain no new training result.

The original worktree is clean; six audited historical artifact hashes are unchanged (preservation.json). Branch: `audit/pretrain-20261016`. No merge to main is part of this repair.
