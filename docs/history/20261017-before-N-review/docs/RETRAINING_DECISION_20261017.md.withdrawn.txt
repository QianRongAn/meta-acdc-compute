# 2026-10-17 scoped retraining decision

**GO: structural native-versus-synthetic-decoy retraining. NO GO: experimentally measured binding or clinical-risk training.**

The 2026-10-16 repair release remains unchanged. New release `data/releases/20261017-structural-v4/` links the same validated arrays, evidence tables and frozen splits; it adds a separate, explicit scope decision and manifest. The public loader accepts the new manifest with its usual training checks enabled. No hash check was disabled and no old manifest was edited. hy4 remains deferred at the user's request; this is machine-reviewed readiness, not independent scientific validation.

On 2026-10-17, the current code/source fingerprints, all 888 graphs, QC, 86 families, donor isolation, external reservation and 199,279 unlabeled VDJdb associations were rechecked successfully. The actual six-layer, hidden-128 model was initialized from scratch and run on training-fold native/displaced/graft examples. Loss and gradients were finite; the checked graph's evaluation logit was identical alone and with other graphs (maximum absolute difference 0). No optimizer step, validation/test model evaluation or checkpoint write occurred. This verifies implementation readiness, not predictive performance or stability after training.

See `evidence/retraining_readiness_20261017.json`, `evidence/structural_training_release_20261017.json` and `evidence/structural_v4_release_manifest.json` for the checks and release fingerprint. The first report describes v3 before the separate v4 decision; its `training_manifest_allowed=false` is the preserved v3 state, superseded for structural training by the v4 decision only.

Run requirements: fresh initialization, CPU on this host, new output directories, optimizer on train only, model selection on validation, calibration/test kept held out, external primary receptor families excluded. Metrics describe discrimination of the specified synthetic controls, never binding probability or clinical accuracy. Preserve failed seeds and negative findings.

The available 2 GB GPU is not a requirement for CPU EGNN retraining. A structure-prediction backend is required for generating the missing consistent external SPR-pair structures; it is not required to retrain on the already rebuilt structural dataset. The external panel remains too small and case-selected to validate broad binding or clinical generalization. Those gaps do not prohibit the narrower structural experiment.

No training job has been started by this readiness assessment.
