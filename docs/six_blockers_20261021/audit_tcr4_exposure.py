from pathlib import Path
import json,sys
sys.path.insert(0,'/home/administrator/docs/meta-acdc-readiness/src')
from meta_acdc.data.build_release import domain
from meta_acdc.data.experimental_evidence import reserved_hits
from meta_acdc.data.independent_spr import author_keys
R=Path('/home/administrator/docs/meta-acdc-readiness');D=R/'docs/six_blockers_20261021';raw=R/'data/external_evidence/20261021-six-blockers'
row=json.loads((D/'tcr4_primary_challenge.json').read_text())['observations'][0]
query=dict(alpha_variable=domain(row['alpha_pdb_construct'])[0],beta_variable=domain(row['beta_pdb_construct'])[0],cdr3_alpha=row['cdr3_alpha'],cdr3_beta=row['cdr3_beta'])
registries={}
for name,path in [('old91','docs/database_binding/affinity_dataset.json'),('clinical','docs/clinical_primary_20261018/clinical_reservations.json'),('combined_prior_reservations','docs/database_binding/combined_reservation_registry.json')]:
 obj=json.loads((R/path).read_text());rows=obj if isinstance(obj,list) else obj['receptors'];registries[name]=reserved_hits(query,rows)
# Project identity fields only; no outcomes or old predictions enter this audit.
frozen=json.loads((R/'docs/independent_spr_20261018/frozen_test.json').read_text())
identity=[{k:x.get(k,'') for k in ['receptor_name','alpha_variable','beta_variable','cdr3_alpha','cdr3_beta']} for x in frozen]
registries['consumed_frozen_test_identity_only']=reserved_hits(query,identity)
study=json.loads((raw/'pubmed_39833157_core.json').read_text())['resultList']['result'][0];key=author_keys(dict(authors=study['authorList']))
source_hits=[dict(pmid=s['pmid'],shared_author_keys=sorted(key&author_keys(s))) for s in json.loads((R/'docs/independent_spr_20261018/development_study_sources.json').read_text()) if key&author_keys(s)]
out=dict(construct_id=row['construct_id'],identity_fields=query,registry_hits=registries,development_source_overlap=source_hits,previous_measurement_pmid='35484264',new_measurement_pmid='39833157',new_receptor_lineages=0,new_fresh_independent_acceptance_records=0,old_consumed_test_identity_projected_for_exposure_audit=True,old_consumed_test_labels_used_for_tuning=False,new_measurements_model_scored=False,training_allowed=False,formal_independent_acceptance_allowed=False,reason='Existing TCR4 lineage is in consumed independent registry. New SPR non-detection/4 C measurement extends evidence only; never count same receptor from new paper as a fresh independent family.')
assert registries['consumed_frozen_test_identity_only']
(D/'tcr4_exposure_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
