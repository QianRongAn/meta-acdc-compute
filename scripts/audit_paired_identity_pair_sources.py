from pathlib import Path
import json,torch
r=Path(__file__).resolve().parents[1];base=r/'docs/fingerprint_paired_identity_20261024';p=json.loads((base/'development/protocol.json').read_text());d=torch.load(r/p['dataset'],map_location='cpu',weights_only=False);rows=[]
for fold in p['folds']:
 positive=[];negative=[]
 for i in fold['validation_indices']:
  for j in fold['validation_indices']:
   if i==j:continue
   a,b=d['records'][i],d['records'][j];entry=dict(query=i,gallery=j,same_pdb=a['pdb']==b['pdb'])
   (positive if a['receptor_id']==b['receptor_id'] else negative).append(entry)
 rows.append(dict(fold_id=fold['fold_id'],n_positive_directed_pairs=len(positive),n_negative_directed_pairs=len(negative),all_positive_same_pdb=all(x['same_pdb'] for x in positive),all_negative_different_pdb=all(not x['same_pdb'] for x in negative),positive_pairs=positive,negative_pairs=negative))
assert all(x['all_positive_same_pdb'] and x['all_negative_different_pdb'] for x in rows)
(base/'pair_source_audit.json').write_text(json.dumps(dict(scope='metadata diagnostic only, not a newly selected model or baseline',all_positive_pairs_same_pdb=True,all_negative_pairs_cross_pdb=True,biological_generalization_blocked=True,reason='current positives are repeat crystal copies, so PDB origin is fully confounded with pair identity in these held-out contexts; no claim of cross-structure receptor generalization',folds=rows),indent=2)+'\n')
print('PAIR_SOURCE_CONFOUND_CONFIRMED')
