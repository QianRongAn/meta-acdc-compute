"""Derive a new release with external families/donors removed; never reshuffle folds."""
from __future__ import annotations
import argparse,csv,json,os
from pathlib import Path
from .release import sha256
from .primary_panel import reservation_exclusions
from .experimental_evidence import filter_associations,assert_reserved_absent


def derive(source,output,registry_path):
    import torch
    # A valid historical manifest is necessary, even though current code/registry has evolved.
    original=json.loads((source/'release_manifest.json').read_text())
    for name,item in original['files'].items():
        if sha256(source/name)!=item['sha256']:raise ValueError('historical release changed: '+name)
    output.mkdir(parents=True,exist_ok=False)
    structures=json.loads((source/'structure_records.json').read_text())
    data=torch.load(source/'graphs.pt',map_location='cpu',weights_only=False)
    registry=json.loads(registry_path.read_text())['receptors']
    keep,excluded,complexes=reservation_exclusions(structures,data['records'],registry)
    remaining=[x for x in structures if x['complex_id'] not in complexes]
    assert_reserved_absent(remaining,registry_path)
    for key in ['graphs','labels','pdbs','peptide_masks','node_keys','records']:data[key]=[data[key][i] for i in keep]
    data['provenance']=dict(data['provenance'],reservation_source_manifest_sha256=sha256(source/'release_manifest.json'),reservation_registry_sha256=sha256(registry_path))
    torch.save(data,output/'graphs.pt')
    (output/'structure_records.json').write_text(json.dumps(remaining,indent=2))
    (output/'graph_records.json').write_text(json.dumps(data['records'],indent=2))
    splits={fold:[i for i,x in enumerate(data['records']) if x['fold']==fold] for fold in ['train','val','calibration','test']}
    # Use exact existing fold names; no reassignment or regeneration of controls.
    splits={fold:[i for i,x in enumerate(data['records']) if x['fold']==fold] for fold in json.loads((source/'splits.json').read_text())}
    if any(not indices or {data['labels'][i] for i in indices}!={0,1} for indices in splits.values()):raise ValueError('empty or one-class fold after reservation')
    (output/'splits.json').write_text(json.dumps(splits,indent=2))
    with (source/'family_map.tsv').open(newline='') as f,(output/'family_map.tsv').open('x',newline='') as out:
        reader=csv.DictReader(f,delimiter='\t');w=csv.DictWriter(out,fieldnames=reader.fieldnames,delimiter='\t');w.writeheader();w.writerows(x for x in reader if x['complex_id'] not in complexes)
    for name in ['iedb.assays.tsv','iedb.assays.tsv.rejected.tsv','clinical_index_cases.tsv','native_qc.tsv','structure_qc.tsv','input_sources.json']:
        if (source/name).exists():os.link(source/name,output/name)
    assoc=filter_associations(source/'vdjdb.associations.tsv',output/'vdjdb.associations.tsv',output/'vdjdb.reserved_quarantine.tsv',registry_path)
    # Verify the retained associations rather than trusting the filter's counters.
    with (output/'vdjdb.associations.tsv').open(newline='') as f:assert_reserved_absent(csv.DictReader(f,delimiter='\t'),registry_path)
    report={'source_release':str(source),'source_manifest_sha256':sha256(source/'release_manifest.json'),'registry_sha256':sha256(registry_path),'excluded_complexes':sorted(complexes),'excluded_graphs':excluded,'n_graphs_before':len(keep)+len(excluded),'n_graphs_after':len(keep),'n_native_after':len(remaining),'fold_counts':{k:len(v) for k,v in splits.items()},'folds_preserved':True,'association_filter':assoc,'training_allowed':False,'reason':'external receptor reservation only; arrays retained byte-for-byte, old release unchanged'}
    (output/'reservation_report.json').write_text(json.dumps(report,indent=2))
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--registry',type=Path,required=True);a=p.parse_args();print(json.dumps(derive(a.source,a.out,a.registry),indent=2))
if __name__=='__main__':main()
