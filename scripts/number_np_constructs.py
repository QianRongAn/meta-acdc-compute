"""Number raw NP deposited constructs with pinned ANARCII; annotation only, no identity admission."""
from pathlib import Path
import argparse,hashlib,json,importlib.metadata
from anarcii import Anarcii
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_boundaries_20261029'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    facts=json.loads((BASE/'final_primary_cells.json').read_text())
    entries={e['rcsb_id'].lower():e for p in (ROOT/'docs/fingerprint_expanded_controls_20261029/sources/metadata_seed').glob('batch*.json') for e in json.loads(p.read_text())['data']['entries'] if e}
    inputs=[]
    for r in facts['rows']:
        for role,chain in r['chains'].items():
            seq=chain['primary_labeled_variable'];matches=[e for e in entries[r['pdb']]['polymer_entities'] if e['entity_poly']['pdbx_seq_one_letter_code_can'].count(seq)==1];assert len(matches)==1
            full=matches[0]['entity_poly']['pdbx_seq_one_letter_code_can'];inputs.append(dict(name=r['name']+'_'+role,pdb=r['pdb'],role=role,entity=matches[0]['rcsb_id'],author_chains=matches[0]['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'],sequence=full,sequence_sha256=hashlib.sha256(full.encode()).hexdigest()))
    fasta=''.join('>'+r['name']+'\n'+r['sequence']+'\n' for r in inputs);(a.out/'inputs.fasta').write_text(fasta)
    model=Anarcii(seq_type='tcr',mode='accuracy',cpu=True,ncpu=4,batch_size=8)
    result=model.number(str(a.out/'inputs.fasta'))
    package=Path(__import__('anarcii').__file__).parent
    provenance=dict(package_version=importlib.metadata.version('anarcii'),parameters=dict(seq_type='tcr',mode='accuracy',cpu=True,ncpu=4,batch_size=8),package_files={str(p.relative_to(package)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(package.rglob('*')) if p.is_file() and p.suffix in ['.pt','.pth','.py','.json']},scheme='IMGT',purpose='independent annotation inference only; not STCRDab evidence, not training or admission')
    (a.out/'numbering_raw.json').write_text(json.dumps(dict(inputs=inputs,results=result,provenance=provenance,identity_admitted=False,model_training_started=False,independent_acceptance_accessed=False),indent=2)+'\n')
    print({k:dict(chain=v.get('chain_type'),start=v.get('query_start'),end=v.get('query_end'),error=v.get('error')) for k,v in result.items()})
if __name__=='__main__':main()
