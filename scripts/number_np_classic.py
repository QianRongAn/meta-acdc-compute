"""Independent classic HMM numbering of unchanged NP raw constructs; no admission."""
from pathlib import Path
import argparse,hashlib,json,importlib.metadata
import anarci
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'docs/fingerprint_np_boundaries_20261029'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--hmmerpath',required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    raw=json.loads((BASE/'numbering/numbering_raw.json').read_text());inputs=raw['inputs'];sequences=[(r['name'],r['sequence']) for r in inputs]
    numbered,details,hits=anarci.anarci(sequences,scheme='imgt',allow={'A','B'},hmmerpath=a.hmmerpath,ncpu=2,assign_germline=True,allowed_species=['mouse'],bit_score_threshold=80)
    package=Path(anarci.__file__).parent
    result=dict(inputs=inputs,numbered=numbered,details=details,hits=hits,provenance=dict(package_version=importlib.metadata.version('anarci'),biopython=importlib.metadata.version('biopython'),scheme='imgt',bit_score_threshold=80,assign_germline=True,allowed_species=['mouse'],hmmer_version='3.4',package_files={str(p.relative_to(package)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(package.rglob('*')) if p.is_file() and p.suffix in ['.py','.hmm','.h3f','.h3i','.h3m','.h3p']}),identity_admitted=False,protocol_frozen=False,model_training_started=False,independent_acceptance_accessed=False)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print([(r['name'],None if not n else (n[0][1],n[0][2],len(n[0][0])),d) for r,n,d in zip(inputs,numbered,details)])
if __name__=='__main__':main()
