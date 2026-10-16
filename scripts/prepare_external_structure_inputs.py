"""Prepare all SPR pairs with the same five-chain protocol; do not submit jobs."""
import argparse,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from meta_acdc.data.release import sha256
from meta_acdc.data.experimental_evidence import AA

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    rows=json.loads((ROOT/'docs/evidence/primary_measurements.json').read_text());jobs=[]
    for row in rows:
        if row['endpoint']!='tcr_binding_spr':continue
        tcr=[ROOT/x['path'] for x in row['sequence_sources']]
        pdb=tcr[0].name.split('_')[0]
        folder=ROOT/'data/external_evidence/20261016-expansion'
        mhc=folder/(pdb+'_entity_1.json');b2m=folder/(pdb+('_entity_3.json' if pdb=='8YIV' else '_entity_2.json'))
        paths=tcr+[mhc,b2m];sources=[json.loads(x.read_text()) for x in paths]
        if 'MHC' not in sources[2]['rcsb_polymer_entity']['pdbx_description'] and 'histocompatibility' not in sources[2]['rcsb_polymer_entity']['pdbx_description']:raise ValueError('MHC identity')
        if 'microglobulin' not in sources[3]['rcsb_polymer_entity']['pdbx_description']:raise ValueError('B2M identity')
        chains=[{'role':role,'sequence':x['entity_poly']['pdbx_seq_one_letter_code_can']} for role,x in zip(['alpha','beta','mhc','b2m'],sources)]
        chains.insert(2,{'role':'peptide','sequence':row['peptide']})
        if any(not set(x['sequence'])<=AA for x in chains):raise ValueError('noncanonical chain')
        fingerprint=hashlib.sha256(json.dumps(chains,sort_keys=True).encode()).hexdigest();jid='external_'+fingerprint[:16]
        fasta=''.join('>'+x['role']+'\n'+x['sequence']+'\n' for x in chains);(a.out/(jid+'.fasta')).write_text(fasta)
        jobs.append({'job_id':jid,'evidence_id':row['evidence_id'],'lineage':row['lineage'],'mhc_allele':row['mhc_allele'],'input_sha256':fingerprint,'fasta_sha256':sha256(a.out/(jid+'.fasta')),'chains':chains,'chain_sources':[{'path':str(x.relative_to(ROOT)),'sha256':sha256(x)} for x in paths],'state':'prepared_not_submitted'})
    if len({x['job_id'] for x in jobs})!=len(jobs):raise ValueError('duplicate pair input')
    (a.out/'jobs.json').write_text(json.dumps({'jobs':jobs,'n_jobs':len(jobs),'submitted':False,'protocol':'Use one frozen prediction backend/version and common seeds for every pair, same template/MSA settings, QC before scoring. Native positives and predicted negatives must not be compared as performance.','backend':None,'seeds':None,'execution_blocker':'No configured verified local five-chain structure prediction backend in this task; inputs prepared without submission.','training_allowed':False},indent=2)+'\n')
    print('Prepared',len(jobs),'five-chain inputs; no jobs submitted')
if __name__=='__main__':main()
