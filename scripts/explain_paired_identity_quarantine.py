"""Describe terminal gaps without imputing residues or changing frozen admission."""
from pathlib import Path
import json
from collections import Counter
from meta_acdc.structure.paired_variable_identity import germline_references,fasta_records
ROOT=Path(__file__).resolve().parents[1]
def main():
    base=ROOT/'docs/fingerprint_paired_identity_20261024';src=base/'sources';registry=json.loads((base/'curation_v3/identity_registry.json').read_text());refs=germline_references(list(src.glob('Homo_sapiens_TR*.fasta'))+list(src.glob('Mus_musculus_TR*.fasta')));audit=json.loads((ROOT/'docs/fingerprint_hard_family_20261024/identity_audit/identity_audit.json').read_text());old={x['index']:x for x in audit['rows']};rows=[]
    for row in registry['rows']:
        if row['status']!='quarantined' or row['index'] not in old:continue
        for role in ['alpha','beta']:
            detail=row.get(role,{})
            if detail.get('boundary_status')!='quarantined' or not detail.get('domain_candidate'):continue
            ann=detail['domain_candidate'];chain=detail['chain'];pdb=old[row['index']]['pdb'];path=src/(pdb+'_'+chain+'.fasta')
            if not path.exists() or not path.read_text().startswith('>'):continue
            gene=fasta_records(path.read_text())[0][0].split('|')[-1];variable=ann['sequence'];hits=[]
            for (species,g),values in refs.items():
                if g!=gene:continue
                for ref in values:
                    aa=ref['amino_acids']
                    for offset in range(1,13):
                        if variable.startswith(aa[offset:offset+10]):hits.append(dict(species=species,reference_header=ref['header'],candidate_missing_N_residues=aa[:offset],n_candidate_missing=offset))
            rows.append(dict(index=row['index'],complex_id=row['complex_id'],role=role,gene=gene,annotated_variable_prefix=variable[:15],terminal_offset_matches=hits,interpretation='exact terminal-offset match supports sequence truncation relative to reference, but does not establish intentional engineering or authorize adding residues' if hits else 'terminal offset not established; allele, terminal engineering, annotation and construct history remain unresolved'))
    counts=Counter('reference_terminal_offset' if row['terminal_offset_matches'] else 'other_unresolved_boundary' for row in rows)
    (base/'TERMINAL_QUARANTINE_EXPLANATIONS.json').write_text(json.dumps(dict(scope='post-freeze source-only explanation; labels and training protocol unchanged',counts=dict(counts),rows=rows),indent=2)+'\n');print(dict(counts))
if __name__=='__main__':main()
