"""Observed ATOM sequence hashes are not canonical paired-variable receptor truth."""
from collections import defaultdict
from pathlib import Path
from meta_acdc.structure.graph import AA3TO1

IDENTITY_SEMANTICS='observed_alpha_beta_sequence_hash; canonical_paired_variable_identity_unverified'

def seqres_sequences(text):
    seq=defaultdict(str);declared={}
    for line in text.splitlines():
        if not line.startswith('SEQRES'):continue
        chain=line[11];declared[chain]=int(line[13:17])
        seq[chain]+=''.join(AA3TO1.get(x,'X') for x in line[19:70].split())
    return {c:s for c,s in seq.items() if len(s)==declared[c] and 'X' not in s}

def is_subsequence(observed,reference):
    it=iter(reference)
    return all(any(x==a for x in it) for a in observed)

def same_construct_conflicts(rows):
    """Only exact complete deposited alpha AND beta sequences may establish equivalence."""
    groups=defaultdict(list)
    for row in rows:
        if row['eligible_construct_evidence']:
            groups[(row['seqres_alpha'],row['seqres_beta'])].append(row)
    return [v for v in groups.values() if len({r['observed_sequence_id'] for r in v})>1]
