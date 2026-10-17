"""Shared single-chain CDR3 candidate associations; not paired-receptor validation.
Confidence strata describe database evidence, never binding labels.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

from meta_acdc.structure.graph import classify_chains, parse_pdb, TCR, PEPTIDE

CDR3_MIN_LEN = 8


def chain_sequences(path: Path) -> dict[str, str]:
    """Full ATOM sequence per chain.

    Appends CA residues in FILE ORDER — crucial for CDR3 regions, where
    insertion codes (resid 100A/100B) mean consecutive residues can share a
    number. Only consecutive duplicate (resid, icode) pairs are skipped
    (alternative conformations).
    """
    last: dict[str, tuple[int, str]] = {}
    seq: dict[str, list[str]] = defaultdict(list)
    with open(path, errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21].strip()
            resid = int(line[22:26])
            icode = line[26].strip()
            if last.get(chain) == (resid, icode):
                continue
            last[chain] = (resid, icode)
            seq[chain].append(line[17:20].strip())
    return {c: "".join(s) for c, s in seq.items()}


def aa3to1(seq: str) -> str:
    table = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
        "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
        "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M",
    }
    out = []
    for i in range(0, len(seq) - 2, 3):
        out.append(table.get(seq[i:i + 3], "X"))
    return "".join(out)


def association_hits(records, rows):
    """Exact CDR3 substring candidates with chain, MHC, species and evidence retained."""
    index = {'alpha':defaultdict(list),'beta':defaultdict(list)}
    for row in rows:
        for kind,column in [('alpha','cdr3_alpha'),('beta','cdr3_beta')]:
            cdr = row.get(column,'')
            if len(cdr) >= CDR3_MIN_LEN:
                index[kind][cdr].append(row)
    for record in records:
        for kind in ['alpha','beta']:
            seq = record[kind]
            # Candidate association only: no paired receptor identity or assay result inferred.
            lengths = {len(cdr) for cdr in index[kind]}
            found = {seq[i:i+n] for n in lengths for i in range(max(0,len(seq)-n+1))}
            for cdr in sorted(found & index[kind].keys()):
                for row in index[kind][cdr]:
                    confidence = row.get('confidence','')
                    if confidence == '': confidence = row.get('score','')
                    yield dict(pdb=record['pdb'],complex_id=record['complex_id'],chain=record[kind+'_chain'],
                               gene='TRA' if kind=='alpha' else 'TRB',cdr3=cdr,
                               vdjdb_epitope=row['peptide'],pdb_peptide=record['peptide'],
                               match='same' if row['peptide']==record['peptide'] else 'different',
                               confidence=confidence,reference_id=row.get('reference_id',''),
                               mhc_allele=row.get('mhc_allele',''),organism=row.get('organism',''),
                               evidence_scope='shared_single_chain_CDR3_candidate_association')


def main() -> int:
    import json
    from collections import Counter
    from meta_acdc.structure.identity import class_i_units
    ap = argparse.ArgumentParser()
    ap.add_argument('--struct-dir',type=Path,default=Path('data/raw/structures'))
    ap.add_argument('--vdjdb',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args = ap.parse_args()
    if args.out.exists() or args.out.with_suffix('.summary.json').exists():
        raise ValueError('refusing to overwrite association evidence')
    records=[]
    with (args.struct_dir/'complexes.tsv').open() as fh:
        pdbs=sorted({r['pdb'].lower() for r in csv.DictReader(fh,delimiter='\t')})
    for pdb in pdbs:
        parts,_=class_i_units(parse_pdb(args.struct_dir/(pdb+'.pdb')))
        for i,(_,meta) in enumerate(parts):
            records.append(dict(meta,pdb=pdb,complex_id=pdb+'_'+str(i)))
    with args.vdjdb.open() as fh: rows=list(csv.DictReader(fh,delimiter='\t'))
    hits=list(association_hits(records,rows))
    columns=['pdb','complex_id','chain','gene','cdr3','vdjdb_epitope','pdb_peptide','match','confidence','reference_id','mhc_allele','organism','evidence_scope']
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('w',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=columns,delimiter='\t');writer.writeheader();writer.writerows(hits)
    summary=dict(scope='candidate_association_not_cross_reactivity_prevalence',n_hits=len(hits),
                 confidence_by_hit=dict(Counter(str(x['confidence']) for x in hits)),
                 different_epitope_confidence_by_hit=dict(Counter(str(x['confidence']) for x in hits if x['match']=='different')),
                 denominator='matched database association rows, not independent paired receptors')
    args.out.with_suffix('.summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2));return 0


if __name__=='__main__':raise SystemExit(main())
