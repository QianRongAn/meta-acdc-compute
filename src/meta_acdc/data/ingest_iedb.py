"""Stream IEDB's two-header CSV/TSV as an assay evidence ledger.

No peptide/allele aggregation: contradictory assays stay separate. MHC assays
are never TCR-binding labels. Output is atomically published after validation.
"""
from __future__ import annotations
import argparse, csv, math, re
from pathlib import Path
try:
    from .schema import DataSource, DatasetStats, InteractionRecord, SCHEMA_HEADER
except ImportError:
    from schema import DataSource, DatasetStats, InteractionRecord, SCHEMA_HEADER


def assay_family(method: str) -> str:
    m = method.lower()
    if any(k in m for k in ('koff','kon','half-life','t1/2','dissociation rate','stability')):
        return 'kinetics_or_stability'
    if re.search(r'\bkd\b|dissociation constant', m):
        return 'binding-kd'
    if 'elut' in m or 'mass spectr' in m:
        return 'presentation-eluted'
    if 'ic50' in m:
        return 'binding-ic50'
    return 'binding-unspecified' if any(k in m for k in ('binding','mhc')) else 'other'


def ic50_score(raw: str, units: str):
    # Values with '<'/'>' are censored, not exact affinities.
    try:
        v = float(raw)
    except (ValueError,TypeError):
        return None, 'non-exact quantitative value'
    u = units.strip().replace('μ','u').replace('µ','u').casefold()
    factors = {'nm':1.0, 'um':1000.0, 'pm':0.001, 'm':1e9}
    if not math.isfinite(v) or v <= 0 or u not in factors:
        return None, 'nonpositive/nonfinite value or unsupported unit'
    return -math.log10(v*factors[u]), 'score=-log10(value in nM); endpoint-specific'


def parse_score(raw):
    return ic50_score(raw,'nM')[0]


def ingest(path: Path, out_path: Path, max_rows: int | None = None, delimiter=None):
    stats = DatasetStats(); families = {}; peptides=set(); alleles=set()
    # Validate headers BEFORE creating an output file. Use bounded sample only.
    with path.open(newline='',encoding='utf-8-sig',errors='strict') as f:
        first=f.readline()
        sep=delimiter or ('\t' if first.count('\t')>first.count(',') else ',')
        f.seek(0);reader=csv.reader(f,delimiter=sep)
        h1,h2=next(reader),next(reader)
        if len(h1)!=len(h2) or len(h1)<8:
            raise ValueError('IEDB requires matching two-row CSV/TSV headers')
        cols={(a.strip().casefold(),b.strip().casefold()):i for i,(a,b) in enumerate(zip(h1,h2))}
        def ix(a,b,required=True):
            v=cols.get((a.casefold(),b.casefold()))
            if required and v is None:raise ValueError(f'missing IEDB column {a}|{b}')
            return v
        im=ix('Assay','Method');ip=ix('Epitope','Name');ia=ix('MHC Restriction','Name')
        iq=ix('Assay','Qualitative Measurement');iv=ix('Assay','Quantitative measurement')
        iu=ix('Assay','Units',False);ih=ix('Host','Name');ir=ix('Reference','PMID')
        iid=ix('Assay ID','IEDB IRI',False)
        iqual=ix('Assay','Quantitative measurement inequality',False)
        out_path.parent.mkdir(parents=True,exist_ok=True)
        tmp=out_path.with_suffix(out_path.suffix+'.tmp')
        rejected=out_path.with_suffix(out_path.suffix+'.rejected.tsv')
        rejected_tmp=rejected.with_suffix(rejected.suffix+'.tmp')
        with tmp.open('w',newline='') as target, rejected_tmp.open('w',newline='') as reject:
            w=csv.writer(target,delimiter='\t');w.writerow(SCHEMA_HEADER)
            rw=csv.writer(reject,delimiter='\t');rw.writerow(['row','assay_id','reason'])
            for n,line in enumerate(reader,1):
                if max_rows is not None and n>max_rows:break
                if len(line)!=len(h1):
                    rw.writerow([n,'','column_count']);continue
                method=line[im].strip();host=line[ih].strip().lower()
                if host and 'homo sapiens' not in host and 'human' not in host:
                    rw.writerow([n,line[iid] if iid is not None else '', 'nonhuman_host']);continue
                if not method or 't cell' in method.lower() or not any(k in method.lower() for k in ['mhc','binding','ic50','mass spectr','elut','koff','kon','stability','dissociation']):
                    rw.writerow([n,line[iid] if iid is not None else '', 'out_of_mhc_assay_scope']);continue
                peptide=line[ip].strip();allele=line[ia].strip()
                if not peptide or not allele:
                    rw.writerow([n,'','missing_peptide_or_allele']);continue
                # Ambiguous allele sets are retained as sets, never first-allele substitution.
                family=assay_family(method);qual=line[iq].strip().lower()
                label=1 if qual.startswith('positive') else 0 if qual.startswith('negative') else None
                raw=line[iv].strip();unit=line[iu].strip() if iu is not None else ''
                qualifier=line[iqual].strip() if iqual is not None else ''
                score,note=(None,'raw quantitative evidence retained')
                if family in {'binding-ic50','binding-kd'} and not qualifier:
                    score,note=ic50_score(raw,unit)
                rec=InteractionRecord(source=DataSource.IEDB,peptide=peptide,mhc_allele=allele,label=label,score=score,
                    reference_id=line[ir].strip(),organism='HomoSapiens',method='iedb:'+family,notes=note,
                    assay_id=line[iid].strip() if iid is not None else f'input-row:{n}',endpoint=family,
                    evidence_status='assay_positive' if label==1 else 'assay_negative' if label==0 else 'assay_unknown',
                    label_scope='mhc_presentation' if family=='presentation-eluted' else 'mhc_assay',
                    raw_value=raw,raw_unit=unit,qualifier=qualifier)
                w.writerow(rec.to_row());stats.n_records+=1;stats.n_labeled+=int(label is not None)
                peptides.add(peptide);alleles.add(allele);families[family]=families.get(family,0)+1
                if n%500000==0:print(f'IEDB source_rows={n} retained={stats.n_records}',flush=True)
        tmp.replace(out_path);rejected_tmp.replace(rejected)
    stats.n_unique_peptides=len(peptides);stats.n_alleles=len(alleles)
    stats.n_unlabeled=stats.n_records-stats.n_labeled
    stats.warnings += [f'assay families: {families}', 'assay-level ledger: conflicting results retained; MHC evidence is not TCR binding',f'delimiter={sep!r}']
    return stats


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--in',dest='in_path',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--max-rows',type=int)
    p.add_argument('--delimiter',choices=['csv','tsv']);a=p.parse_args()
    print(ingest(a.in_path,a.out,a.max_rows,{'csv':',','tsv':'\t'}.get(a.delimiter)).summary())
    return 0
if __name__=='__main__':raise SystemExit(main())
