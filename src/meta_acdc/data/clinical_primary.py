"""Clinical primary-source curation: never turn cohort outcomes into receptor labels."""
from __future__ import annotations
from itertools import product
from pathlib import Path
import hashlib,io,json,re,zipfile
from xml.etree import ElementTree as ET

CODONS=dict(zip(map("".join,product("TCAG",repeat=3)),
 "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
AA=set("ACDEFGHIKLMNPQRSTVWY")

def translate_fragment(dna,cdr3):
    """Publish frame-zero translation only; retain partial terminal codon explicitly."""
    if not isinstance(dna,str) or not re.fullmatch("[ACGT]+",dna):raise ValueError("noncanonical DNA")
    aa="".join(CODONS[dna[i:i+3]] for i in range(0,len(dna)-2,3))
    if "*" in aa:raise ValueError("internal stop in variable fragment")
    if not cdr3 or not set(cdr3)<=AA:raise ValueError("noncanonical CDR3")
    return dict(dna=dna,amino_acids=aa,frame=0,trailing_nt=dna[len(aa)*3:],cdr3_matches=cdr3 in aa,
                scope="published_variable_fragment; not complete engineered chain")

def read_workbook(path):
    ns={"s":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as z:
        strings=[]
        if "xl/sharedStrings.xml" in z.namelist():
            strings=["".join(e.itertext()) for e in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("s:si",ns)]
        names=[x for x in z.namelist() if re.fullmatch(r"xl/worksheets/sheet[0-9]+.xml",x)]
        if len(names)!=1:raise ValueError("expected one explicitly reviewed sheet")
        rows=[]
        for row in ET.fromstring(z.read(names[0])).findall(".//s:sheetData/s:row",ns):
            values={}
            for c in row.findall("s:c",ns):
                if c.find("s:f",ns) is not None:raise ValueError("formula cannot be silently evaluated")
                v=c.find("s:v",ns);v=v.text if v is not None else None
                if c.get("t")=="s":v=strings[int(v)]
                elif c.get("t")=="inlineStr":v="".join(c.find("s:is",ns).itertext())
                values[re.match("[A-Z]+",c.get("r"))[0]]=v
            rows.append((int(row.get("r")),values))
        return rows

def foy_receptors(path):
    rows=read_workbook(path)
    if rows[0][1].get("J")!="TCR alpha variable seq " or rows[0][1].get("K")!="TCR beta variable seq ":
        raise ValueError("unreviewed sequence source headers")
    out=[]
    for line,x in rows[1:]:
        a=translate_fragment(x["J"],x["H"]);b=translate_fragment(x["K"],x["I"])
        ok=a["cdr3_matches"] and b["cdr3_matches"]
        out.append(dict(receptor_name="Foy:"+x["G"],lineage="Foy:"+x["G"],patient_id=x["B"],
            cohort_id="NCT03970382:36356599",pmid="36356599",excel_row=line,
            mhc_allele=x["C"].replace("HLA-A","HLA-A*").replace("HLA-B","HLA-B*").replace("HLA-C","HLA-C*"),
            peptide=x["F"],target_gene=x["D"],mutation=x["E"],
            cdr3_alpha=x["H"],cdr3_beta=x["I"],alpha_variable=a["amino_acids"],beta_variable=b["amino_acids"],
            alpha_fragment=a,beta_fragment=b,identity_status="variable_pair_verified" if ok else "quarantined_CDR3_DNA_conflict",
            complete_engineered_sequence=False,clinical_training_allowed=False,clinical_risk_label=None,
            role="retrospective_reserved_candidate",source_locator="Supplementary Table 2!A%d:K%d"%(line,line)))
    if len(out)!=37 or len({x["receptor_name"] for x in out})!=37:raise ValueError("unreviewed receptor count/duplicates")
    return out

def e7_construct(text):
    lines=text.split("Supplemental Figure 1.")[0].splitlines()
    segments=[];current=[]
    for line in lines:
        line=line.strip()
        if re.fullmatch("[ACGT]+",line):current.append(line)
        elif current:segments.append("".join(current));current=[]
    if current:segments.append("".join(current))
    if len(segments)!=5:raise ValueError("unreviewed E7 sequence layout")
    dna="".join(segments)
    if len(dna)%3 or not dna.startswith("ATG"):raise ValueError("invalid full E7 ORF")
    aa="".join(CODONS[dna[i:i+3]] for i in range(0,len(dna),3))
    if not aa.endswith("*") or "*" in aa.rstrip("*"):raise ValueError("invalid E7 stop")
    for cdr in ["CASVDGNNRLAF","CASSLGWRGGRYNEQF"]:
        if cdr not in aa:raise ValueError("E7 CDR3 mismatch")
    return dict(receptor_name="HPV16-E7-optimized",lineage="HPV16-E7",dna=dna,protein=aa.rstrip("*"),terminal_stop_codons=len(aa)-len(aa.rstrip("*")),segments=segments,
        sequence_scope="complete published preclinical beta-furinP2A-alpha coding insert; vector backbone not included",
        alpha_variable="".join(CODONS[segments[3][i:i+3]] for i in range(0,len(segments[3]),3)),
        beta_variable="".join(CODONS[segments[0][i:i+3]] for i in range(0,len(segments[0]),3)),
        cdr3_alpha="CASVDGNNRLAF",cdr3_beta="CASSLGWRGGRYNEQF",complete_published_coding_insert=True,
        clinical_lot_sequence_verified=False,clinical_product_link="reference-linked optimized E7 receptor; no public lot-level sequencing certificate",
        primary_pmid="29669936",clinical_pmid="33558725",trial="NCT02858310",
        engineering=["mouse constants","alpha T48C","beta S57C","alpha TM LSVMGLRIL -> LLVIVLRIL","beta-alpha order","furin-P2A"],
        source_typo="Fifth block headed Beta chain in source; follows alpha variable and corresponds to alpha constant per Methods. Raw text retained.",
        clinical_training_allowed=False,clinical_risk_label=None,role="retrospective_reserved_candidate")

def validate_bundle(bundle):
    rec=bundle["receptors"];patients=bundle["patients"]
    if len({p["patient_id"] for p in patients})!=len(patients):raise ValueError("duplicate patient")
    by={p["patient_id"]:p for p in patients}
    for x in rec:
        if x.get("clinical_training_allowed") is not False or x.get("clinical_risk_label") is not None:
            raise ValueError("clinical outcomes cannot be receptor risk labels")
        if x.get("cohort_id")=="NCT03970382:36356599":
            if x["patient_id"] not in by:raise ValueError("unlinked patient")
            p=by[x["patient_id"]]
            if x["receptor_name"] not in p["receptors"]:raise ValueError("patient/product join mismatch")
            for kind in ["alpha","beta"]:
                fragment=x[kind+"_fragment"]
                check=translate_fragment(fragment["dna"],x["cdr3_"+kind])
                if check!=fragment or x[kind+"_variable"]!=check["amino_acids"]:raise ValueError("translated sequence changed")
            match=all(x[k]["cdr3_matches"] for k in ["alpha_fragment","beta_fragment"])
            if (x["identity_status"]=="variable_pair_verified") != match:raise ValueError("identity quarantine bypass")
    for p in patients:
        expected=sorted(x["receptor_name"] for x in rec if x.get("patient_id")==p["patient_id"])
        if sorted(p["receptors"])!=expected or p["n_receptors"]!=len(expected):raise ValueError("combination attribution mismatch")
        if p.get("per_receptor_causality_established") is not False:raise ValueError("unsupported causal attribution")
    return True
