"""Check original observed J/terminal intervals in whole-run SAM; never fill missing bases."""
from pathlib import Path
import argparse, collections, gzip, hashlib, json, sys
sys.path[:0] = [".", "src"]
from scripts.validate_np105_raw import aligned_bases
from scripts.audit_np105_whole_variants import haplotype_umis, CDR3
ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "docs/fingerprint_np105_recovery_20261029"
PREFIX = ROOT / "docs/fingerprint_np105_raw_20261029"
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--sam",type=Path,required=True);ap.add_argument("--out",type=Path,required=True);a=ap.parse_args()
    classic=json.loads(gzip.decompress((OLD/"evidence/whole_classic_numbering.json.gz").read_bytes()))
    prior=json.loads((PREFIX/"raw_sequence_validation.json").read_text())["domains"]
    prefix_inputs=json.loads((PREFIX/"raw_classic_numbering.json").read_text())["inputs"]
    bases=json.loads(gzip.decompress((OLD/"evidence/whole_base_support.json.gz").read_bytes()))
    paired={x["barcode"]:x for x in bases["cells"]};cells={x["barcode"]:x for x in classic["cells"]}
    targets={}
    for role,start,end in [("TRA",274,336),("TRB",296,349)]:
        ref=prior[role];observed=next(x["ORF_nt"] for x in prefix_inputs if x["role"]==role)
        sequence=observed[ref["ORF_start_nt0"]+start:ref["ORF_start_nt0"]+end]
        assert len(sequence)==end-start
        targets[role]=dict(start=start,end=end,observed_sequence=sequence,offset=ref["ORF_start_nt0"],cdr3_start=ref["ORF_start_nt0"]+ref["nucleotides"].index(CDR3[role]))
    calls=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(set)))
    linked=collections.defaultdict(set)
    with a.sam.open() as stream:
        for line in stream:
            if line.startswith("@"):continue
            f=line.rstrip().split("\t");flag=int(f[1])
            if flag & (4|256|2048) or int(f[4])<30 or f[2] not in targets:continue
            bc=f[0].split(".BC_")[1].split(".UMI_")[0]
            if bc not in paired:continue
            umi=f[0].split(".UMI_")[1];role=f[2];target=targets[role];key=(bc,role)
            aligned={i:(base,q) for i,base,q in aligned_bases(int(f[3]),f[5],f[9],f[10])}
            for pos in range(target["start"],target["end"]):
                value=aligned.get(target["offset"]+pos)
                if value and value[1]>=30:calls[key][umi][pos].add(value[0])
            if all(aligned.get(target["cdr3_start"]+i,(None,0))[0]==base and aligned.get(target["cdr3_start"]+i,(None,0))[1]>=30 for i,base in enumerate(CDR3[role])):linked[key].add(umi)
    rows=[]
    for bc,quality in sorted(paired.items()):
        row=dict(barcode=bc,original_canonical_Q30_pair=quality["observed_canonical_pair_supported"],roles={})
        for role,target in targets.items():
            positions=list(range(target["start"],target["end"]));umis=haplotype_umis(calls[(bc,role)],positions,target["observed_sequence"])
            full=[x for x in cells[bc][role] if x["complete"]]
            observations=[]
            for x in full:
                start=x["boundary"]["start0"]*3
                actual=x["ORF_nt"][start+target["start"]:start+target["end"]]
                observations.append(dict(name=x["name"],length=len(actual),exact_observed_target=actual==target["observed_sequence"]))
            assembled=bool(observations) and all(x["exact_observed_target"] for x in observations)
            row["roles"][role]=dict(independently_assembled_all_complete_ORFs_have_exact_target=assembled,assembly_observations=observations,unambiguous_Q30_UMIs_cover_all_positions=len(umis),also_exact_CDR3_linked_UMIs=len(umis & linked[(bc,role)]),terminal_supported=assembled and bool(umis & linked[(bc,role)]))
        row["both_terminal_targets_supported"]=all(x["terminal_supported"] for x in row["roles"].values());rows.append(row)
    summary=dict(paired_cells=60,original_canonical_Q30_pairs=54,both_terminal_targets_supported=sum(x["both_terminal_targets_supported"] for x in rows),original_canonical_Q30_pairs_both_terminals_supported=sum(x["original_canonical_Q30_pair"] and x["both_terminal_targets_supported"] for x in rows),by_role={role:sum(x["roles"][role]["terminal_supported"] for x in rows) for role in targets})
    report=dict(summary=summary,targets=targets,cells=rows,SAM_bytes=a.sam.stat().st_size,SAM_SHA256=hashlib.file_digest(a.sam.open("rb"),"sha256").hexdigest(),source_whole_MD5_verified=True,reference_bases_added=0,new_identity_admissions=0,protocol_frozen=False,independent_acceptance_accessed=False,interpretation="Observed terminal support only; UMI linkage is not genomic phase or proof against barcode/UMI collisions; original54/26/55 gates unchanged")
    a.out.write_text(json.dumps(report,indent=2)+"\n");print(json.dumps(summary,indent=2))
if __name__=="__main__":main()
