"""Sequence-linked primary-paper adjudication of explicit metadata conflicts."""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import xml.etree.ElementTree as ET

from meta_acdc.structure.paired_variable_identity import germline_references, validate_domain, pair_key
from meta_acdc.structure.primary_role_evidence import primary_roles
from meta_acdc.structure.receptor_identity_audit import seqres_sequences

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'docs/fingerprint_construct_followup_20261027'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def linked_roles(sequences, alpha_motifs, beta_motifs):
    """All published motifs must jointly select exactly one distinct author chain."""
    selected = {}
    for role, motifs in [('alpha', alpha_motifs), ('beta', beta_motifs)]:
        if not motifs or any(not m for m in motifs):
            raise ValueError('missing published sequence anchors')
        candidates = [chain for chain, seq in sequences.items() if all(m in seq for m in motifs)]
        if len(candidates) != 1:
            raise ValueError('published sequence anchors do not uniquely select '+role)
        selected[role] = candidates[0]
    if selected['alpha'] == selected['beta']:
        raise ValueError('same author chain cannot establish paired roles')
    return selected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--include-natural', action='store_true')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    articles = {stem: ET.parse(BASE / 'sources/secondary' / (stem+'_article.xml')).getroot()
                for stem in ['6at6', '7n5']}
    for stem, pmcid in [('6at6', 'PMC5847591'), ('7n5', 'PMC9399087')]:
        assert any(a.text == pmcid for a in articles[stem].iter('article-id'))
    table = next(t for t in articles['6at6'].iter('table-wrap') if t.get('id') == 'Tab1')
    kfj5 = next([''.join(cell.itertext()).strip() for cell in row]
                for row in table.findall('.//tr') if ''.join(row[0].itertext()).strip() == 'KFJ5')
    alpha_kfj, beta_kfj = kfj5[3:6], kfj5[10:13]
    assert kfj5[1] == '4*01' and kfj5[6] == '28*01'
    methods = next(sec for sec in articles['7n5'].iter('sec')
                   if sec.findtext('title') == 'TCR and peptide/MHC synthesis')
    method_text = ''.join(methods.itertext())
    assert 'Ser110α' in method_text and 'Cys110α' in method_text
    availability = ' '.join(articles['7n5'].itertext())
    for code in ['7N5C', '7N5P', '7N4K']:
        assert code in availability
    # Transcribed from the rendered primary Figure 1a; retained with source image hash.
    motifs = {'6at6': (alpha_kfj, beta_kfj),
              '7n5c': (['ILSGGCNYKLT'], ['ASSFGREQY']),
              '7n5p': (['ILSGGSNYKLT'], ['ASSFGREQY'])}
    natural = None
    if args.include_natural:
        natural = json.loads((BASE/'natural_control_primary_receipt.json').read_text())
        source = BASE/'sources/natural/AS_supplement.pdf'
        if source.exists():
            assert sha(source) == natural['supplement_sha256']
        clone = next(c for c in natural['clones'] if c['name'] == 'AS8.4')
        assert clone['TRAV'] == 'TRAV21' and clone['TRBV'] == 'TRBV9'
        motifs['8cx4'] = ([clone['CDR3_alpha']], [clone['CDR3_beta']])
    image = BASE / 'sources/retry/7n5_figure1.jpg'
    assert image.read_bytes().startswith((b'\xff\xd8', b'\x89PNG\r\n\x1a\n')), 'primary figure is not an image'
    old = json.loads((ROOT / 'docs/fingerprint_cross_source_extension_20261024/sequence_curation_merged.json').read_text())
    references = ROOT / 'docs/fingerprint_paired_identity_20261024/sources'
    refs = germline_references(list(references.glob('Homo_sapiens_TR*.fasta')) + list(references.glob('Mus_musculus_TR*.fasta')))
    annotation_dirs = [ROOT / 'docs/fingerprint_cross_source_extension_20261024/sources/stcrdab_all',
                       ROOT / 'docs/fingerprint_cross_source_extension_20261024/sources/stcrdab',
                       ROOT / 'docs/fingerprint_cross_source_20261024/sources/stcrdab']
    rows, receipts = [], []
    for pdb in motifs:
        coordinate = BASE / ('sources/secondary/7n5p.pdb' if pdb == '7n5p' else f'sources/tertiary/{pdb}.pdb')
        if pdb == '8cx4':
            coordinate = BASE/'sources/natural_coordinate/8cx4.pdb'
        assert coordinate.read_text().startswith('HEADER')
        sequences = seqres_sequences(coordinate.read_text())
        roles = linked_roles(sequences, *motifs[pdb])
        original = next(r for r in old['rows'] if r['pdb'] == pdb and
                        r['alpha_chain'] == roles['alpha'] and r['beta_chain'] == roles['beta'])
        assert original['status'] == 'sequence_supported_pending_primary_header' and not original['reasons']
        row = copy.deepcopy(original)
        declarations = primary_roles(coordinate)
        for role in ['alpha', 'beta']:
            chain = roles[role]
            assert sequences[chain] == row[role]['full_construct']
            annotations = [d for d in annotation_dirs if (d/(pdb+'.html')).exists() and
                           (d/(pdb+'_'+row[role]['annotation_chain']+'.fasta')).exists()]
            fasta = next(d/(pdb+'_'+row[role]['annotation_chain']+'.fasta') for d in annotations
                         if sha(d/(pdb+'_'+row[role]['annotation_chain']+'.fasta')) == row[role]['annotation_sha256'])
            html = fasta.parent / (pdb+'.html')
            validated = validate_domain(row[role]['full_construct'], row[role]['domain_candidate'],
                                        fasta.read_text(), html.read_text(), refs)
            assert validated['variable'] == row[role]['variable']
            row[role]['primary_role_evidence'] = dict(
                author_chain=chain, deposited_declaration=declarations[chain],
                adjudication='independent primary-paper sequence anchors uniquely identify role; literal conflict retained',
                paper_motifs=motifs[pdb][0 if role == 'alpha' else 1],
                full_construct_sha256=hashlib.sha256(sequences[chain].encode()).hexdigest())
        row.update(status='identity_admitted_structure_QC_pending', primary_header_sha256=sha(coordinate),
                   primary_header=str(coordinate.relative_to(ROOT)),
                   adjudication_scope=('additional natural-clone sequence receipt' if pdb == '8cx4' else 'three explicit sequence-linked paper receipts')+'; general lexical parser unchanged')
        assert pair_key(row['alpha']['variable'], row['beta']['variable']) == row['receptor_id']
        rows.append(row)
        stem = '6at6' if pdb == '6at6' else '7n5'
        if pdb == '8cx4':
            assert row['V_context'] == 'human|TRAV21|TRBV9'
            receipts.append(dict(pdb=pdb,alpha_chain=roles['alpha'],beta_chain=roles['beta'],
                construct_name='AS8.4, patient K',primary_doi=natural['doi'],
                supplement_sha256=natural['supplement_sha256'],locator=natural['locator'],
                factual_transcript='natural_control_primary_receipt.json',
                paired_roles_uniquely_sequence_linked=True,full_variable_boundary_policy_unchanged=True,
                raw_conflict_not_erased=True))
            continue
        receipts.append(dict(pdb=pdb, alpha_chain=roles['alpha'], beta_chain=roles['beta'],
                             construct_name={'6at6': 'KFJ5, unbound', '7n5c': '6218alphaC', '7n5p': '6218'}[pdb],
                             primary_doi=row['source_study_key'], primary_xml_sha256=sha(BASE / f'sources/secondary/{stem}_article.xml'),
                             locator='Table 1 and protein-expression methods' if pdb == '6at6' else 'Figure 1a; TCR and peptide/MHC synthesis; data availability; SI Table 1 p12',
                             primary_figure_sha256=sha(image) if pdb != '6at6' else None,
                             paired_roles_uniquely_sequence_linked=True,
                             full_variable_boundary_policy_unchanged=True, raw_conflict_not_erased=True))
    mutant = next(r for r in rows if r['pdb'] == '7n5c')
    parent = next(r for r in rows if r['pdb'] == '7n5p')
    a, b = mutant['alpha']['variable'], parent['alpha']['variable']
    assert len(a) == len(b)
    differences = [dict(variable_index0=i, parent=y, engineered=x)
                   for i, (x, y) in enumerate(zip(a, b)) if x != y]
    assert len(differences) == 1 and differences[0]['parent'] == 'S' and differences[0]['engineered'] == 'C'
    assert mutant['beta']['variable'] == parent['beta']['variable']
    result = dict(rows=rows, primary_receipts=receipts, alpha_mutation=differences,
                  natural_independent_control=False,
                  reason='7N5C is S110C-engineered 6218 from the same lineage as 7N5P; 6AT6 is unbound KFJ5 and shares its complete variable identity with same-study 6AVF',
                  protocol_frozen=False, model_training_started=False, independent_acceptance_accessed=False)
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(n_roles_adjudicated=len(rows), mutation=differences, boundaries_changed=False), indent=2))


if __name__ == '__main__':
    main()
