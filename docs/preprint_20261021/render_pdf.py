"""Render the editable Markdown author-review drafts using ReportLab."""
from pathlib import Path
import argparse,hashlib,html,json,re,sys,textwrap
from importlib.metadata import version
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER,TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image,Table,TableStyle,KeepTogether
HERE=Path(__file__).resolve().parent
STYLES=getSampleStyleSheet()
STYLES.add(ParagraphStyle(name='PaperTitle',fontName='Times-Bold',fontSize=18,leading=22,spaceAfter=16))
STYLES.add(ParagraphStyle(name='PaperBody',fontName='Times-Roman',fontSize=10.5,leading=14.3,spaceAfter=8))
STYLES.add(ParagraphStyle(name='PaperH2',fontName='Times-Bold',fontSize=13,leading=17,spaceBefore=13,spaceAfter=7,keepWithNext=True))
STYLES.add(ParagraphStyle(name='PaperH3',fontName='Times-Bold',fontSize=11,leading=14,spaceBefore=10,spaceAfter=6,keepWithNext=True))
STYLES.add(ParagraphStyle(name='PaperCaption',fontName='Times-Roman',fontSize=9,leading=12,spaceAfter=12))
STYLES.add(ParagraphStyle(name='TableCell',fontName='Times-Roman',fontSize=8,leading=10))
STYLES.add(ParagraphStyle(name='CodeLine',fontName='Courier',fontSize=7.5,leading=10,spaceAfter=3))
WIDTH=A4[0]-42*mm

def markup(s):
    s=html.escape(s)
    s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s)
    s=re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',r'<link href="\2" color="#26547c">\1</link>',s)
    return s

def footer(c,doc):
    c.saveState();c.setStrokeColor(colors.HexColor('#b8c1ca'));c.line(21*mm,19*mm,A4[0]-21*mm,19*mm)
    c.setFont('Helvetica',8);c.setFillColor(colors.HexColor('#52616c'))
    c.drawString(21*mm,A4[1]-14*mm,'AUTHOR REVIEW DRAFT | Meta-ACDC case study')
    c.drawString(21*mm,14*mm,'21 October 2026 | Not submitted | Author details pending')
    c.drawRightString(A4[0]-21*mm,14*mm,str(doc.page));c.restoreState()

def parse(path):
    lines=path.read_text(encoding='utf-8').splitlines();story=[];i=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:i+=1;continue
        if line.startswith('# '):story.append(Paragraph(markup(line[2:]),STYLES['PaperTitle']));i+=1;continue
        if line.startswith('## '):story.append(Paragraph(markup(line[3:]),STYLES['PaperH2']));i+=1;continue
        if line.startswith('### '):story.append(Paragraph(markup(line[4:]),STYLES['PaperH3']));i+=1;continue
        image=re.fullmatch(r'!\[(.+)\]\(([^)]+)\)',line)
        if image:
            from PIL import Image as PILImage
            p=HERE/image[2];w,h=PILImage.open(p).size;im=Image(str(p),width=WIDTH,height=WIDTH*h/w)
            story.append(KeepTogether([im,Spacer(1,4),Paragraph(markup(image[1]),STYLES['PaperCaption'])]));i+=1;continue
        if line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                cells=[x.strip() for x in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r':?-+:?',x) for x in cells):rows.append(cells)
                i+=1
            n=len(rows[0]);widths=[WIDTH/n]*n
            if n==7:widths=[WIDTH*x for x in [.19,.16,.12,.10,.10,.165,.165]]
            elif n==4:widths=[WIDTH*x for x in [.13,.21,.10,.56]]
            body=[[Paragraph(markup(x),STYLES['TableCell']) for x in row] for row in rows]
            t=Table(body,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e7eef3')),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#657887')),('LINEBELOW',(0,-1),(-1,-1),.4,colors.HexColor('#9ba9b3')),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('LEFTPADDING',(0,0),(-1,-1),4),('RIGHTPADDING',(0,0),(-1,-1),4)]));story.extend([KeepTogether([t]),Spacer(1,10)]);continue
        if lines[i].startswith('    '):
            while i<len(lines) and lines[i].startswith('    '):
                fragments=[Paragraph(markup(fragment),STYLES['CodeLine']) for fragment in textwrap.wrap(lines[i][4:],width=96,break_long_words=True,break_on_hyphens=False)]
                story.append(KeepTogether(fragments))
                i+=1
            story.append(Spacer(1,7));continue
        parts=[]
        while i<len(lines) and lines[i].strip() and not lines[i].startswith(('#','|','![','    ')):
            parts.append(lines[i].strip());i+=1
        if parts:
            text=' '.join(parts)
            para=Paragraph(markup(text),STYLES['PaperBody'])
            if text.startswith('Table S'):para.keepWithNext=True
            story.append(KeepTogether([para]) if re.match(r'^\d+\. ',text) else para)
    return story

def render(out):
    out=Path(out)
    if out.exists() and any(out.iterdir()):raise ValueError('choose a fresh PDF output directory')
    out.mkdir(parents=True,exist_ok=True);outputs={}
    for name in ['manuscript','supplement']:
        p=out/(name+'.pdf');doc=SimpleDocTemplate(str(p),pagesize=A4,leftMargin=21*mm,rightMargin=21*mm,topMargin=22*mm,bottomMargin=24*mm,title='Meta-ACDC: '+name+' author-review draft',author='Single author; identifying information pending',allowSplitting=True)
        doc.build(parse(HERE/(name+'.md')),onFirstPage=footer,onLaterPages=footer)
        outputs[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    receipt=dict(python=sys.version,reportlab=version('reportlab'),inputs={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['manuscript.md','supplement.md','render_pdf.py']},outputs=outputs,status='author-review draft; not submitted')
    (out/'render_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args();render(a.out)
