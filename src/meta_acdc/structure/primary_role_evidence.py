"""Literal depositor role evidence; aliases do not repair missing sequence or atoms."""
from pathlib import Path
import re

def classify_role(molecule,synonym=''):
 raw=molecule+' '+synonym;text=raw.upper().replace('Α','ALPHA').replace('Β','BETA');text=re.sub(r'\bALFA\b','ALPHA',text)
 if re.search(r'MICROGLOBULIN|\bB2M\b',text):return 'b2m'
 tcr=bool(re.search(r'\bTCR\b|T[- ]CELL|T[- ]LYMPHOCYTE|ANTIGEN RECEPTOR',text))
 if not tcr and re.search(r'\bHLA\b|\bMHC\b|CLASS I',text):return 'mhc'
 a=bool(re.search(r'\b(?:V[- ]?)?ALPHA[0-9.]*\b|\bTRAC\b',text)) or (tcr and 'LIGHT CHAIN' in text)
 b=bool(re.search(r'\b(?:V[- ]?)?BETA[0-9.]*\b|\bTRBC[12]?\b',text)) or (tcr and 'HEAVY CHAIN' in text)
 if a and b:return 'conflicting_alpha_beta'
 if a:return 'alpha'
 if b:return 'beta'
 return 'unspecified'

def primary_roles(path:Path):
 parts=[]
 for line in path.read_text().splitlines():
  if line.startswith('COMPND'):parts.append(line[10:].strip())
  if line.startswith(('ATOM  ','HETATM')):break
 entities=[];current={}
 for token in ' '.join(parts).split(';'):
  if ':' not in token:continue
  key,value=token.split(':',1);key=key.strip();value=value.strip()
  if key=='MOL_ID':
   if current:entities.append(current)
   current=dict(mol_id=value)
  else:current[key.lower()]=value
 if current:entities.append(current)
 out={}
 for entity in entities:
  molecule=entity.get('molecule','');synonym=entity.get('synonym','');kind=classify_role(molecule,synonym)
  for chain in entity.get('chain','').split(','):
   chain=chain.strip()
   if not chain:continue
   if chain in out:raise ValueError('duplicate primary author chain')
   out[chain]=dict(role=kind,molecule=molecule,synonym=synonym,mol_id=entity['mol_id'],basis='literal primary COMPND molecule plus synonym; ALFA spelling normalized; BETA-2 alone is not microglobulin')
 return out
