"""Relocate the unmodified Ubuntu SPAdes package; change paths only, never assembly logic."""
from pathlib import Path
import os,sys,runpy
root=Path(os.environ['NP105_SPADES_ROOT']).resolve();binary=root/'usr/libexec/spades';shared=root/'usr/share/spades';sys.path.insert(0,str(binary));import spades_init
def relocated_init():
    spades_init.spades_home=str(shared);spades_init.bin_home=str(binary);spades_init.python_modules_home=str(shared);spades_init.ext_python_modules_home=str(shared)
    for d in spades_init.source_dirs:sys.path.append(str(shared/'spades_pipeline'/d))
    spades_init.spades_version=(shared/'VERSION').read_text().strip()
spades_init.init=relocated_init
runpy.run_path(str(binary/'spades.py'),run_name='__main__')
