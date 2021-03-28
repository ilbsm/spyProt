from spyprot import MMCIFfile
from spyprot.cif2Wanda import run_cif2Wanda
import tempfile
import os

dir = os.path.join(tempfile.gettempdir(), 'test_pdb2wanda')
os.makedirs(dir, exist_ok=True)
pdbid = '6xj0'
p = MMCIFfile(dir, pdbid).download()

if os.path.exists(p):
    (bridges, gaps_dict, str_begin_dict) = run_cif2Wanda(p,
                                                         dir,
                                                         pdbid)  # return tuple of list of tuples: (bridge type, chain_name, bridge_begin, bridge_end) and dict with gaps
    print(bridges)
    print(gaps_dict)
    print(str_begin_dict)
