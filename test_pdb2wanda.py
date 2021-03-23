from spyprot import PdbFile
from spyprot.pdb2Wanda import run_pdb2Wanda
import tempfile
import os

dir = os.path.join(tempfile.gettempdir(), 'test_pdb2wanda')
os.makedirs(dir, exist_ok=True)
pdbid = '1j85'
p = PdbFile(dir, pdbid).download()

if os.path.exists(p):
    with open(p, 'r') as file:
        content = file.read()
    data = content.split('\n')

    (bridges, gaps_dict, str_begin_dict) = run_pdb2Wanda(data,
                                                         dir,
                                                         pdbid)  # return tuple of list of tuples: (bridge type, chain_name, bridge_begin, bridge_end) and dict with gaps
    print(bridges)
    print(gaps_dict)
    print(str_begin_dict)