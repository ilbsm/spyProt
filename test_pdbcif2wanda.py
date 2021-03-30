from spyprot import PdbFile
from spyprot.pdb2Wanda import run_pdb2Wanda
from spyprot import MMCIFfile
from spyprot.cif2Wanda import run_cif2Wanda
import tempfile
import os


def test_pdb2wanda():
    dir = os.path.join(tempfile.gettempdir(), 'test_pdb2wanda')
    try:
        os.makedirs(dir)
    except:
        pass
    pdbid = '6xj0'
    p = PdbFile(dir, pdbid).download()

    if os.path.exists(p):
        (bridges, gaps_dict, str_begin_dict) = run_pdb2Wanda(p,
                                                             dir,
                                                             pdbid)  # return tuple of list of tuples: (bridge type, chain_name, bridge_begin, bridge_end) and dict with gaps
        assert bridges == [('SS', 'A', 121, 151)]
        assert gaps_dict == {'A': ''}
        assert str_begin_dict == {'A': 2}


def test_cif2wanda():
    dir = os.path.join(tempfile.gettempdir(), 'test_pdb2wanda')
    try:
        os.makedirs(dir)
    except:
        pass
    pdbid = '6xj0'
    p = MMCIFfile(dir, pdbid).download()

    if os.path.exists(p):
        (bridges, gaps_dict, str_begin_dict) = run_cif2Wanda(p,
                                                             dir,
                                                             pdbid)  # return tuple of list of tuples: (bridge type, chain_name, bridge_begin, bridge_end) and dict with gaps
        assert bridges == [('SS', 'A', 121, 151)]
        assert gaps_dict == {'A': ''}
        assert str_begin_dict == {'A': 2}
