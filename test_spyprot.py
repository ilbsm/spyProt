import warnings

from spyprot import getCath
from spyprot.fetchPDBinfo import fetchPDBinfo
from spyprot.fetchAnnotations import PfamAnnotation, ECAnnotation


# def test_get_annotations():
#     a = getAnnotations("1j85", "A")
#     assert str(a.get()[0]) =="['1j85', 'A', '2', '142', 'PF00588.18', 'SpoU_methylase', 'SpoU rRNA Methylase family', '4.3E-32']"


def test_get_ec_annotations():
    a = ECAnnotation("1cak", "A")
    assert str(a.get()) == "[('4.2.1.1', 'Carbonic anhydrase.')]"


def test_get_pfam_annotations():
    a = PfamAnnotation("1j85", "A")
    assert str(a.get()) == "[('1j85', 'A', 'SpoU_methylase', 'SpoU rRNA Methylase family', 'PF00588')]"


def test_fetch_pdb_info():
    from Bio import BiopythonExperimentalWarning
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonExperimentalWarning)
        a = fetchPDBinfo("1uak", "A")
        assert str(a.getCalfaBreaks()) == '[162, 172]'


def test_cddf_parser():
    d = getCath("1j85", "A")
    assert d[0]['DOMAIN'] == '1j85A00' and d[0]['CATHCODE'] == '3.40.1280.10' and d[0]['CLASS'] == 'Alpha Beta'
