import pytest

from spyprot import getCath
from spyprot.fetchPDBinfo import fetchPDBinfo
from spyprot.fetchAnnotations import PfamAnnotation, ECAnnotation


# def test_get_annotations():
#     a = getAnnotations("1j85", "A")
#     assert str(a.get()[0]) =="['1j85', 'A', '2', '142', 'PF00588.18', 'SpoU_methylase', 'SpoU rRNA Methylase family', '4.3E-32']"


@pytest.mark.skip(reason="EC server is not reliable")
def test_get_ec_annotations():
    a = ECAnnotation("1cak", "A")
    assert str(a.get()) == "[('4.2.1.1', 'carbonic anhydrase.'), ('4.2.1.69', 'cyanamide hydratase.')]"


def test_get_pfam_annotations():
    a = PfamAnnotation("1j85", "A")
    assert str(a.get()) == "[('1j85', 'A', 'SpoU_methylase', 'SpoU rRNA Methylase family', 'PF00588', 'K31', 100)]"
    a = PfamAnnotation("3bdg", "B")
    result = a.get()
    assert len(result) == 1
    assert result[0][:6] == ('3bdg', 'B', 'Alk_phosphatase', 'Alkaline phosphatase', 'PF00245', 'S31')
    assert 0 <= result[0][6] <= 100


def test_pfam_region_alignment():
    assert PfamAnnotation.get_region('XXACDEFGYY', 'ACDEFG') == ('ACDEFG', 3, 8)


def test_fetch_pdb_info():
    import tempfile
    a = fetchPDBinfo("1uak", "A")
    assert str(a.getCalfaBreaks()) == '[162, 172]'

@pytest.mark.skip(reason="CDD server is down")
def test_cddf_parser():
    d = getCath("1j85", "A")
    assert d[0]['DOMAIN'] == '1j85A00' and d[0]['CATHCODE'] == '3.40.1280.10' and d[0]['CLASS'] == 'Alpha Beta'
