from getAnnotations import getPfamAnnotation, getAnnotations, getECAnnotation


def test_get_annotations():
    a = getAnnotations("1j85", "A")
    assert str(a.get()[0]) =="['1j85', 'A', '2', '142', 'PF00588.18', 'SpoU_methylase', 'SpoU rRNA Methylase family', '4.3E-32']"


def test_get_ec_annotations():
    a = getECAnnotation("1cak", "A")
    assert str(a.get()) == "[('4.2.1.1', 'Carbonic anhydrase.')]"


def test_get_pfam_annotations():
    a = getPfamAnnotation("1j85", "A")
    assert str(a.get()) == "[('1j85', 'A', 'SpoU_methylase', 'SpoU rRNA Methylase family', 'PF00588')]"


if __name__ == '__main__':
    test_get_annotations()
    test_get_ec_annotations()
    test_get_pfam_annotations()