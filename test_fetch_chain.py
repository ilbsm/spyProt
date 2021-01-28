from spyprot.fetchChainInfo import SimilarChains, UniqueChains, IdenticalChains, ReleasedPDBs
from datetime import date, datetime, timedelta


def test_fetchChainInfo_IdenticalChains():
    f = IdenticalChains('6wm4', 'Y')
    assert str(f.get()) == "['X', 'Y', 'Z']"


def test_fetchChainInfo_UniqueChains():
    f = UniqueChains('6wm4')
    assert str(f.get()) == "['R', 'O', 'T', 'A', 'D', 'X', 'G', 'N', 'U', 'V', '0', '1', 'Q', 'P', 'H', 'K', 'S']"


def test_fetchChainInfo_UniqueChains_2():
    f = UniqueChains('2jlo')
    assert f.get() == []


def test_fetchChainInfo_UniqueChains_3():
    f = UniqueChains('6zj3', only_rna=True)
    assert str(f.get()) == "['LA', 'LB', 'LG', 'LF', 'LJ', 'LK', 'LL', 'LM', 'LN', 'LO', 'S2', 'S3', 'S4', 'LC', 'LD', 'LE', 'LH', 'LI', 'S5', 'S1']"
    f = UniqueChains('6zj3', only_prot=True)
    assert len(f.get())==78
    f = UniqueChains('6zj3', only_prot=False, only_rna=False)
    assert len(f.get()) == 98


def test_fetchChainInfo_SimilarChains():
    sim = SimilarChains(pdb='1j85', chain='A')
    assert str(sim.get()) == "[('1J85', 'A'), ('3E5Y', 'A'), ('3N4J', 'A'), ('3N4K', 'A'), ('4JAK', 'A'), ('4JAL', 'A'), ('4KDZ', 'A'), ('4KGN', 'A'), ('4PZK', 'A'), ('5CO4', 'A'), ('6AHW', 'A'), ('6QH8', 'A'), ('6QKV', 'A')]"


def test_fetchChainInfo_ReleasedProteins():
    from_date = datetime(2020, 11, 10).date()
    to_date = datetime(2020, 11, 18).date()
    res = ReleasedPDBs(from_date, to_date).get()
    assert str(res).startswith("[('6y9d', 'A'), ('6wv7', 'A'), ('6wyg', 'A'), ('6yve', 'AAA'), ('6z2x', 'C'), ('6z2x', 'E'),")


def test_fetchChainInfo_ReleasedProteins_2():
    from_date = "2020-11-18"
    res = ReleasedPDBs(from_date).get()
    assert len(res) == 490
    assert str(res).startswith("[('6y9d', 'A'), ('6yve', 'AAA'), ('6zj6', 'AAA'), ('6xif', 'A'), ('6xif', 'B'), ('6xif', 'I'), ('6xox', 'A'),")
    res = ReleasedPDBs(from_date, only_rna=True).get()
    assert str(res) == "[('6wvj', 'R'), ('6vem', 'A'), ('6wbr', 'B'), ('6hpj', 'A'), ('6wc0', 'B')]"
    res = ReleasedPDBs(from_date, only_rna=False, only_prot=False).get()
    assert len(res) == 495
