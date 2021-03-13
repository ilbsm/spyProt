import pytest

from spyprot.fetchChainInfo import SimilarChains, UniqueChains, IdenticalChains, ReleasedPDBs, SearchException
from datetime import datetime


def test_fetchChainInfo_IdenticalChains():
    f = IdenticalChains('6wm4', 'Y')
    assert str(f.get()) == "['X', 'Y', 'Z']"


def test_fetchChainInfo_UniqueChains():
    f = UniqueChains('6wm4')
    assert str(f.get()) == "['R', 'O', 'H', 'K', 'S', 'T', 'A', 'D', 'X', 'G', 'N', 'U', 'V', '0', '1', 'Q', 'P']"


def test_fetchChainInfo_UniqueChains_2():
    f = UniqueChains('2jlo')
    assert f.get() == []


def test_fetchChainInfo_UniqueChains_3():
    f = UniqueChains('6zj3', only_rna=True)
    assert str(
        f.get()) == "['S1', 'S2', 'S3', 'S4', 'S5', 'LA', 'LB', 'LC', 'LD', 'LE', 'LF', 'LG', 'LH', 'LI', 'LJ', 'LK', 'LL', 'LM', 'LN', 'LO']"
    f = UniqueChains('6zj3', only_prot=True)
    assert len(f.get()) == 78
    f = UniqueChains('6zj3', only_prot=False, only_rna=False)
    assert len(f.get()) == 98


def test_fetchChainInfo_SimilarChains():
    sim = SimilarChains(pdb='1j85', chain='A')
    assert ('3N4J', 'A') in sim.get()
    assert len(sim.get()) >= 12
    sim = SimilarChains('6lt7', 'A', identity=30)
    assert len(sim.get()) >= 3
    sim = SimilarChains(pdb='7css', chain='A')
    assert len(sim.get()) == 0


def test_fetchChainInfo_SimilarChains_fail():
    with pytest.raises(SearchException):
        sim = SimilarChains(seq='AFHAGAGOANBAG')
        sim.get()
    with pytest.raises(SearchException):
        sim = SimilarChains(pdb='1j85aa', chain='A')
        sim.get()
    with pytest.raises(SearchException):
        sim = SimilarChains(pdb='7djq', chain='C')
        sim.get()


def test_fetchChainInfo_ReleasedProteins():
    from_date = datetime(2020, 11, 10).date()
    to_date = datetime(2020, 11, 18).date()
    res = ReleasedPDBs(from_date, to_date).get()
    assert str(res).startswith(
        "[('5rw2', 'A'), ('5rw3', 'A'), ('5rw4', 'A'), ('5rw5', 'A'), ('5rw6', 'A'), ('5rw7', 'A'), ('5rw8', 'A'), ('5rw9', 'A'), ('5rwa', 'A'),")


def test_fetchChainInfo_ReleasedProteins_2():
    from_date = "2020-11-18"
    res = ReleasedPDBs(from_date).get()
    assert len(res) == 490
    assert str(res).startswith(
        "[('6hpj', 'B'), ('6kml', 'A'), ('6kml', 'B'), ('6kmq', 'A'), ('6kmq', 'B'), ('6l9k', 'A'), ('6l9k', 'Q'), ('6l9l', 'A'), ('6l9l', 'B'), ('6l9l', 'C'),")
    res = ReleasedPDBs(from_date, only_rna=True).get()
    assert str(res) == "[('6hpj', 'A'), ('6vem', 'A'), ('6wbr', 'B'), ('6wc0', 'B'), ('6wvj', 'R')]"
    res = ReleasedPDBs(from_date, only_rna=False, only_prot=False).get()
    assert len(res) == 495
