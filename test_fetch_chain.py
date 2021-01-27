from spyprot.fetchChainInfo import SimilarChains, UniqueChains, IdenticalChains


def test_fetchChainInfo_IdenticalChains():
    f = IdenticalChains('6wm4', 'Y')
    assert str(f.get()) == "['X', 'Y', 'Z']"


def test_fetchChainInfo_UniqueChains():
    f = UniqueChains('6wm4')
    assert str(f.get()) == "['R', 'O', 'T', 'A', 'D', 'X', 'G', 'N', 'U', 'V', '0', '1', 'Q', 'P', 'H', 'K', 'S']"


def test_fetchChainInfo_UniqueChains_2():
    f = UniqueChains('2jlo')
    assert f.get() == []


def test_fetchChainInfo_SimilarChains():
    sim = SimilarChains(pdb='1j85', chain='A')
    assert str(sim.get()) == "[('1J85', 'A'), ('3E5Y', 'A'), ('3N4J', 'A'), ('3N4K', 'A'), ('4JAK', 'A'), ('4JAL', 'A'), ('4KDZ', 'A'), ('4KGN', 'A'), ('4PZK', 'A'), ('5CO4', 'A'), ('6AHW', 'A'), ('6QH8', 'A'), ('6QKV', 'A')]"

