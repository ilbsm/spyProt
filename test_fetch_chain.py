from spyprot.fetchChainInfo import SimilarChains, UniqueChains

def test_fetchChainInfo_FetchUniqueChains():
    f = UniqueChains('6wm4')
    assert str(f.get()) == "['R', 'O', 'T', 'B', 'D', 'Y', 'G', 'N', 'U', 'V', '0', '1', 'Q', 'P', 'H', 'K', 'S']"


def test_fetchChainInfo_FetchUniqueChains2():
    f = UniqueChains('2jlo')
    assert f.get() == []


def test_fetchChainInfo_FetchSimilarChains():
    sim = SimilarChains(pdb='1j85', chain='A')
    assert str(sim.get()) == "[('6QKV', 'A'), ('6AHW', 'A'), ('5CO4', 'B'), ('4PZK', 'A'), ('4KGN', 'H'), ('4KDZ', 'A'), ('4JAK', 'A'), ('4JAL', 'A'), ('3N4J', 'A'), ('3E5Y', 'A'), ('1J85', 'A'), ('3N4K', 'A'), ('6QH8', 'D')]"

