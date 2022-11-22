import pytest

from spyprot.fetchChainInfo import SimilarChains, UniqueChains, IdenticalChains, ReleasedPDBs, SearchException, \
    UniprotInfo, IdenticalChainsAndEntityId, PdbSequence, UniprotSearch
from datetime import datetime


def test_fetchChainInfo_IdenticalChains():
    f = IdenticalChains('6wm4', 'Y')
    assert str(f.get()) == "['X', 'Y', 'Z']"


def test_fetchChainInfo_IdenticalChainsAndEntityId():
    f = IdenticalChainsAndEntityId('7dzp', 'A')
    ent_id, chains = f.get()
    assert ent_id == None
    assert chains == None
    f = IdenticalChainsAndEntityId('6wm4', 'Y')
    ent_id, chains = f.get()
    assert int(ent_id) == 9
    assert str(chains) == "['X', 'Y', 'Z']"


def test_fetchChainInfo_UniqueChains():
    f = UniqueChains('6wm4')
    assert str(sorted(f.get())) == str(
        sorted(['R', 'O', 'H', 'K', 'S', 'T', 'A', 'D', 'X', 'G', 'N', 'U', 'V', '0', '1', 'Q', 'P']))


def test_fetchChainInfo_UniqueChains_2():
    f = UniqueChains('2jlo')
    assert f.get() == []


def test_fetchChainInfo_UniqueChains_3():
    f = UniqueChains('6zj3', only_rna=True)
    assert str(sorted(f.get())) == str(sorted(
        ['S1', 'S2', 'S3', 'S4', 'S5', 'LA', 'LB', 'LC', 'LD', 'LE', 'LF', 'LG', 'LH', 'LI', 'LJ', 'LK', 'LL', 'LM',
         'LN', 'LO']))
    f = UniqueChains('6zj3', only_prot=True)
    assert len(f.get()) == 78
    f = UniqueChains('6zj3', only_prot=False, only_rna=False)
    assert len(f.get()) == 98


def test_fetchChainInfo_SimilarChains():
    sim = SimilarChains(pdb='1j85', chain='A')
    assert str(
        sim.get()) == "[('1J85', 'A'), ('3N4J', 'A'), ('3N4K', 'A'), ('4JAK', 'A'), ('4JAL', 'A'), ('4KDZ', 'A'), ('6QH8', 'A'), ('7D51', 'A'), ('7E3S', 'A')]"
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
    assert len(res) == 487
    assert str(res).startswith(
        "[('6hpj', 'B'), ('6kml', 'A'), ('6kml', 'B'), ('6kmq', 'A'), ('6kmq', 'B'), ('6l9k', 'A'), ('6l9k', 'Q'), ('6l9l', 'A'), ('6l9l', 'B'), ('6l9l', 'C'),")
    res = ReleasedPDBs(from_date, only_rna=True).get()
    assert str(res) == "[('6hpj', 'A'), ('6vem', 'A'), ('6wbr', 'B'), ('6wc0', 'B'), ('6wvj', 'R')]"
    res = ReleasedPDBs(from_date, only_rna=False, only_prot=False).get()
    assert len(res) == 492


def test_UniprotInfo_sequence():
    s = UniprotInfo('Q57812').get_sequence()
    assert s == 'MPLVGFMKEKKRATFYLYKNIDGRKLRYLLHKLENVENVDIDTLRRAIEAEKKYKRSITLTEEEEVIIQRLGKSANLLLNCELVKLDEGERA'


def test_UniprotInfo_title():
    res = UniprotInfo('C6KTF9').get_title()
    assert res == 'Erythrocyte membrane protein 1, PfEMP1'

def test_UniprotInfo_pdbs():
    res = UniprotInfo('Q15149').get_pdbs()
    assert str(res) == "[('1MB8', 'A', '175-403'), ('2N03', 'A', '4403-4606'), ('2ODU', 'A', '410-640'), ('2ODV', 'A', '410-640'), ('3F7P', 'A/B', '175-403'), ('3PDY', 'A/B', '653-858'), ('3PE0', 'A/B', '750-1028'), ('4GDO', 'A/B/C/D/E/F', '1492-1530'), ('4Q58', 'A/B', '175-400'), ('4Q59', 'A/B', '175-400'), ('5J1F', 'A/B', '860-928,999-1116'), ('5J1G', 'A/B', '1114-1343'), ('5J1H', 'A/B', '860-928,999-1116'), ('5J1I', 'A/B', '1114-1482')]"
    res = UniprotInfo('Q03661').get_pdbs()
    assert str(res) == "[('6QSZ', 'B/D/F/H/J/L/N/P', '1443-1458')]"
    res = UniprotInfo('P42945').get_pdbs()
    assert str(res) == "[('5WLC', 'LM', '1-1769'), ('5WYJ', 'AE', '1-808'), ('5WYK', 'AE', '1-1769'), ('6KE6', 'AE', '1-1769'), ('6LQP', 'AE', '1-1769'), ('6LQQ', 'AE', '1-1769'), ('6LQR', 'AE', '1-1769'), ('6LQS', 'AE', '1-1769'), ('6LQT', 'AE', '1-1769'), ('6LQU', 'AE', '1-1769'), ('6LQV', 'AE', '1-1769'), ('6ND4', 'M', '1-1769'), ('6ZQA', 'UJ', '1-1769'), ('6ZQB', 'UJ', '1-1769'), ('6ZQC', 'UJ', '1-1769'), ('6ZQD', 'UJ', '1-1769'), ('6ZQE', 'UJ', '1-1769'), ('7AJT', 'UJ', '1-1769'), ('7AJU', 'UJ', '1-1769'), ('7D4I', 'AE', '1-1769'), ('7D5S', 'AE', '1-1769'), ('7D5T', 'AE', '1-1769'), ('7D63', 'AE', '1-1769')]"


def test_UniprotInfo_length():
    lengths = {'O60281': 2723, 'Q15149': 4684, 'Q03661': 1658, 'P42945': 1769}
    for uniid in lengths.keys():
        res = UniprotInfo(uniid).get_overall_length()
        assert res == lengths[uniid]


def test_fetchChainInfo_PdbSequence():
    seqs = ['101m', '102l', '102m']
    res = PdbSequence(seqs).get()
    assert str(res) == "{'102m': 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKAGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG', '102l': 'MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL', '101m': 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG'}"


def test_UniprotSearch():
    uni = UniprotSearch(
        ['xref_pdb', 'xref_pfam', 'organism_name', 'organism_id', 'length', 'protein_name', 'reviewed', 'gene_names', 'lineage', 'lineage_ids'], accessions=['A0A3B5KNW3'])
    res = uni.get()
    assert str(res) == "{'A0A3B5KNW3': ['', '', 'Takifugu rubripes (Japanese pufferfish) (Fugu rubripes)', '31033', '543', 'Interleukin-1 receptor accessory protein-like', 'unreviewed', 'LOC105417273', 'Takifugu (genus), Tetraodontidae (family), Tetradontoidea (superfamily), Tetraodontoidei (suborder), Tetraodontiformes (order), Eupercaria (no rank), Percomorphaceae (no rank), Euacanthomorphacea (no rank), Acanthomorphata (no rank), Ctenosquamata (no rank), Eurypterygia (no rank), Neoteleostei (no rank), Euteleosteomorpha (cohort), Clupeocephala (no rank), Osteoglossocephalai (no rank), Teleostei (infraclass), Neopterygii (subclass), Actinopteri (class), Actinopterygii (superclass), Euteleostomi (no rank), Teleostomi (no rank), Gnathostomata (no rank), Vertebrata (no rank), Craniata (subphylum), Chordata (phylum), Deuterostomia (no rank), Bilateria (no rank), Eumetazoa (no rank), Metazoa (kingdom), Opisthokonta (no rank), Eukaryota (superkingdom), cellular organisms (no rank)', '31032 (genus), 31031 (family), 32517 (superfamily), 31028 (suborder), 31022 (order), 1489922 (no rank), 1489872 (no rank), 123369 (no rank), 123368 (no rank), 123367 (no rank), 123366 (no rank), 123365 (no rank), 1489388 (cohort), 186625 (no rank), 1489341 (no rank), 32443 (infraclass), 41665 (subclass), 186623 (class), 7898 (superclass), 117571 (no rank), 117570 (no rank), 7776 (no rank), 7742 (no rank), 89593 (subphylum), 7711 (phylum), 33511 (no rank), 33213 (no rank), 6072 (no rank), 33208 (kingdom), 33154 (no rank), 2759 (superkingdom), 131567 (no rank)']}"


def test_UniprotSearch_Array():
    uni = UniprotSearch(
        ['xref_pdb', 'xref_pfam', 'organism_name', 'organism_id', 'length', 'protein_name', 'reviewed', 'gene_names', 'lineage', 'lineage_ids'], accessions=['A0A3B5KNW3', 'A0A815ZZQ3'])
    res = uni.get(as_dict=False)
    assert str(res) == "[['A0A815ZZQ3', '', 'PF13517;', 'Adineta ricciae (Rotifer)', '249248', '773', 'Hypothetical protein', 'unreviewed', 'XAT740_LOCUS46458', 'Adineta (genus), Adinetidae (family), Adinetida (order), Bdelloidea (subclass), Eurotatoria (class), Rotifera (phylum), Gnathifera (no rank), Spiralia (no rank), Protostomia (no rank), Bilateria (no rank), Eumetazoa (no rank), Metazoa (kingdom), Opisthokonta (no rank), Eukaryota (superkingdom), cellular organisms (no rank)', '104781 (genus), 104780 (family), 104779 (order), 44578 (subclass), 2816136 (class), 10190 (phylum), 2697496 (no rank), 2697495 (no rank), 33317 (no rank), 33213 (no rank), 6072 (no rank), 33208 (kingdom), 33154 (no rank), 2759 (superkingdom), 131567 (no rank)'], ['A0A3B5KNW3', '', '', 'Takifugu rubripes (Japanese pufferfish) (Fugu rubripes)', '31033', '543', 'Interleukin-1 receptor accessory protein-like', 'unreviewed', 'LOC105417273', 'Takifugu (genus), Tetraodontidae (family), Tetradontoidea (superfamily), Tetraodontoidei (suborder), Tetraodontiformes (order), Eupercaria (no rank), Percomorphaceae (no rank), Euacanthomorphacea (no rank), Acanthomorphata (no rank), Ctenosquamata (no rank), Eurypterygia (no rank), Neoteleostei (no rank), Euteleosteomorpha (cohort), Clupeocephala (no rank), Osteoglossocephalai (no rank), Teleostei (infraclass), Neopterygii (subclass), Actinopteri (class), Actinopterygii (superclass), Euteleostomi (no rank), Teleostomi (no rank), Gnathostomata (no rank), Vertebrata (no rank), Craniata (subphylum), Chordata (phylum), Deuterostomia (no rank), Bilateria (no rank), Eumetazoa (no rank), Metazoa (kingdom), Opisthokonta (no rank), Eukaryota (superkingdom), cellular organisms (no rank)', '31032 (genus), 31031 (family), 32517 (superfamily), 31028 (suborder), 31022 (order), 1489922 (no rank), 1489872 (no rank), 123369 (no rank), 123368 (no rank), 123367 (no rank), 123366 (no rank), 123365 (no rank), 1489388 (cohort), 186625 (no rank), 1489341 (no rank), 32443 (infraclass), 41665 (subclass), 186623 (class), 7898 (superclass), 117571 (no rank), 117570 (no rank), 7776 (no rank), 7742 (no rank), 89593 (subphylum), 7711 (phylum), 33511 (no rank), 33213 (no rank), 6072 (no rank), 33208 (kingdom), 33154 (no rank), 2759 (superkingdom), 131567 (no rank)']]"
