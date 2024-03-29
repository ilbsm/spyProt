import os
import tempfile

from spyprot.ilbsm_database_downloader import ILBSMDatabaseDownloader, AlphaKnotDatabaseDownloader


def test_knotprot_download():
    SEARCH_STRING = 'https://knotprot.cent.uw.edu.pl/browse/?set=True&bridgeType=probab&knotTypes=-31&array=0&raw=1'
    URL_KNOT_MATRIX = 'https://knotprot.cent.uw.edu.pl/static/knot_data/{0}/{1}/{0}_{1}.png'
    URL_KNOTOID_SUB_RAW = 'https://knotprot.cent.uw.edu.pl/static/knot_data//{0}/{1}/{0}_{1}_data.txt'
    URL_CHAIN_XYZ = 'https://knotprot.cent.uw.edu.pl/chains/{0}/{1}/chain.xyz.txt'
    with tempfile.TemporaryDirectory() as out_dir:
        dd = ILBSMDatabaseDownloader(SEARCH_STRING, [URL_KNOT_MATRIX,
                                                     URL_KNOTOID_SUB_RAW,
                                                     URL_CHAIN_XYZ], out_dir, create_separate_dirs=False)
        dd.get_all()
        assert len(os.listdir(out_dir)) >= 96


def test_lassoprot_download():
    SEARCH_STRING = 'https://lassoprot.cent.uw.edu.pl/browse/?lassoType=L3&set=0&bridgeType=ssbridge%2Camide%2Cester%2Cthioester%2Cothers&array=0&raw=1'
    URL_ALL_FILES = 'https://lassoprot.cent.uw.edu.pl/files/{0}_{1}_all.tar.bz2'
    with tempfile.TemporaryDirectory() as out_dir:
        dd = ILBSMDatabaseDownloader(SEARCH_STRING, [URL_ALL_FILES], out_dir, create_separate_dirs=False)
        dd.get_all()
        assert len(os.listdir(out_dir)) >= 53


def test_genus_download():
    SEARCH_STRING = 'https://genus.fuw.edu.pl/browse/?moltag=hydrolase%2Fpeptide&set=True&is_rna=&raw=1'
    URL_BONDS_PROTEIN = 'https://genus.fuw.edu.pl/file/{0}/{1}/{0}_{1}.chimera'
    with tempfile.TemporaryDirectory() as out_dir:
        dd = ILBSMDatabaseDownloader(SEARCH_STRING, [URL_BONDS_PROTEIN], out_dir, create_separate_dirs=False)
        dd.get_all()
        assert len(os.listdir(out_dir)) >= 95


def test_alphaknot_download():
    SEARCH_STRING = 'https://alphaknot.cent.uw.edu.pl/browse/?cats=&v=&adv=MOUSE&organisms=&knotTypes=4_1&knotTypes=5_1&knotTypes=6_1&array=0&raw=1'
    URL_KNOT_MODEL = 'https://alphaknot.cent.uw.edu.pl/model_file/{3}/{0}-F{1}-v{2}.cif'
    with tempfile.TemporaryDirectory() as out_dir:
        dd = AlphaKnotDatabaseDownloader(SEARCH_STRING, [URL_KNOT_MODEL], out_dir, create_separate_dirs=False)
        dd.get_all()
        assert len(os.listdir(out_dir)) >= 42
