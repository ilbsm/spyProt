from spyprot.fetchChainInfo import PdbFile, MMCIFfile
import tempfile
import os


def test_download_pdb_files():
    dir = os.path.join(tempfile.gettempdir(), 'test_pdb')
    try:
        os.makedirs(dir)
    except:
        pass
    p = PdbFile(dir, "1j85", "A").download()
    p = PdbFile(dir, "1k36", "A").download()
    p = PdbFile(dir, "1k36", "A", 'CA').download()
    p = MMCIFfile(dir, "1k36", "A").download()

    # Download and extract chain 5
    p = PdbFile(dir, "6p5i", "5").download()
    # Download and extract only C3' atoms from chain 5
    p = PdbFile(dir, "6p5i", "5", 'C3\'').download()

    m = MMCIFfile(dir, "1k36", "A").download()
    # Download and extract chain 5
    m = MMCIFfile(dir, "6p5i", "5").download()
    # Download and extract only C3' atoms from chain 5

    MMCIFfile(dir, "3bjx", "A").download()
    MMCIFfile(dir, "1et4", "A").download()

    PdbFile(dir, "4v7m", "B3").download()
    # Bundle but not translation
    # PdbFile(dir, "6az3", "2").download()
    PdbFile(dir, "6az3", "2", 'C3\'').download()

# test_download_pdb_files()
