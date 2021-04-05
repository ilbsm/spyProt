from spyprot import fetchPDBinfo
from spyprot.fetchChainInfo import PdbFile, MMCIFfile, PdbMetaData
import tempfile
import os


def test_download_pdb_files():
    dir = os.path.join(tempfile.gettempdir(), 'test_pdb')
    os.makedirs(dir, exist_ok=True)
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


#test_download_pdb_files()
# dir = os.path.join(tempfile.gettempdir(), 'test_pdb')
# outf = PdbFile(dir, "4v7m").download()
# print(outf)
#
# PdbFile(dir, "6az3", "2", 'C3\'').download()
#
#
# dir = os.path.join(tempfile.gettempdir(), 'test_pdb')
# #p = PdbFile(dir, "7bl1", "AAA").download()
# p = PdbFile(dir, "7bl1").download()
dir = os.path.join(tempfile.gettempdir(), 'test_pdb')
pdb = '1j85'

m = MMCIFfile(dir, pdb, 'A')
m.download()
print(m.get_xyz_list())
print(m.get_breaks())

print(m.get_seq_one_letter_code())
print(m.get_meta_pubmed())

#pdb = '6sj9'
#m = MMCIFfile(dir, pdb, 'A')
#m.download()
#print(m.get_xyz_list())
#print(m.get_breaks())



z = fetchPDBinfo(pdbcode=pdb, work_dir=dir)
print(z.getPubtitlePubmed())
print(z.getSeqOneLetterCode())
print(z.getPdbCreationDate())

print(m.get_pdb_creation_date())

print(z.getChains())
z.setChain('A')
print(z.getCalfaBreaks())

pdbm = PdbMetaData(pdb)
print(pdbm)
