from spyprot import fetchPDBinfo, getCoordinates
from spyprot.fetchChainInfo import PdbFile, MMCIFfile, PdbMetaData
import tempfile
import os
import filecmp


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


def test_cif_parser_vs_xml():
    pdb_chains = ['4v7m_B3', '6az3_B', '6az3_e', '1j85_A', '6ki1_A']
    dir = os.path.join(tempfile.gettempdir(), 'test_cif_xml')
    os.makedirs(dir, exist_ok=True)
    from Bio import BiopythonWarning
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        for el in pdb_chains:
            pdb = el.split('_')[0]
            chain = el.split('_')[1]
            m = MMCIFfile(dir, pdb, chain)
            m.download()
            m.get_pdb_data()
            m.save_xyz(dir + '/' + el + '_cif.xyz')
            m.save_pdb(dir + '/' + el + '_cif.pdb')
            c = getCoordinates(pdbcode=pdb, work_dir=dir)
            c.getCalfa(chain=chain, output=dir + '/' + el + '_xml.xyz')
            c.getCalfaPdbFormat(chain=chain, output=dir + '/' + el + '_xml.pdb')
            assert filecmp.cmp(dir + '/' + el + '_cif.xyz', dir + '/' + el + '_xml.xyz', shallow=False) == True
            assert filecmp.cmp(dir + '/' + el + '_cif.pdb', dir + '/' + el + '_xml.pdb', shallow=False) == True

            z = fetchPDBinfo(pdbcode=pdb, chain=chain, work_dir=dir)
            ch_r = z.getOrderedChains()
            ch_l = m.get_chains()
            if ch_l!=ch_r:
                print('Problem ' + pdb + ': ' + str(len(ch_l)) + '!=' + str(len(ch_r)) + '\n' + str(ch_l) + '\n' + str(ch_r))
            assert z.getMissingArray() == m.get_missing_array()
            assert z.getMissing() == m.get_missing()
            assert z.getCAlen() == m.get_ca_len()
            assert z.getSeqOneLetterCode() == m.get_seq_one_letter_code()
            assert z.getSeqLength() == m.get_seq_len()
            assert z.getPdbCreationDate() == m.get_pdb_creation_date()
            pm = m.get_meta_pubmed()
            zpm = z.getPubtitlePubmed()
            for i in range(len(zpm)):
                l = str(zpm[i]).replace(' ','').lower()
                r = str(pm[i]).replace(' ','').lower()
                if l!=r:
                    print('problem ' + pdb + ' ' + chain + ': ' + l + '!=' +r)
                assert l==r


def test_pdb_parser_vs_xml():
    # Will not work correctly for bundles
    pdb_chains = ['6az3_B', '6az3_e', '1j85_A', '6ki1_A']
    dir = os.path.join(tempfile.gettempdir(), 'test_pdb_xml')
    os.makedirs(dir, exist_ok=True)
    from Bio import BiopythonWarning
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        for el in pdb_chains:
            pdb = el.split('_')[0]
            chain = el.split('_')[1]
            m = PdbFile(dir, pdb, chain)
            m.download()
            m.get_pdb_data()
            m.save_xyz(dir + '/' + el + '_cif.xyz')
            m.save_pdb(dir + '/' + el + '_cif.pdb')
            c = getCoordinates(pdbcode=pdb, work_dir=dir)
            c.getCalfa(chain=chain, output=dir + '/' + el + '_xml.xyz')
            c.getCalfaPdbFormat(chain=chain, output=dir + '/' + el + '_xml.pdb')
            assert filecmp.cmp(dir + '/' + el + '_cif.xyz', dir + '/' + el + '_xml.xyz', shallow=False) == True
            assert filecmp.cmp(dir + '/' + el + '_cif.pdb', dir + '/' + el + '_xml.pdb', shallow=False) == True

            z = fetchPDBinfo(pdbcode=pdb, chain=chain, work_dir=dir)
            ch_r = z.getOrderedChains()
            ch_l = m.get_chains()
            if ch_l!=ch_r:
                print('Problem ' + pdb + ': ' + str(len(ch_l)) + '!=' + str(len(ch_r)) + '\n' + str(ch_l) + '\n' + str(ch_r))
            assert z.getMissingArray() == m.get_missing_array()
            assert z.getMissing() == m.get_missing()
            assert z.getCAlen() == m.get_ca_len()
            assert z.getSeqOneLetterCode() == m.get_seq_one_letter_code()
            assert z.getSeqLength() == m.get_seq_len()



#dir = os.path.join(tempfile.gettempdir(), 'test_pdb')
#PdbFile(dir, "4v7m", "B3").download()
#p = PdbFile(dir, "1j85", "A")
#p.download()
#p.get_pdb_creation_date()
#p.get_pdb_data()
#test_download_pdb_files()
#test_pdb_parser_vs_xml()
#test_cif_parser_vs_xml()

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
# dir = os.path.join(tempfile.gettempdir(), 'test_pdb')
# #pdb = '1x0g'
# pdb = '6ki1'
# chain = 'A'
#
# m = MMCIFfile(dir, pdb, chain)
# m.download()
# print(m.get_xyz_list())
# #print(m.get_xyz_list(preserve_seqid=True))
# print(m.get_breaks())
#
# print(m.get_seq_one_letter_code())
# print(m.get_meta_pubmed())
#
# print(m.get_missing())
# broken_res = m.get_breaks()
# missing_array = m.get_missing_array()
# str_len = m.get_ca_len()
# seq_len = m.get_seq_len()
# missing = ", ".join(m.get_missing())
# starting_index = m.get_first_residue_id()
# print(broken_res)
# print(missing_array)
# print(str_len)
# print(seq_len)
# print(missing)
# print(starting_index)
#
# #pdb = '6sj9'
# #m = MMCIFfile(dir, pdb, 'A')
# #m.download()
# #print(m.get_xyz_list())
# #print(m.get_breaks())
#
#
#
# z = fetchPDBinfo(pdbcode=pdb, chain=chain, work_dir=dir)
# print(z.getChains())
# print(z.getPubtitlePubmed())
# print(z.getSeqOneLetterCode())
#
# #z.setChain(chain)
# print(z.getCalfaBreaks())
#
# print(z.getMissing())
# broken_res = z.getCalfaBreaks()
# missing_array = z.getMissingArray()
# sequence = z.getSeqOneLetterCode()
# str_len = int(z.getCAlen())
# seq_len = z.getSeqLength()
# missing = ", ".join(z.getMissing())
# date = z.getPdbCreationDate()
# starting_index = z.getFirstResidueIndex()
# print(broken_res)
# print(missing_array)
# print(str_len)
# print(seq_len)
# print(missing)
# print(starting_index)
#
# c = getCoordinates(pdb, work_dir=dir)
# c.getCalfaPdbFormat(chain=chain, output=dir + '/' + pdb + '__.pdb', preserve_seqid=True)
#
#
# pdbm = PdbMetaData(pdb)
# print(pdbm)
