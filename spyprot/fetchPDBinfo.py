#!/usr/bin/env python
import urllib2 as urllib
from itertools import count, groupby
from os import path
import tempfile

from lxml import etree


def convertToRanges(L):
    G = (list(x) for _, x in groupby(L, lambda x, c=count(): next(c) - x))
    return ["-".join(map(str, (g[0], g[-1])[:len(g)])) for g in G]


class XmlPdbParser:
    def __init__(self, pdbcode, chain, work_dir=tempfile.gettempdir(), localfile=False, atom="CA"):
        self.pdb = pdbcode
        self.chain = chain
        self.atom = atom
        self.work_dir = work_dir
        self.localfile = localfile
        self.xmlgzfile = self.download_if_not_exist()
        self.root = etree.parse(self.xmlgzfile, etree.XMLParser(ns_clean=True)).getroot()
        self.NS = self.root.nsmap

    def download_if_not_exist(self):
        pdb_xml_file = path.join(self.work_dir, self.pdb + '.xml.gz')
        if not path.isfile(pdb_xml_file):
            f = urllib.urlopen('https://www.rcsb.org/pdb/files/' + self.pdb + '.xml.gz')
            fw = open(pdb_xml_file, "wb")
            fw.write(f.read())
            fw.close()
            f.close()
        return pdb_xml_file


class getCoordinates(XmlPdbParser):
    def __init__(self, pdbcode, work_dir=tempfile.gettempdir(), localfile=False, atom="CA"):
        XmlPdbParser.__init__(self, pdbcode, None, work_dir, localfile, atom=atom)
        l = self.root.xpath(
            "//PDBx:atom_siteCategory/PDBx:atom_site[PDBx:auth_atom_id=\"" + atom + "\"][PDBx:pdbx_PDB_model_num='1'][PDBx:group_PDB='ATOM']/PDBx:auth_asym_id",
            namespaces=self.NS)
        self.chains = set()
        self.chains.update(e.text for e in l)

    def getChainIndexes(self):
        return self.chains

    def getCalfa(self, chain, output='', preserve_seqid=False):
        l = self.root.xpath(
            "//PDBx:atom_siteCategory/PDBx:atom_site[not(PDBx:label_alt_id/text()) or PDBx:label_alt_id='A'][PDBx:auth_atom_id=\"" + self.atom + "\"][PDBx:pdbx_PDB_model_num='1'][PDBx:group_PDB='ATOM' or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='MSE') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='ORN') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='PCA') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='DGL')][PDBx:auth_asym_id='" + chain + "']",
            namespaces=self.NS)
        # mozna by pewnie zrobic inaczej or hetatom and (xxx or xxx or xxx ....)
        o = ''
        # r = []
        prev_seqid = 999999999999999
        first_resid = -1
        for e in l:
            d = e.find("PDBx:label_seq_id", namespaces=self.NS)
            if d.text:
                first_resid = int(d.text)
                break
        prev_seqid = first_resid - 9999999999999
        for e in l:
            x = e.find("PDBx:Cartn_x", namespaces=self.NS)
            y = e.find("PDBx:Cartn_y", namespaces=self.NS)
            z = e.find("PDBx:Cartn_z", namespaces=self.NS)
            seqid = e.find("PDBx:label_seq_id", namespaces=self.NS)
            if not seqid.text:
                continue
            if int(seqid.text) == prev_seqid:
                prev_seqid = int(seqid.text)
                print(prev_seqid)
                continue
            prev_seqid = int(seqid.text)
            l = [seqid, x, y, z]
            if None not in l:
                if preserve_seqid:
                    new_seqid = int(l[0].text)
                else:
                    new_seqid = int(l[0].text) - first_resid + 1
                o += "%4d %8.3f %8.3f %8.3f\n" % (new_seqid, float(l[1].text), float(l[2].text), float(l[3].text))
        if output != '':
            f = open(output, "w")
            f.write(o)
            f.close()
        else:
            return o

    def getCalfaPdbFormat(self, chain, output='', preserve_seqid=False):
        l = self.root.xpath(
            "//PDBx:atom_siteCategory/PDBx:atom_site[not(PDBx:label_alt_id/text()) or PDBx:label_alt_id='A'][PDBx:auth_atom_id=\"" + self.atom + "\"][PDBx:pdbx_PDB_model_num='1'][PDBx:group_PDB='ATOM' or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='MSE') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='ORN') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='PCA') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='DGL')][PDBx:auth_asym_id='" + chain + "']",
            namespaces=self.NS)
        o = ''
        r = []
        prev_seqid = 999999999999999
        first_resid = -1
        for e in l:
            d = e.find("PDBx:label_seq_id", namespaces=self.NS)
            if d.text:
                first_resid = int(d.text)
                break
        prev_seqid = first_resid - 9999999999999

        for e in l:
            x = e.find("PDBx:Cartn_x", namespaces=self.NS)
            y = e.find("PDBx:Cartn_y", namespaces=self.NS)
            z = e.find("PDBx:Cartn_z", namespaces=self.NS)
            seqname = e.find("PDBx:label_comp_id", namespaces=self.NS)
            seqid = e.find("PDBx:label_seq_id", namespaces=self.NS)
            true_seqid = e.find("PDBx:auth_seq_id", namespaces=self.NS)
            if not seqid.text:
                continue
            if preserve_seqid and not true_seqid.text:
                continue
            if int(seqid.text) == prev_seqid:
                prev_seqid = int(seqid.text)
                continue
            prev_seqid = int(seqid.text)
            bfactor = e.find("PDBx:B_iso_or_equiv", namespaces=self.NS)
            l = (seqid, seqname, chain, seqid, x, y, z, bfactor)
            if None not in l:
                if preserve_seqid:
                    new_seqid = int(true_seqid.text)
                else:
                    new_seqid = int(seqid.text) - first_resid + 1
                l = (new_seqid, str(seqname.text), str(chain), new_seqid, float(x.text), float(y.text), float(z.text),
                     float(bfactor.text))
                if self.atom == 'CA':
                    o += "ATOM%7d  CA%5s%2s%4d%12.3f%8.3f%8.3f  1.00%6.2f           C\n" % l
                else:
                    o += "ATOM%7d  C3'%4s%2s%4d%12.3f%8.3f%8.3f  1.00%6.2f           C\n" % l
        if output != '':
            f = open(output, "w")
            f.write(o)
            f.close()
        else:
            return o


class fetchPDBinfo(XmlPdbParser):
    def __init__(self, pdbcode, chain='A', work_dir=tempfile.gettempdir(), localfile=False, atom="CA"):
        XmlPdbParser.__init__(self, pdbcode, chain, work_dir, localfile, atom=atom)
        self.codification = {"ALA": 'A',
                             "CYS": 'C',
                             "ASP": 'D',
                             "GLU": 'E',
                             "PHE": 'F',
                             "GLY": 'G',
                             "HIS": 'H',
                             "ILE": 'I',
                             "LYS": 'K',
                             "LEU": 'L',
                             "MET": 'M',
                             "MSE": 'M',
                             "ASN": 'N',
                             "PYL": 'O',
                             "PRO": 'P',
                             "GLN": 'Q',
                             "ARG": 'R',
                             "SER": 'S',
                             "THR": 'T',
                             "SEC": 'U',
                             "VAL": 'V',
                             "TRP": 'W',
                             "5HP": 'E',
                             "ABA": 'A',
                             "AIB": 'A',
                             "BMT": 'T',
                             "CEA": 'C',
                             "CGU": 'E',
                             "CME": 'C',
                             "CRO": 'X',
                             "CSD": 'C',
                             "CSO": 'C',
                             "CSS": 'C',
                             "CSW": 'C',
                             "CSX": 'C',
                             "CXM": 'M',
                             "DAL": 'A',
                             "DAR": 'R',
                             "DCY": 'C',
                             "DGL": 'E',
                             "DGN": 'Q',
                             "DHI": 'H',
                             "DIL": 'I',
                             "DIV": 'V',
                             "DLE": 'L',
                             "DLY": 'K',
                             "DPN": 'F',
                             "DPR": 'P',
                             "DSG": 'N',
                             "DSN": 'S',
                             "DSP": 'D',
                             "DTH": 'T',
                             "DTR": 'X',
                             "DTY": 'Y',
                             "DVA": 'V',
                             "FME": 'M',
                             "HYP": 'P',
                             "KCX": 'K',
                             "LLP": 'K',
                             "MLE": 'L',
                             "MVA": 'V',
                             "NLE": 'L',
                             "OCS": 'C',
                             "ORN": 'A',
                             "PCA": 'E',
                             "PTR": 'Y',
                             "SAR": 'G',
                             "SEP": 'S',
                             "STY": 'Y',
                             "TPO": 'T',
                             "TPQ": 'F',
                             "TYS": 'Y',
                             "TYR": 'Y'}
        self._getCAmissing()

    def getCalfaBreaks(self, preserve_seqid=False):
        l = self.root.xpath(
            "//PDBx:atom_siteCategory/PDBx:atom_site[not(PDBx:label_alt_id/text()) or PDBx:label_alt_id='A'][PDBx:auth_atom_id=\"" + self.atom + "\"][PDBx:pdbx_PDB_model_num='1'][PDBx:group_PDB='ATOM' or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='MSE') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='ORN') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='PCA') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='DGL')][PDBx:auth_asym_id='" + self.chain + "']",
            namespaces=self.NS)
        # mozna by pewnie zrobic inaczej or hetatom and (xxx or xxx or xxx ....)
        o = []
        # r = []
        prev_seqid = 999999999999999
        first_resid = -1
        for e in l:
            d = e.find("PDBx:label_seq_id", namespaces=self.NS)
            if d is not None and d.text:
                first_resid = int(d.text)
                break
        prev_seqid = first_resid - 9999999999999
        for e in l:
            x = e.find("PDBx:Cartn_x", namespaces=self.NS)
            y = e.find("PDBx:Cartn_y", namespaces=self.NS)
            z = e.find("PDBx:Cartn_z", namespaces=self.NS)
            seqid = e.find("PDBx:label_seq_id", namespaces=self.NS)
            if not seqid.text:
                continue
            if int(seqid.text) == prev_seqid:
                prev_seqid = int(seqid.text)
                print(prev_seqid)
                continue
            prev_seqid = int(seqid.text)
            l = [seqid, x, y, z]
            if None not in l:
                if preserve_seqid:
                    new_seqid = int(l[0].text)
                else:
                    new_seqid = int(l[0].text) - first_resid + 1
                o.append((new_seqid, float(l[1].text), float(l[2].text), float(l[3].text)))
        # check chain breaks
        eps = 4.2 * 4.2
        brk = []
        for i in range(1, len(o)):
            i1 = o[i - 1][0]
            i2 = o[i][0]
            x1 = o[i - 1][1]
            y1 = o[i - 1][2]
            z1 = o[i - 1][3]

            x2 = o[i][1]
            y2 = o[i][2]
            z2 = o[i][3]

            d = (x1 - x2) ** 2 + (y2 - y1) ** 2 + (z1 - z2) ** 2
            if d > eps:
                brk.append(i1)
                brk.append(i2)
        return brk

    def getFirstResidueIndex(self):
        return 1

    def getTrueFirstResidueIndex(self):
        l = self.root.xpath(
            "//PDBx:atom_siteCategory/PDBx:atom_site[not(PDBx:label_alt_id/text()) or PDBx:label_alt_id='A'][PDBx:auth_atom_id=\"" + self.atom + "\"][PDBx:pdbx_PDB_model_num='1'][PDBx:group_PDB='ATOM' or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='MSE') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='ORN') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='PCA') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='DGL')][PDBx:auth_asym_id='" + self.chain + "']/PDBx:auth_seq_id",
            namespaces=self.NS)
        for e in l:
            if e.text:
                # print e.text
                return int(e.text)
        return 1  # could be not true anymore

    def getSeqOneLetterCode(self):
        l = self.root.xpath(
            "//PDBx:atom_siteCategory/PDBx:atom_site[PDBx:auth_atom_id=\"" + self.atom + "\"][PDBx:auth_asym_id='" + self.chain + "'][PDBx:pdbx_PDB_model_num='1'][PDBx:group_PDB='ATOM' or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='MSE') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='ORN') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='PCA') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='DGL')][not(PDBx:label_alt_id/text()) or PDBx:label_alt_id='A']/PDBx:label_comp_id",
            namespaces=self.NS)
        k = list(self.codification.keys())
        seq = ""
        for e in l:
            s = e.text
            if self.atom == 'CA':
                if s in k:
                    seq += self.codification[s]
                else:
                    seq += "X"
            else:
                seq += s
        return seq

    def getMissingArray(self):
        return self.missing_array

    def _getCAmissing(self):
        seq = []
        l = self.root.xpath(
            "//PDBx:atom_siteCategory/PDBx:atom_site[PDBx:auth_atom_id=\"" + self.atom + "\"][PDBx:auth_asym_id='" + self.chain + "'][PDBx:pdbx_PDB_model_num='1'][PDBx:group_PDB='ATOM' or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='MSE') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='ORN') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='PCA') or (PDBx:group_PDB='HETATM' and PDBx:label_comp_id='DGL')][not(PDBx:label_alt_id/text()) or PDBx:label_alt_id='A']/PDBx:label_seq_id",
            namespaces=self.NS)
        for e in l:
            if e.text: seq.append(int(e.text))
        if len(seq) > 0:  # TODO make better function for getting missing - it fails e.g. in case of 3bpd_C
            first = seq[0]
        for e in range(len(seq)):
            seq[e] = seq[e] - first + 1  # przenumerowanie od 1 z zachowaniem przerw

        self.crystlen = len(seq)
        if len(seq) > 0:
            self.missing = [x for x in range(seq[0], seq[-1] + 1) if x not in seq]
            self.missing_array = self.missing
            self.missing = convertToRanges(self.missing)
            self.seq_idx = seq
            self.seq_len = seq[-1]
        else:
            self.missing = []
            self.missing_array = self.missing
            self.seq_len = 0

    def getSeqLength(self):
        return self.seq_len

    def getPdbCreationDate(self):
        self.date = self.root.xpath(
            "//PDBx:pdbx_database_statusCategory/PDBx:pdbx_database_status/PDBx:recvd_initial_deposition_date",
            namespaces=self.NS)[0].text
        #        self.date = self.root.xpath("//PDBx:database_PDB_revCategory/PDBx:database_PDB_rev/PDBx:date_original",namespaces=self.NS)[0].text
        try:
            return self.date
        except:
            return None

    def getCAlen(self):
        return self.crystlen

    def getMissing(self):
        return self.missing

    def getPubtitlePubmed(self):
        pubmed = self.root.xpath("//PDBx:citationCategory/PDBx:citation[1]/PDBx:pdbx_database_id_PubMed",
                                 namespaces=self.NS)
        doi = self.root.xpath("//PDBx:citationCategory/PDBx:citation[1]/PDBx:pdbx_database_id_DOI", namespaces=self.NS)
        title = self.root.xpath("//PDBx:citationCategory/PDBx:citation[1]/PDBx:title", namespaces=self.NS)
        desc = self.root.xpath("//PDBx:structCategory/PDBx:struct/PDBx:title", namespaces=self.NS)
        src = self.root.xpath("//PDBx:entity_src_genCategory/PDBx:entity_src_gen/PDBx:pdbx_gene_src_scientific_name",
                              namespaces=self.NS)
        key = self.root.xpath("//PDBx:struct_keywordsCategory/PDBx:struct_keywords/PDBx:pdbx_keywords",
                              namespaces=self.NS)
        molecutag = self.root.xpath("//PDBx:entityCategory/PDBx:entity/PDBx:pdbx_description", namespaces=self.NS)

        doi = doi[0].text if doi else None
        molecutag = molecutag[0].text if molecutag else None
        key = key[0].text.capitalize() if key else None
        src = src[0].text if src else None
        pubmed = pubmed[0].text if pubmed else None
        title = title[0].text if title else None
        desc = desc[0].text.capitalize() if desc else None

        return (doi, pubmed, desc, title, src, key, molecutag)


if __name__ == "__main__":
    #    pdb=argv[1].lower()
    r = fetchPDBinfo('6CZR', '1a', '/tmp', localfile=False, atom="C3'")
    # r = fetchPDBinfo('/tmp','1j85', 'A', localfile=False, atom="CA")
    print(r.getSeqLength())
    print(r.getCAlen())
    # b = getCoordinates("/tmp", "1uak")
    # b.getCalfaPdbFormat("A", output=path.join("/tmp", "t.pdb"))
    # a = fetchPDBinfo("/tmp", "1uak","A")
    # print(a.getCalfaBreaks())
    # sequence = a.getSeqOneLetterCode()
    # print("seq",sequence)
    # str_len = int(a.getCAlen())
    # print("str_len",str_len)
    # seq_len = a.getSeqLength()
    # print("seq_len",seq_len)
    # missing = ", ".join(a.getMissing())
    # print("miss",missing)
    # date = a.getPdbCreationDate()
    # print("date",date)
    # structure_start = a.getFirstResidueIndex()
    # print("start",structure_start)
    # doi, pubmed, pubtitle, prottitle, source, keyword, molname =  a.getPubtitlePubmed()
    # print("rest","\n".join([doi, pubmed, pubtitle, prottitle, source, keyword, molname]))
