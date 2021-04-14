import glob
import os
import csv
import shutil
import tarfile
import urllib.request, urllib.error, urllib.parse
import datetime
import json
import warnings
from os import makedirs, path
from itertools import count, groupby

import requests
from Bio.PDB import MMCIFIO, Select, PDBParser, PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
from mysolr import Solr
import logging, sys

from requests import HTTPError

from spyprot.common import _gunzip

PDBE_SOLR_URL = "https://www.ebi.ac.uk/pdbe/search/pdb"
UNLIMITED_ROWS = 10000000
DEBUG = False

ALLOWED_HETATM = ['MSE', 'ORN', 'PCA', 'DGL']

SEQ_CODE = {"ALA": 'A',
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


class SpyprotWarning(Warning):
    """Spyprot warning."""
    pass


class ProteinFile:
    '''
       Base classs for download of PDB/CIF files from PDB Database and filter by chain and filter out HOH
       Supports also PDB Bundles when there are many subchains for a given protein
    '''

    def __init__(self, dir, pdbcode, chain=None, atom='CA', filter_by_atom=None, preserve_seqid=False):
        self.pdbcode = pdbcode.lower()
        self.chain = chain
        self.dir = dir
        self.atom = atom
        self.filter_by_atom = filter_by_atom
        self.preserve_seqid = preserve_seqid
        self.pdbdata = []
        self.crystlen = -1
        self.missing_array = []
        self.missing = []
        self.seq_idx = []
        self.seq_len = -1
        ###
        self.parser = None
        self.structure = None
        self.data_file = None
        self.io = None
        self.ext = ''
        try:
            makedirs(self.dir)
        except OSError as e:
            pass

    def _get_structure(self):
        if not self.structure:
            self.structure = self._get_parser().get_structure(self.pdbcode, self.data_file)
        return self.structure

    def _get_parser(self):
        return self.parser

    @staticmethod
    def _get_file(url, file):
        if not path.isfile(file) or path.getsize(file) == 0:
            response = urllib.request.urlopen(url)
            if file.endswith('gz'):
                with open(file, 'wb') as myfile:
                    myfile.write(response.read())
            else:
                html = response.read().decode("UTF-8")
                with open(file, 'w') as myfile:
                    myfile.write(html)

    def get_pdb_data(self):
        first_resid = -1
        for ch in self._get_structure().get_chains():
            if ch.get_id() == self.chain:
                for residue in ch.get_residues():
                    if self.atom in residue.child_dict:
                        first_resid = self._get_seqid_from_residue(residue)
                        break
        prev_seqid = first_resid - 9999999999999
        for ch in self._get_structure().get_chains():
            if ch.get_id() == self.chain:
                for residue in ch.get_residues():
                    if self.atom in residue.child_dict:
                        # Filter HETATM other than those in ALLOWED_HETATM
                        if not residue.id[0].strip() == '' or (
                                residue.id[0].strip().startswith('H_') and residue.resname not in ALLOWED_HETATM):
                            continue
                        seq_id = self._get_seqid_from_residue(residue)
                        if seq_id == prev_seqid:
                            prev_seqid = seq_id
                            print(prev_seqid)
                            continue
                        elif prev_seqid > 0 and seq_id - first_resid != prev_seqid:
                            warnings.warn(
                                'Gap in residues for %s %s between: %i and %i ' % (
                                self.pdbcode, self.chain, prev_seqid, seq_id), SpyprotWarning,
                            )
                        coords = [residue.child_dict[self.atom].coord[i] for i in
                                  range(len(residue.child_dict[self.atom].coord))]
                        line = [seq_id] + coords
                        if None not in line:
                            if self.preserve_seqid:
                                new_seqid = line[0]
                            else:
                                new_seqid = line[0] - first_resid + 1
                        prev_seqid = new_seqid
                        self.pdbdata.append(
                            [new_seqid] + line[1:] + [residue.resname, residue.child_dict[self.atom].bfactor])
        self.crystlen = len(self.pdbdata)
        if self.crystlen == 0:
            return []
        seq = [x[0] for x in self.pdbdata]
        self.missing = [x for x in range(seq[0], seq[-1] + 1) if x not in seq]
        self.missing_array = self.missing
        self.missing = ProteinFile.convertToRanges(self.missing)
        self.seq_idx = seq
        self.seq_len = seq[-1]
        return self.pdbdata

    # Overridden in subclasses
    def _get_seqid_from_residue(self, residue):
        pass

    def get_ca_len(self):
        return self.crystlen

    def get_seq_len(self):
        return self.seq_len

    def get_missing(self):
        return self.missing

    def get_missing_array(self):
        return self.missing_array

    def get_xyz_list(self):
        if not self.pdbdata:
            self.get_pdb_data()
        return [el[:4] for el in self.pdbdata]

    def save_xyz(self, output):
        o = ''
        for l in self.get_xyz_list():
            o += "%4d %8.3f %8.3f %8.3f\n" % (l[0], float(l[1]), float(l[2]), float(l[3]))
        if output != '':
            with open(output, "w") as f:
                f.write(o)

    def save_pdb(self, output):
        if not self.pdbdata:
            self.get_pdb_data()
        o = ''
        for l in self.pdbdata:
            line = (l[0], l[4], self.chain, l[0], float(l[1]), float(l[2]), float(l[3]), float(l[5]))
            if self.atom == 'CA':
                o += "ATOM%7d  CA%5s%2s%4d%12.3f%8.3f%8.3f  1.00%6.2f           C\n" % line
            else:
                o += "ATOM%7d  C3'%4s%2s%4d%12.3f%8.3f%8.3f  1.00%6.2f           C\n" % line
        if output != '':
            with open(output, "w") as f:
                f.write(o)

    def get_chains(self):
        return [ch.get_id() for ch in self._get_structure().get_chains()]

    def get_breaks(self):
        # check chain breaks
        eps = 4.2 * 4.2
        brk = []
        o = self.get_xyz_list()
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

    def get_first_residue_id(self):
        for ch in self._get_structure().get_chains():
            if ch.get_id() == self.chain:
                for residue in ch.get_residues():
                    return residue.get_id()[1]
        return 0

    def get_residue_list(self):
        residues = {}
        for ch in self._get_structure().get_chains():
            if ch.get_id() == self.chain:
                for residue in ch.get_residues():
                    if residue.__dict__['resname'] != 'HOH':
                        residues[residue.get_id()] = residue.__dict__['resname']
        return residues

    def get_seq_one_letter_code(self):
        seq = ''
        k = list(SEQ_CODE.keys())
        for ch in self._get_structure().get_chains():
            if ch.get_id() == self.chain:
                for residue in ch.get_residues():
                    if (residue.id[0].strip() == '' and self.atom in residue.child_dict.keys()) or (
                            residue.id[0].strip().startswith('H_') and residue.resname in ALLOWED_HETATM) and (
                            residue.child_dict[self.atom].altloc.strip() == '' or residue.child_dict[
                        self.atom].altloc.strip() == 'A'):
                        if self.atom == 'CA':
                            if residue.resname in k:
                                seq += SEQ_CODE[residue.resname]
                            else:
                                seq += "X"
                        else:
                            seq += residue.resname
        return seq

    def filter_by_chain(self):
        self.io.set_structure(self._get_structure())
        self.out_file = self.data_file.replace(self.ext, "_" + self.chain + self.ext)
        self.io.save(self.out_file, ChainAndResidueSelect(self.chain))
        del self.io

    def filter_by_atom(self):
        self.io.set_structure(self._get_structure())
        self.out_file = self.out_file.replace(self.ext, "_" + self.atom + self.ext)
        self.io.save(self.out_file, ChainAndAtomSelect(chain=self.chain, atom=self.atom))
        del self.io

    @staticmethod
    def convertToRanges(L):
        G = (list(x) for _, x in groupby(L, lambda x, c=count(): next(c) - x))
        return ["-".join(map(str, (g[0], g[-1])[:len(g)])) for g in G]


class ChainAndResidueSelect(Select):
    def __init__(self, chain, model=0, residue_out="HOH"):
        self.chain = chain
        self.model = model
        self.residue_out = residue_out

    def accept_model(self, model):
        if model is None or (model is not None and model.id == self.model):
            return True
        else:
            return False

    def accept_chain(self, chain):
        if chain is None or (chain is not None and chain.id == self.chain):
            return True
        else:
            return False

    def accept_residue(self, residue):
        if residue.resname != self.residue_out and str(residue._id[0]).strip() == '':
            return True
        else:
            return False


class ChainSelect(Select):
    def __init__(self, chain, model=1, residue_out="HOH"):
        self.chain = chain
        self.model = model
        self.residue_out = residue_out

    def accept_chain(self, chain):
        if chain is None or (chain is not None and chain.id == self.chain):
            return True
        else:
            return False

    def accept_residue(self, residue):
        if residue.resname != self.residue_out and str(residue._id[0]).strip() == '':
            return True
        else:
            return False


class ChainAndAtomSelect(Select):
    def __init__(self, chain, model=1, residue_out="HOH", atom='CA'):
        self.chain = chain
        self.model = model
        self.residue_out = residue_out
        self.atom = atom

    def accept_chain(self, chain):
        if chain is None or (chain is not None and chain.id == self.chain):
            return True
        else:
            return False

    def accept_residue(self, residue):
        if residue.resname != self.residue_out and str(residue._id[0]).strip() == '':
            return True
        else:
            return False

    def accept_atom(self, atom):
        if atom is None or (atom is not None and atom.get_name() == self.atom):
            return True
        else:
            return False


class PdbFile(ProteinFile):
    '''
       Download PDB files from RCSB PDB Database and filter by chain
       Supports also PDB Bundles when there are many subchains for a given protein

       example:
       PdbFile("/tmp", "1j85", "A").download()

       Parameters
       ==========
       path: string - where to store PDB file
       code: string - PDB ID
       chain: string - if empty return full PDB file or list of translated PDB files in case of PDB Bundles

    '''

    def __init__(self, dir, pdbcode, chain=None, atom='CA', filter_by_atom=None, preserve_seqid=False):
        super().__init__(dir=dir, pdbcode=pdbcode, chain=chain, atom=atom, filter_by_atom=filter_by_atom,
                         preserve_seqid=preserve_seqid)
        self.data_file = path.join(self.dir, self.pdbcode + '.pdb')
        self.parser = PDBParser()
        self.io = PDBIO()
        self.ext = '.pdb'

    def download(self):
        try:
            ProteinFile._get_file('https://files.rcsb.org/view/' + self.pdbcode.upper() + '.pdb', self.data_file)
            if self.chain is not None:
                self.filter_by_chain()
        except urllib.error.HTTPError as e:
            warnings.warn(
                "%s %s ..trying to download from PDB Bundle"
                % (self.pdbcode, str(self.chain)), SpyprotWarning,
            )
            response = urllib.request.urlopen(
                'https://files.rcsb.org/pub/pdb/compatible/pdb_bundle/' + self.pdbcode.lower()[
                                                                          1:3] + '/' + self.pdbcode.lower() + '/' + self.pdbcode.lower() + '-pdb-bundle.tar.gz')
            tar = tarfile.open(fileobj=response, mode="r|gz")
            tar.extractall(self.dir)
            tar.close()
            mapFile, mapChain, mapFileChain = self.parsePdbBundleChainIdFile(
                self.dir + '/' + self.pdbcode.lower() + '-chain-id-mapping.txt')
            # chain is set
            if self.chain:
                newChain = mapChain.get(self.chain)
                pdbBundleFile = self.dir + '/' + mapFile.get(self.chain)
                if newChain != self.chain:
                    self.parsePdbAndTranslateChain(pdbBundleFile, self.chain, newChain)
                else:
                    shutil.move(self.dir + '/' + mapFile.get(self.chain), self.data_file)
                    self.filter_by_chain()
                fileList = glob.glob(os.path.join(self.dir, self.pdbcode + "*bundle*.pdb"))
                for filePath in fileList:
                    try:
                        os.remove(filePath)
                    except OSError:
                        warnings.warn("Error while deleting file", SpyprotWarning)
            else:
                self.out_files = []
                file_num = 1
                for file in sorted(mapFileChain.keys()):
                    mapChainFiltered = {key: mapChain[key] for key in mapFileChain.get(file)}
                    self.parsePdbAndTranslateAllChains(self.dir + '/' + file, mapChainFiltered, file_num)
                    file_num = file_num + 1
                self.data_file = self.out_files
        if self.filter_by_atom:
            self.filter_by_atom()
        return self.data_file

    @staticmethod
    def parsePdbBundleChainIdFile(chainFile):
        with open(chainFile, encoding='utf-8') as fp:
            line = fp.readline()
            cnt = 1
            files = []
            mapChain = {}
            mapFileChain = {}
            mapFile = {}
            actualFile = ''
            while line:
                line = line.strip().rstrip()
                cnt += 1
                if line.find('pdb-bundle') >= 0:
                    actualFile = line[:-1]
                    files.append(actualFile)
                    line = fp.readline()
                    continue
                elif actualFile != '' and line != '':
                    mapping = line.split()
                    key = mapping[1].strip()
                    val = mapping[0].strip()
                    mapChain[key] = val
                    mapFile[key] = actualFile
                    if not actualFile in mapFileChain.keys():
                        mapFileChain[actualFile] = []
                    mapFileChain[actualFile].append(key)
                line = fp.readline()
            return mapFile, mapChain, mapFileChain

    # Remapping chain for PDB bundles with many subchains
    def parsePdbAndTranslateChain(self, pdbFileIn, chain, newChain):
        # print chain + '->' + newChain
        self.data_file = self.data_file.replace(".pdb", "_" + self.chain + ".pdb")
        with open(pdbFileIn, "r", encoding='utf-8') as infile, open(self.data_file, "w", encoding='utf-8') as outfile:
            reader = csv.reader(infile)
            for i, line in enumerate(reader):
                if line[0].find('ATOM') == 0 or line[0].find('HETATM') == 0:
                    newLine = list(line[0])
                    if newLine[21] == newChain:
                        if len(chain) > 2:
                            newLine += chain
                            newLine[21] = "%"
                        elif len(chain) > 1:
                            newLine[20] = chain[0]
                            newLine[21] = chain[1]
                        else:
                            newLine[21] = chain
                        outfile.write("".join(newLine) + "\n")
                else:
                    outfile.write(line[0] + "\n")

        # Remapping chain for PDB bundles with many subchains

    def parsePdbAndTranslateAllChains(self, pdbFileIn, mapChainFiltered, file_num):
        # print chain + '->' + newChain
        outf = self.data_file.replace(".pdb", "_bundle_" + str(file_num) + ".pdb")
        self.out_files.append(outf)
        with open(pdbFileIn, "r", encoding='utf-8') as infile, open(outf, "w", encoding='utf-8') as outfile:
            reader = csv.reader(infile)
            for i, line in enumerate(reader):
                if line[0].find('ATOM') == 0 or line[0].find('HETATM') == 0:
                    newLine = list(line[0])
                    for chain in mapChainFiltered.keys():
                        newChain = mapChainFiltered.get(chain)
                        if newLine[21] == newChain:
                            if len(chain) > 2:
                                newLine += chain
                                newLine[21] = "%"
                            elif len(chain) > 1:
                                newLine[20] = chain[0]
                                newLine[21] = chain[1]
                            else:
                                newLine[21] = chain
                            outfile.write("".join(newLine) + "\n")
                else:
                    outfile.write(line[0] + "\n")

    # Unreliable
    # def get_pdb_creation_date(self):
    #    return self.get_parser().header['deposition_date']
    def _get_seqid_from_residue(self, residue):
        return int(residue.get_id()[1])


class MMCIFfile(ProteinFile):
    def __init__(self, dir, pdbcode, chain=None, atom='CA', filter_by_atom=None, preserve_seqid=False):
        super().__init__(dir=dir, pdbcode=pdbcode, chain=chain, atom=atom, filter_by_atom=filter_by_atom,
                         preserve_seqid=preserve_seqid)
        self.data_file = path.join(self.dir, self.pdbcode + ".cif")
        self.parser = MMCIFParser()
        self.io = MMCIFIO()
        self.ext = '.cif'

    def download_only(self):
        try:
            if not path.isfile(self.data_file) or path.getsize(self.data_file) == 0:
                ProteinFile._get_file('http://www.ebi.ac.uk/pdbe/entry-files/download/' + self.pdbcode + '.cif.gz',
                                      self.data_file + '.gz')
                _gunzip(self.data_file + '.gz')
                os.remove(self.data_file + '.gz')
        except Exception as e:
            ProteinFile._get_file('http://www.ebi.ac.uk/pdbe/entry-files/download/' + self.pdbcode + '.cif',
                                  self.data_file)

    def download(self):
        self.download_only()
        if self.chain is not None:
            self.filter_by_chain()
        else:
            self.out_file = self.data_file
        if self.filter_by_atom:
            self.filter_by_atom()
        return self.out_file

    def _get_seqid_from_residue(self, residue):
        return int(self._get_parser()._mmcif_dict["_atom_site.label_seq_id"][
                       residue.child_dict[self.atom].serial_number - 1])

    def get_pdb_creation_date(self):
        # return self._get_parser()._mmcif_dict['_pdbx_database_status.recvd_initial_deposition_date'][0]
        return self.get_par_from_dict('_pdbx_database_status.recvd_initial_deposition_date')

    def get_meta_pubmed(self):
        pubmed_id = self.get_par_from_dict('_citation.pdbx_database_id_PubMed')
        doi = self.get_par_from_dict('_citation.pdbx_database_id_DOI')
        title = self.get_par_from_dict('_citation.title')
        desc = self.get_par_from_dict('_struct.title')
        src = self.get_par_from_dict('_entity_src_gen.pdbx_gene_src_scientific_name')
        key = self.get_par_from_dict('_struct_keywords.pdbx_keywords')
        molecutag = self.get_par_from_dict('_entity.pdbx_description')
        return (doi, pubmed_id, desc, title, src, key, molecutag)

    def get_par_from_dict(self, par):
        # return self._get_parser()._mmcif_dict[par][0] if par in self._get_parser()._mmcif_dict else None
        if par in self._get_parser()._mmcif_dict:
            res = self._get_parser()._mmcif_dict[par]
            if type(res) == str:
                return res
            elif type(res) == list:
                return res[0]
            else:
                return 'coooooooo'
        else:
            return None

    def get_seq_one_letter_code_can(self):
        model_id = -1
        for model in self._get_structure().child_list:
            for ch in model.child_list:
                if ch.id == self.chain:
                    model_id = model.id
        return self.parser._mmcif_dict['_entity_poly.pdbx_seq_one_letter_code'][model_id].replace('\n', '')


class SearchException(Exception):
    pass


class SequenceException(Exception):
    pass


class PDBeSolrSearch:
    def __init__(self):
        self.solr = Solr(PDBE_SOLR_URL, version=4)
        logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                            format='LOG|%(asctime)s|%(levelname)s  %(message)s', datefmt='%d-%b-%Y %H:%M:%S')
        logging.getLogger("requests").setLevel(logging.WARNING)

    @staticmethod
    def join_with_AND(selectors):
        return " AND ".join(
            ["%s:%s" % (k, v) for k, v in selectors]
        )

    @staticmethod
    def join_with_OR(selectors):
        return " OR ".join(
            ["%s:%s" % (k, v) for k, v in selectors]
        )

    def exec_query(self, field_list, query_details):
        try:
            query = {"rows": UNLIMITED_ROWS, "fl": field_list, "q": self.join_with_AND(query_details)}
            response = self.solr.search(**query)
            documents = response.documents
            if DEBUG:
                logging.getLogger().debug("Found %d matching entities in %d entries." % (
                    len(documents), len({rd["pdb_id"] for rd in documents})))
            return documents
        except Exception as e:
            raise SearchException('%s: error getting query response from: %s for: %s\n %s' % (
                self.__class__.__name__, PDBE_SOLR_URL, str(query), str(e)))

    def get(self):
        return self.results


class IdenticalChains(PDBeSolrSearch):
    '''
       Find identical chains to a given one

       example:
       a = IdenticalChains("2jlo",chain="A").get()

       Parameters
       ==========
       pdbcode: string - PDB ID
    '''

    def __init__(self, pdbcode, chain):
        super().__init__()

        self.results = []
        documents = self.exec_query("pdb_id,entity_id,chain_id,assembly_composition", [('pdb_id', pdbcode)])

        for i in range(len(documents)):
            chain_id = documents[i]['chain_id']
            if chain in chain_id:
                self.results = sorted(chain_id)
                break


class SimilarChains(PDBeSolrSearch):
    '''
       Find similar chains to a given one with sequence identity given as parameter

       example:
       a = SimilarChains("2jlo",chain="A",identity=90).get()

       Parameters
       ==========
       pdbcode: string - PDB ID
       chain: string
       seq: string - sequence to compare against - instead od pdbcode and chain
       identity: int - (percentage) of sequence identity
    '''

    def __init__(self, pdb=None, chain='A', seq=None, identity=40):
        super().__init__()
        self.chain = chain
        self.pdb = pdb.lower() if pdb else None
        self.identity = float(identity / 100)
        self.identifiers = []
        self.results = []

        try:
            self.seq = seq if seq else self.get_seq()
            self.get_similar()
            self.translate_enity_ids_to_chains()
        except (urllib.error.URLError, HTTPError, SequenceException) as he:
            raise SearchException('Problem looking for similar chains to: ' + pdb + ' ' + chain + ': ' + str(he))

    def get_seq(self):
        documents = self.exec_query("pdb_id,entity_id,chain_id,molecule_sequence", [('pdb_id', self.pdb)])
        for i in range(len(documents)):
            chain_id = documents[i]['chain_id']
            if self.chain in chain_id:
                return documents[i]['molecule_sequence']
        raise SequenceException('SimilarChains: Error getting sequence from PDBe for: %s, %s' % (self.pdb, self.chain))

    def get_similar(self):
        query = """
{
  "query": {
    "type": "terminal",
    "service": "sequence",
    "parameters": {
      "evalue_cutoff": 1,
      "identity_cutoff": %s,
      "target": "pdb_protein_sequence",
      "value": "%s"
    }
  },
  "return_type": "polymer_entity",
  "request_options": {
    "pager": {
      "start": 0,
      "rows": 100
    },
    "scoring_strategy": "sequence",
    "sort": [
      {
        "sort_by": "score",
        "direction": "desc"
      }
    ]
  }
}               
        """ % (self.identity, self.seq)
        url = "https://search.rcsb.org/rcsbsearch/v1/query?json={0}".format(urllib.parse.quote(query))
        try:
            response = requests.get(url)
            if response.status_code == 204 and response.text == '':
                return
            if response.status_code != 200:
                raise SequenceException(response.text)
            result = response.json()
            if 'result_set' in result:
                for el in result['result_set']:
                    self.identifiers.append(el['identifier'])
        except Exception as er:
            raise SearchException('SimilarChains: Error in response from RCSB search for: %s, %s URL:\n%s\n%s' % (
                self.pdb, self.chain, urllib.parse.unquote(url), str(er)))

    def translate_enity_ids_to_chains(self):
        if not self.identifiers:
            return
        ident_list = self.identifiers[0]
        for entr_ent in self.identifiers[1:]:
            ident_list += ' OR ' + entr_ent.lower()
        documents = self.exec_query("pdb_id,entity_id,chain_id,assembly_composition",
                                    [('entry_entity', '(' + ident_list + ')')])
        for i in range(len(documents)):
            pid = documents[i]['pdb_id']
            chain_id = documents[i]['chain_id']
            self.results.append((pid.upper(), sorted(chain_id)[0]))

    def get(self):
        return sorted(self.results)


class ReleasedPDBs(PDBeSolrSearch):
    '''
          Download list of pdb ids released in range from_date - to_date (or in day from_date).

          Return empty list if None found
          example:
          a = ReleasedProteins("2020-10-10", "2020-11-10").get()

          Parameters
          ==========
          from_date: datetime or string in format YYYY-MM-DD
          to_date: datetime or string in format YYYY-MM-DD
          uniq_chains: if True return a list of tuples containing (PDBCode, Chain), else return just PDBCodes.
       '''

    def __init__(self, from_date, to_date='', uniq_chains=True, only_prot=True, only_rna=False):
        super().__init__()
        if type(from_date) is datetime.date:
            from_date = from_date.strftime("%Y-%m-%d")
        if to_date == '':
            to_date = from_date
        if type(to_date) is datetime.date:
            to_date = to_date.strftime("%Y-%m-%d")

        if only_rna:
            molecule_type = 'RNA'
        elif not only_prot:
            molecule_type = '(Protein OR RNA)'
        else:
            molecule_type = 'Protein'

        self.results = []
        documents = self.exec_query("pdb_id,entity_id,chain_id,molecule_type",
                                    [
                                        ('release_date', '[' + from_date + 'T00:00:00Z TO ' + to_date + 'T23:59:59Z]'),
                                        ('molecule_type', molecule_type)
                                    ])
        for i in range(len(documents)):
            pid = documents[i]['pdb_id']
            chain_id = documents[i]['chain_id']
            if uniq_chains:
                self.results.append((pid, sorted(chain_id)[0]))
            else:
                if pid not in self.results:
                    self.results.append(pid)
        self.results = sorted(self.results)


class UniqueChains(PDBeSolrSearch):
    '''
       Find a list of unique chains for a given PDB id

       example:
       a = UniqChains("2jlo").get()

       Parameters
       ==========
       pdbcode: string - PDB ID
    '''

    def __init__(self, pdbcode, only_prot=True, only_rna=False):
        super().__init__()

        if only_rna:
            molecule_type = 'RNA'
        elif not only_prot:
            molecule_type = '(Protein OR RNA)'
        else:
            molecule_type = 'Protein'
        self.results = []
        documents = self.exec_query("pdb_id,entity_id,chain_id,assembly_composition,molecule_type",
                                    [
                                        ('pdb_id', pdbcode),
                                        ('molecule_type', molecule_type)
                                    ])
        for i in range(len(documents)):
            chain_id = documents[i]['chain_id']
            self.results.append(sorted(chain_id)[0])


class PdbMetaData(PDBeSolrSearch):
    '''
       Retrieve meta data for a given structure

       example:
       a = IdenticalChains("2jlo",chain="A").get()

       Parameters
       ==========
       pdbcode: string - PDB ID
    '''

    def __init__(self, pdbcode):
        super().__init__()

        self.results = []
        documents = self.exec_query(
            "pdb_id,entity_id, chain_id, assembly_composition, deposition_date, number_of_protein_chains, "
            "pubmed_id, citation_doi, citation_title, title, mutation, mutation_type, assembly_composition, enzyme_systematic_name, enzyme_name, molecule_sequence",
            [('pdb_id', pdbcode)])

        for i in range(len(documents)):
            chain_id = documents[i]['chain_id']
            self.results = sorted(chain_id)
            # break

            # doi = doi[0].text if doi else None
            # molecutag = molecutag[0].text if molecutag else None
            # key = key[0].text.capitalize() if key else None
            # src = src[0].text if src else None
            # pubmed = pubmed[0].text if pubmed else None
            # title = title[0].text if title else None
            # desc = desc[0].text.capitalize() if desc else None
            #
            # return (doi, pubmed, desc, title, src, key, molecutag)
