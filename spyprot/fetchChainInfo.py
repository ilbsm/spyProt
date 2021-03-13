import glob
import os
import csv
import shutil
import tarfile
import urllib2 as urllib
import datetime
import json
from os import makedirs, path

import requests
from Bio.PDB import MMCIFIO, Select, PDBParser, PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
from mysolr import Solr
import logging, sys

from requests import HTTPError

PDBE_SOLR_URL = "https://www.ebi.ac.uk/pdbe/search/pdb"
UNLIMITED_ROWS = 10000000
DEBUG = False


class ProteinFile:
    '''
       Base classs for download of PDB/CIF files from PDB Database and filter by chain and filter out HOH
       Supports also PDB Bundles when there are many subchains for a given protein
    '''

    def __init__(self, dir, pdbcode, chain=None, atom=None):
        self.pdbcode = pdbcode
        self.chain = chain
        self.dir = dir
        self.atom = atom


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
       chain: string

    '''

    def download(self):
        self.out_file = self.dir + '/' + self.pdbcode + '.pdb'
        try:
            makedirs(self.dir)
        except OSError as e:
            pass
        try:
            response = urllib.urlopen('https://files.rcsb.org/view/' + self.pdbcode.upper() + '.pdb')
            html = response.read().decode("UTF-8")
            with open(self.out_file, 'w') as myfile:
                myfile.write(html)
            if self.chain is not None:
                self.filter_by_chain()
        except urllib.HTTPError as e:
            print(self.pdbcode + " " + self.chain + " ...trying to download from PDB Bundle")
            response = urllib.urlopen('https://files.rcsb.org/pub/pdb/compatible/pdb_bundle/' + self.pdbcode.lower()[
                                                                                                1:3] + '/' + self.pdbcode.lower() + '/' + self.pdbcode.lower() + '-pdb-bundle.tar.gz')
            tar = tarfile.open(fileobj=response, mode="r|gz")
            tar.extractall(self.dir)
            tar.close()
            mapFile, mapChain = self.parsePdbBundleChainIdFile(
                self.dir + '/' + self.pdbcode.lower() + '-chain-id-mapping.txt')
            newChain = mapChain.get(self.chain)
            pdbBundleFile = self.dir + '/' + mapFile.get(self.chain)
            if newChain != self.chain:
                self.parsePdbAndTranslateChain(pdbBundleFile, self.chain, newChain)
            else:
                shutil.move(self.dir + '/' + mapFile.get(self.chain), self.out_file)
                self.filter_by_chain()
            fileList = glob.glob(os.path.join(self.dir, self.pdbcode + "*bundle*.pdb"))
            for filePath in fileList:
                try:
                    os.remove(filePath)
                except OSError:
                    print("Error while deleting file")
        if self.atom:
            self.filter_by_atom()
        return self.out_file

    @staticmethod
    def parsePdbBundleChainIdFile(chainFile):
        with open(chainFile) as fp:
            line = fp.readline()
            cnt = 1
            files = []
            mapChain = {}
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
                line = fp.readline()
            return mapFile, mapChain

    # Remapping chain for PDB bundles with many subchains
    def parsePdbAndTranslateChain(self, pdbFileIn, chain, newChain):
        # print chain + '->' + newChain
        self.out_file = self.out_file.replace(".pdb", "_" + self.chain + ".pdb")
        with open(pdbFileIn, "r") as infile, open(self.out_file, "w") as outfile:
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

    def filter_by_chain(self):
        parser = PDBParser()
        structure = parser.get_structure(self.pdbcode, self.out_file)
        io = PDBIO()
        io.set_structure(structure)
        self.out_file = self.out_file.replace(".pdb", "_" + self.chain + ".pdb")
        io.save(self.out_file, ChainAndResidueSelect(self.chain))
        del io

    def filter_by_atom(self):
        parser = PDBParser()
        structure = parser.get_structure(self.pdbcode, self.out_file)
        io = PDBIO()
        io.set_structure(structure)
        self.out_file = self.out_file.replace(".pdb", "_" + self.atom + ".pdb")
        io.save(self.out_file, ChainAndAtomSelect(chain=self.chain, atom=self.atom))
        del io


class MMCIFfile(ProteinFile):
    def download_only(self):
        self.cif_file = path.join(self.dir, self.pdbcode + ".cif")
        try:
            makedirs(self.dir)
        except OSError as e:
            pass
        response = urllib.urlopen('http://www.ebi.ac.uk/pdbe/entry-files/download/' + self.pdbcode + '.cif')
        html = response.read().decode("UTF-8")
        with open(self.cif_file, 'w') as myfile:
            myfile.write(html)

    def download(self):
        self.download_only()
        if self.chain is not None:
            self.filter_by_chain()
        else:
            self.out_file = self.cif_file
        if self.atom:
            self.filter_by_atom()
        return self.out_file

    def filter_by_chain(self):
        parser = MMCIFParser()
        structure = parser.get_structure(self.pdbcode, self.cif_file)
        io = MMCIFIO()
        io.set_structure(structure)
        self.out_file = self.cif_file.replace(".cif", "_" + self.chain + ".cif")
        io.save(self.out_file, ChainAndResidueSelect(self.chain))

    def filter_by_atom(self):
        parser = MMCIFParser()
        structure = parser.get_structure(self.pdbcode, self.cif_file)
        io = MMCIFIO()
        io.set_structure(structure)
        self.out_file = self.out_file.replace(".cif", "_" + self.atom + ".cif")
        io.save(self.out_file, ChainAndAtomSelect(chain=self.chain, atom=self.atom))
        del io

    def get_first_residue_id(self):
        parser = MMCIFParser()
        structure = parser.get_structure(self.pdbcode, self.cif_file)
        for ch in structure.get_chains():
            if ch.get_id() == self.chain:
                for residue in ch.get_residues():
                    return residue.get_id()[1]
        return 0

    def get_residue_list(self):
        parser = MMCIFParser()
        structure = parser.get_structure(self.pdbcode, self.cif_file)
        residues = {}
        for ch in structure.get_chains():
            if ch.get_id() == self.chain:
                for residue in ch.get_residues():
                    if residue.__dict__['resname'] != 'HOH':
                        residues[residue.get_id()] = residue.__dict__['resname']
        return residues


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
        # res_decoded = []
        # for res in self.results:
        #     res_decoded.append(str(res))
        # return res_decoded


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
        PDBeSolrSearch.__init__(self)

        self.results = []
        documents = self.exec_query("pdb_id,entity_id,chain_id,assembly_composition", [('pdb_id', pdbcode)])

        for i in range(len(documents)):
            chain_id = documents[i]['chain_id']
            if chain in sorted(chain_id):
                chains_to_add = []
                for ch in chain_id:
                    chains_to_add.append(str(ch))
                self.results = sorted(chains_to_add)
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
        PDBeSolrSearch.__init__(self)
        self.chain = chain
        self.pdb = pdb.lower() if pdb else None
        self.identity = float(identity / 100)
        self.identifiers = []
        self.results = []

        try:
            self.seq = seq if seq else self.get_seq()
            self.get_similar()
            self.translate_enity_ids_to_chains()
        except (urllib.URLError, HTTPError, SequenceException) as he:
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
        url = "https://search.rcsb.org/rcsbsearch/v1/query?json={0}".format(urllib.quote(query))
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
            self.pdb, self.chain, urllib.unquote(url), str(er)))

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
            self.results.append((str(pid.upper()), str(sorted(chain_id)[0])))

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
        PDBeSolrSearch.__init__(self)
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
                self.results.append((str(pid), str(sorted(chain_id)[0])))
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
        PDBeSolrSearch.__init__(self)

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
            self.results.append(str(sorted(chain_id)[0]))
