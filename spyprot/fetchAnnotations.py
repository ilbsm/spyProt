#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Pawel Rubach 2019
# based on:
# Copyright Michal Jamroz, 2014, jamroz@chem.uw.edu.pl
import datetime
import urllib.request, urllib.error, urllib.parse, gzip, re
import os
from os import path
import urllib.request, json

import re

from requests import Session
from Bio import pairwise2
from Bio.Align import substitution_matrices

REFRESH_FILE_INTERVAL = 7 * 24 * 3600  # Refresh every week

PDB_CHAIN_ENZYME = ["ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_enzyme.tsv.gz",
                    "pdb_chain_enzyme.tsv.gz"]
ENZYME_DAT = ["https://ftp.expasy.org/databases/enzyme/enzyme.dat", "enzyme.dat"]

PDB_CHAIN_PFAM = ["ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_pfam.tsv.gz",
                  "pdb_chain_pfam.tsv.gz"]
PFAM_DESC = ["http://pfam.xfam.org/families?output=text", "pdb_chain_pfam_desc"]

DEFAULT_DATA_FILE_PATH = path.join(path.expanduser("~"), ".local", "spyprot")

PDBe = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_to_pfam/"

Uni = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/"


class AnnotationBase:
    '''
       Base class for retrieving annotations for a given PDBID and Chain

       Annotations are parsed from files downloaded from databases and stored locally in:
       $HOME/.local/spyprot/
    '''

    def __init__(self, files, pdb, chain='A', data_file_path=DEFAULT_DATA_FILE_PATH,
                 refresh_file_interval=REFRESH_FILE_INTERVAL):
        self.data_file_path = data_file_path
        if not path.exists(data_file_path):
            os.makedirs(data_file_path)
        self.refresh_file_interval = refresh_file_interval
        self.first_file = self.download_if_not_exist(files[0])
        self.sec_file = self.download_if_not_exist(files[1])
        self.pdb = pdb.lower()
        self.chain = chain
        self.results = []
        self.d = []
        self.identifiers = []

    def get(self):
        return self.results

    @staticmethod
    def touch(path):
        with open(path, 'a'):
            os.utime(path, None)

    @staticmethod
    def download_file(url, file, parent_dir=DEFAULT_DATA_FILE_PATH):
        enz_file = path.join(parent_dir, file)
        AnnotationBase.touch(enz_file + ".flg")
        try:
            f = urllib.request.urlopen(url)
            fw = open(enz_file + "_NEW", "wb")
            fw.write(f.read())
            fw.close()
            if path.isfile(enz_file):
                os.rename(enz_file, enz_file + "_OLD")
            os.rename(enz_file + "_NEW", enz_file)
            if path.isfile(enz_file + "_OLD"):
                os.remove(enz_file + "_OLD")
        except Exception as e:
            print("Problem updating file: " + enz_file + " leaving old one" + str(e))
        if path.isfile(enz_file + ".flg"): os.remove(enz_file + ".flg")

    def download_if_not_exist(self, file_locs):
        enz_file = path.join(self.data_file_path, file_locs[1])
        if not path.isfile(enz_file) or (not path.isfile(enz_file + ".flg") and path.isfile(
                enz_file) and datetime.datetime.now().timestamp() - os.stat(
            enz_file).st_mtime > self.refresh_file_interval):
            print('downloading file: ' + enz_file + " from: " + file_locs[0])
            AnnotationBase.download_file(file_locs[0], file_locs[1], parent_dir=self.data_file_path)
        return enz_file


class ECAnnotation(AnnotationBase):
    '''Retrieve EC annotations for a given PDBID and Chain

       Annotations are parsed from files downloaded from expasy and stored locally in:
       $HOME/.local/spyprot/

       Parameters
       ==========
       pdbcode: string
       chain: string

       Return
       ======
       results: list of tuples of EC Code and EC name
       '''

    def __init__(self, pdb, chain='A'):
        super().__init__([PDB_CHAIN_ENZYME, ENZYME_DAT], pdb, chain)
        q = re.compile(r'^' + self.pdb + '\t' + self.chain + '\t(.*)\t(?P<ec>.*)$')

        with gzip.open(self.first_file) as s:
            for line in s.readlines():
                z = q.match(line.decode('utf-8'))
                if z:
                    self.identifiers.append(z.group('ec'))
            for e in self.identifiers:
                o = self._getNameLocal(e)
                if o == "" and not e.endswith("-"):
                    o = self._getName(e)
                self.results.append((e, o))

    def _getName(self, ec):
        o = ""
        try:
            f = urllib.request.urlopen("https://enzyme.expasy.org/EC/" + ec + ".txt")
            q = re.compile(r'^DE.{3}(?P<name>.*)$')
            for line in f.readlines():
                z = q.match(line.decode('utf-8'))
                if z:
                    o = z.group('name')
                    break

        except Exception as e:
            print("Problem getting EC name: " + ec + "\n" + str(e))
        return o

    def _getNameLocal(self, ec):
        if ec.endswith("-"):
            return ""
        o = ""
        try:
            with open(self.sec_file, encoding='utf-8') as f:
                q = re.compile(r'^ID.{3}(' + ec + ')$')
                id_fnd = False
                for line in f.readlines():
                    z = q.match(line)
                    if z and not id_fnd:
                        q = re.compile(r'^DE.{3}(?P<name>.*)$')
                        id_fnd = True
                        continue
                    if z and id_fnd:
                        o = z.group('name')
                        break

        except Exception as e:
            print("Problem getting EC name: " + ec + "\n" + str(e))
        return o


class MyValidationError(Exception):
    def __init__(self, message, *args):
        self.message = message
        super(MyValidationError, self).__init__(message, *args)


class PfamAnnotation(AnnotationBase):
    '''Retrieve PFAM annotations for a given PDBID and Chain

       Annotations are parsed from files downloaded from pfam and stored locally in:
       $HOME/.local/spyprot/

       Parameters
       ==========
       pdbcode: string
       chain: string

       Return
       ======
       results: list of tuples containing attributes: pdbid, chain, pfam_short, pfam_desc, pfam accession code, knot_type, coverage
       '''

    def __init__(self, pdb, chain='A'):
        super().__init__([PDB_CHAIN_PFAM, PFAM_DESC], pdb, chain)
        self.uniprot = ""
        self.pfam_domains = {}
        self.entangled_regions = {}
        self.coverage = {}
        self.domains = {}
        self.pfam_mappings = {}

        with urllib.request.urlopen(Uni + pdb) as url:
            data = json.loads(url.read().decode())
            self.uniprot = (list(data[pdb]['UniProt'].keys())[0])

        with urllib.request.urlopen(PDBe + self.uniprot) as url:
            data = json.loads(url.read().decode())
            if data[self.uniprot]['Pfam'].keys() != []:
                for i in data[self.uniprot]['Pfam'].keys():
                    self.pfam_domains[i] = data[self.uniprot]['Pfam'][
                        i]  # pfam_domain --> {pfam: {mappinngs, identifier, description}}
                    self.pfam_mappings[i] = [self.pfam_domains[i]['mappings'][0]['unp_start'],
                                             self.pfam_domains[i]['mappings'][0]['unp_end']]

                self.entangled_regions = PfamAnnotation.get_entangled_region(pdb, chain)
                self.coverage, self.domains = PfamAnnotation.get_entangled_domains(self.entangled_regions,
                                                                                   self.pfam_mappings)

                for e in set(self.domains.values()):
                    a, b, c = self._get_names(e)
                    k = self.which_knot(self.domains, e)
                    for i in k:
                        self.results.append((pdb, chain, self.pfam_domains[e]['name'],
                                             self.pfam_domains[e]['description'], e, i, self.coverage[i][e]))

    @staticmethod
    def which_knot(ent_domains, search):
        kn = []
        for knot, domain in ent_domains.items():
            if domain == search:
                kn.append(knot)
        return kn

    def _get_names(self, pfam):
        q = re.compile(r'^' + pfam + '\t(?P<pfam_short>.*)\t(?P<pfam_desc>.*)$')
        with open(self.sec_file, encoding='utf-8') as sd:
            for line in sd:
                z = q.match(line)
                if z:
                    return z.group('pfam_short'), z.group('pfam_desc'), pfam
        return "", "", ""

    @staticmethod
    def get_region(seq1, seq2):
        matrix = substitution_matrices.load("BLOSUM62")
        try:
            alignments = pairwise2.align.localds(seq1, seq2, matrix, -10,
                                                 -0.5)  # czemu akurat te dane --> czy zmiana mialaby wplyw na alignment i przypisaną domene?
            if alignments:
                region = alignments[0][1]
                protein = alignments[0][0]
                i = 0
                j = 0
                if region[0] == "-":
                    while region[i] == "-":
                        i += 1
                    region = region[i:]
                if region[-1] == "-":
                    while region[::-1][j] == "-":
                        j += 1
                    j = j * (-1)
                    region = region[:j]
                return region, i + 1, len(protein) - j * (-1)
        except Exception as e:
            raise MyValidationError("error: problem with the alignment.\n" + str(e))

    @staticmethod
    def fasta2string(fasta):  # fasta[1] takes the first sequence (there can be more than one, eg. structure 2fg7 C)
        out = "".join(fasta[1].split('\n')[1:])
        return out

    @staticmethod
    def nc_cuts_from_knotprot(pdb_id, chain_id):  # find N and C cut of the entanglement region
        candidate_files = ["_nr_chains_knotted_N_C", "_red_chains_knotted_N_C",
                           "_all_chains_knotted_N_C"]  # ??czym sa elementy
        entanglement = {}
        for file in candidate_files:
            with Session() as s:
                site = s.get("https://knotprot.cent.uw.edu.pl/{}".format(file),
                             allow_redirects=False)  # po wydrukowaniu site --> lista struktur i info
                for i in site.text.split('\n'):
                    if pdb_id in i and chain_id in i:  # PDB id;Chain id;Chain Length;{Y,A,T}={Published,Artifact,Not published};{K,S}={Knot,Slipknot};Main knot type (e.g. 31=3.1);N_cut;C_cut;Knots range
                        entanglement_type = "{}{}".format(i.split(';')[4], i.split(';')[5])
                        n_cut = int(i.split(';')[6])
                        c_cut = int(i.split(';')[7])
                        entanglement[entanglement_type] = [n_cut, c_cut]  # e.g. {'K31': [85, 89]}
        if entanglement:
            return entanglement
        else:
            raise MyValidationError("error: PDB id not present in KnotProt API files.")

    @staticmethod
    def seq_structure_from_knotprot(pdb_id,
                                    chain_id):  # returns the sequence of the structure that was used to compute entanglement in KnotProt; it includes only resolved amino acids! no gaps!
        with Session() as s:
            site = s.get("https://knotprot.cent.uw.edu.pl/sequence/{}/{}/".format(pdb_id, chain_id),
                         allow_redirects=False)
            if site.text:
                return site.text
            else:
                raise MyValidationError("error: Couldn't find structure sequence in KnotProt.")

    @staticmethod
    def seq_from_uniprot(
            pdb_id):  # returns fasta of the whole protein, not only the sequence from the structure, sometimes there is more than one sequence
        with Session() as s:
            site = s.get("https://www.uniprot.org/uniprot/?query=database:(type:pdb {})&format=fasta".format(pdb_id),
                         allow_redirects=False)
            if site.text:
                return site.text.split(">")
            else:
                raise MyValidationError("error: Couldn't find protein sequence in UniProt based on PDB id.")

    @staticmethod
    def get_entangled_region(pdb_id, chain_id):
        protein_id = PfamAnnotation.seq_from_uniprot(pdb_id)
        seq_protein = PfamAnnotation.fasta2string(protein_id)  # fasta file
        seq_structure = PfamAnnotation.seq_structure_from_knotprot(pdb_id,
                                                                   chain_id)  # returns the sequence of the structure that was used to compute entanglement in KnotProt
        knotprot_ends = PfamAnnotation.nc_cuts_from_knotprot(pdb_id, chain_id)
        entangled_regions = {}
        for knot in knotprot_ends:
            knotprot_N_end = knotprot_ends[knot][0]
            knotprot_C_end = knotprot_ends[knot][1]
            seq_structure_with_gaps, start_region, stop_region = PfamAnnotation.get_region(seq_protein, seq_structure)
            if seq_structure_with_gaps:
                entangled_region_start = start_region + knotprot_N_end
                entangled_region_end = stop_region - knotprot_C_end
                if not entangled_region_start or not entangled_region_end:
                    entangled_regions[knot] = ["error", "error"]
                else:
                    entangled_regions[knot] = [entangled_region_start, entangled_region_end]
        if not entangled_regions:
            raise MyValidationError("errror: Problem with the localization of the entangled region.")
        return entangled_regions

    @staticmethod
    def get_entangled_domains(dict_knot_ranges,
                              dict_domain_ranges):  # dict_domain_ranges={'PF01746': [22, 221]}; dict_knot_ranges={'K31': [89, 132]}
        coverage = {}
        domains = {}
        for knot in dict_knot_ranges:
            knot_start = dict_knot_ranges[knot][0]
            knot_end = dict_knot_ranges[knot][1]
            coverage[knot] = {}
            domains_joined = ''
            if knot_start != "error":
                knot_length = float(knot_end) - float(knot_start) + 1
                for domain in dict(sorted(dict_domain_ranges.items(), key=lambda item: item[1])):
                    domain_start = dict_domain_ranges[domain][0]
                    domain_end = dict_domain_ranges[domain][1]
                    if domain_start >= knot_start and domain_start <= knot_end and domain_end >= knot_end:
                        perc = (knot_end - domain_start + 1) / knot_length * 100
                        coverage[knot][domain] = round(perc, 1)
                        domains_joined += "-" + domain.split("_")[0]
                    elif domain_start <= knot_start and domain_end >= knot_end:
                        perc = 100
                        coverage[knot][domain] = round(perc, 1)
                        domains_joined += "-" + domain.split("_")[0]
                    elif domain_start <= knot_start and domain_end >= knot_start:
                        perc = (domain_end - knot_start + 1) / knot_length * 100
                        coverage[knot][domain] = round(perc, 1)
                        domains_joined += "-" + domain.split("_")[0]
                    elif domain_start >= knot_start and domain_end <= knot_end:
                        perc = (domain_end - domain_start + 1) / knot_length * 100
                        coverage[knot][domain] = round(perc, 1)
                        domains_joined += "-" + domain.split("_")[0]
            domains[knot] = domains_joined[1:]
        return coverage, domains
