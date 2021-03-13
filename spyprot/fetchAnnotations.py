#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Pawel Rubach 2019
# based on:
# Copyright Michal Jamroz, 2014, jamroz@chem.uw.edu.pl
import datetime
import urllib.request, urllib.error, urllib.parse, gzip, re
import os
from os import path

REFRESH_FILE_INTERVAL = 7 * 24 * 3600  # Refresh every week

PDB_CHAIN_ENZYME = ["ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_enzyme.tsv.gz",
                    "pdb_chain_enzyme.tsv.gz"]
ENZYME_DAT = ["https://ftp.expasy.org/databases/enzyme/enzyme.dat", "enzyme.dat"]

PDB_CHAIN_PFAM = ["ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_pfam.tsv.gz",
                  "pdb_chain_pfam.tsv.gz"]
PFAM_DESC = ["http://pfam.xfam.org/families?output=text", "pdb_chain_pfam_desc"]

DEFAULT_DATA_FILE_PATH = path.join(path.expanduser("~"), ".local", "spyprot")


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

    def download_if_not_exist(self, file_locs):
        enz_file = path.join(self.data_file_path, file_locs[1])
        if not path.isfile(enz_file) or (not path.isfile(enz_file + ".flg") and path.isfile(
                enz_file) and datetime.datetime.now().timestamp() - os.stat(
                enz_file).st_mtime > self.refresh_file_interval):
            print('downloading file: ' + enz_file + " from: " + file_locs[0])
            AnnotationBase.touch(enz_file + ".flg")
            try:
                f = urllib.request.urlopen(file_locs[0])
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
       results: list of tuples containing attributes: pdbid, chain, pfam_short, pfam_desc and pfam accession code
       '''

    def __init__(self, pdb, chain='A'):
        super().__init__([PDB_CHAIN_PFAM, PFAM_DESC], pdb, chain)
        q = re.compile(r'^' + self.pdb + '\t' + self.chain + '\t(.*)\t(?P<pfam>.*)\t(.*)$', re.M)
        with gzip.open(self.first_file) as s:
            for line in s.readlines():
                line = line.decode('utf-8')
                z = q.match(line)
                if z:
                    self.identifiers.append(z.group('pfam'))
            for e in self.identifiers:
                a, b, c = self._getNames(e)
                self.results.append((pdb, chain, a, b, c))

    def _getNames(self, pfam):
        q = re.compile(r'^' + pfam + '\t(?P<pfam_short>.*)\t(?P<pfam_desc>.*)$')
        with open(self.sec_file, encoding='utf-8') as sd:
            for line in sd:
                z = q.match(line)
                if z:
                    return z.group('pfam_short'), z.group('pfam_desc'), pfam
        return "", "", ""

# class getAnnotations:
#     '''Retrieve PFAM annotations for a given PDBID and Chain directly from RCSB
#
#        Parameters
#        ==========
#        pdbcode: string
#        chain: string
#
#        Return
#        ======
#        results: list of tuples containing attributes: pdbid, chain, ..., ..., pfam accession code, pfam_short, pfam_desc, ...
#        '''
#     def __init__(self, pdb, chain='A'):
#         self.pdb = pdb
#         self.chain = chain
#         url = "http://www.pdb.org/pdb/rest/hmmer?structureId="+pdb.upper()
#         r = urllib.request.urlopen(url)
#         data = r.read()
#         self.hmmer = etree.fromstring(data)
#
#     def get(self):
#         out=[]
#         for r in self.hmmer:
#             d = list(r.values())
#             d[0] = d[0].lower()
#             out.append(d)
#         #out[0] = out[0].lower()
#         return out
