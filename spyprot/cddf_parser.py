#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib.error
import urllib.request
from os import path
from re import compile
from spyprot import Command
import os

DEFAULT_DATA_FILE_PATH = path.join(path.expanduser("~"), ".local", "spyprot")
hcath_tmp = DEFAULT_DATA_FILE_PATH + '/cath-domain-description-file.txt_tmp'
hcath = DEFAULT_DATA_FILE_PATH + '/cath-domain-description-file.txt'
cath_uri = \
    'http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-description-file.txt'


# egrep -v "^NAME|^DLEN|^S|^VER|^NSE|^FOR|^DSEQ" CathDomainDescriptionFile_org >CathDomainDescriptionFile


def getCath(pdbcode, chain):
    '''Retrieve CATH annotations for a given PDBID and Chain
       Supported attributes:
       CLASS
       CATHCODE
       DOMAIN
       ARCH
       TOPOL
       HOMOL

       Annotations are parsed from a cath-domain-description file downloaded from cathdb and stored locally in:
       $HOME/.local/spyprot/

       Parameters
       ==========
       pdbcode: string
       chain: string

       Return
       ======
       results: list of CATH annotations
       '''
    if not path.isfile(hcath):
        if not path.exists(DEFAULT_DATA_FILE_PATH):
            os.makedirs(DEFAULT_DATA_FILE_PATH)
        urllib.request.urlretrieve(cath_uri, hcath_tmp)
        cmd = Command('egrep -v "^NAME|^DLEN|^S|^VER|^NSE|^FOR|^DSEQ" ' + hcath_tmp + '>' + hcath)
        cmd.run(1000)
        os.remove(hcath_tmp)
    pdbcode = pdbcode.lower()
    chain = chain
    record = compile(r'^DOMAIN    ' + pdbcode + chain)
    results = []
    with open(hcath, 'r', encoding='utf-8') as cddf:
        found = False
        for line in cddf:
            if record.match(line):
                found = True
                recorddd = {'pdbcode': pdbcode, 'chain': chain}
            elif 'ENDSEG' in line:
                if found:
                    results.append(recorddd)
                found = False
            if found:
                key = line[:10].strip()
                if key not in recorddd:  # only one line..
                    recorddd[key] = line[10:].strip()
    return results


if __name__ == "__main__":
    d = getCath("1j85", "A")
    print(d)
