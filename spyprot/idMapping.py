import requests
import json
from os.path import isfile
from tqdm import tqdm  # completely optional; we might exclude it


############ Functions used to map IDs from one DB to another
##############################################################
# USAGE:
# generate_mappings(just_pdb_to_uniprot=False) | (dumps to file; checks version and updates automatically; ignores if current)
# PDB_Uniprot("P10970")                        | (checks mappings in the dump files, accepts PDB ID, PDB ID + " " + Chain, or Uniprot ID; no need to specify)


# Help functions
# Download file with replacament if already same filename locally
def download_file(url, output_file, headers=""):
    """Downloads file from the web. Doesn't check for overwrites."""
    with open(output_file, 'w') as new_local_file:
        new_local_file.write(requests.get(url, headers=headers).content.decode('utf-8'))


def json_dump(input_object, output_file):
    """Dumps the input Python object to a .json file using the JSON module. Doesn't check for overwrites."""
    try:
        with open(output_file, 'w') as outf:
            json.dump(input_object, outf)
        return True
    except:
        print("Something went wrong when dumping file")
        return False


def json_load(input_file):
    """Loads a previously dumped Python object from a .json file using the JSON module."""
    try:
        with open(input_file, 'r') as inpf:
            loaded = json.load(inpf)
        return loaded
    except:
        print("Couldn't load.")
        return False


##### PDB <=> Uniprot
# Check if local version of the mapping file is most recent and replace if not
# The returned value says whether the local file is up to date with the remote file
# None means: we don't know because we couldn't fint the online file.
# If no local file we return Fale.
# If you just want to check, best to set no_download to True
def PDB_Uniprot_update_list(local_list_file_path="pdb_chain_uniprot.lst", no_download=False):
    """Two-way mapping between PDB ID and Uniprot IDs.
        1. It checks whether the local version of the mapping file is the most recent and replaces it otherwise.
        2. The returned value (T/F) reports whether the local file is up to date with the remote file.
       If it returns None it means that the online file could not be located. If no local file found, return False.
       Used with the no_download=True option can serve to check version of the file without updating."""
    try:
        # Constant URL for periodically updated PDB<=>Uniprot mappings
        remote_url = "http://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst"
        # Load first 100 bytes of the online file to extract its header
        remote_url_headers = {"Range": "bytes=0-100"}
        # Extracting the header (first line)
        remote_file_content = requests.get(remote_url, headers=remote_url_headers).content.decode('utf-8')
        remote_file_header = remote_file_content.split("\n")[0]
    except:
        print("Cannot establish connection. No remote file under this address.")
        return None

    # Check if local file exists
    try:
        # If it does, extract local file's header
        with open(local_list_file_path, 'r') as local:
            local_file_header = local.readline().strip()

        # Compare the local file's header to the remote file's header
        if local_file_header == remote_file_header:
            print("The local mapping file is up to date.")
            return True
        else:
            print("There is a more up to date mapping file available online. Replacing it.")
            if not no_download:
                download_file(remote_url, local_list_file_path)
            return False

    except:
        # No file found. Downloading from remote.
        print("No local file found.")
        if not no_download:
            print("Downloading the remote file")
            download_file(remote_url, local_list_file_path)
        return False


# Parse mapping and put it into a dictionary with "PDBID Chain" key/value pattern
# PDB => Uniprot by default; set false to have it the other way round
def parse_mapping(local_list_file_path="pdb_chain_uniprot.lst", pdb_to_uniprot=True):
    try:
        with open(local_list_file_path, 'r') as raw_mapping:
            # Mapping as a list of lists, each line already divided into columns
            # Already without a header
            raw_mapping = [x.strip().split("\t") for x in raw_mapping.readlines()][1:]
    except:
        print("No local mapping file found. Download it first.")
        return False

    # Create a dictionary with the mapping, depending on the direction
    mapping = dict()
    if pdb_to_uniprot:
        print("Parsing the PDB => Uniprot mapping")
        for line in tqdm(raw_mapping):
            mapping[line[0] + " " + line[1]] = line[2]
    else:
        print("Parsing the Uniprot => PDB mapping")
        for line in tqdm(raw_mapping):
            if line[2] in mapping.keys():
                mapping[line[2]].append(line[0] + " " + line[1])
            else:
                mapping[line[2]] = [line[0] + " " + line[1]]

    return mapping


# Generates mappings if a more recent version can be found and saves them to file
# Returns True if generated anything, returns False if nothing has been generated
def generate_mappings(local_file="pdb_chain_uniprot.lst", just_pdb_to_uniprot=True,
                      pdb_to_uni_mapping_path="pdb_to_uni_mapping.json",
                      uni_to_pdb_mapping_path="uni_to_pdb_mapping.json"):
    """Generates .json dump files containing Python dictionaries with PDB ID=>Uniprot ID and Uniprot ID=>PDB ID mappings."""
    # Check if up to date without downloading first
    if not PDB_Uniprot_update_list(local_list_file_path=local_file, no_download=True):
        # If not, download
        PDB_Uniprot_update_list(local_list_file_path=local_file)

        # Generate new mappings
        json_dump(parse_mapping(pdb_to_uniprot=True), pdb_to_uni_mapping_path)

        if not just_pdb_to_uniprot:
            json_dump(parse_mapping(pdb_to_uniprot=False), uni_to_pdb_mapping_path)

        return True

    # If mapping source file up to date but one of the generated mapping .json files is missing, generate new mapping
    elif not isfile(pdb_to_uni_mapping_path) or not isfile(uni_to_pdb_mapping_path):
        if not isfile(pdb_to_uni_mapping_path):
            json_dump(parse_mapping(pdb_to_uniprot=True), pdb_to_uni_mapping_path)

        if not isfile(uni_to_pdb_mapping_path) and not just_pdb_to_uniprot:
            json_dump(parse_mapping(pdb_to_uniprot=False), uni_to_pdb_mapping_path)

        return True

    else:
        print("Everything is up to date; no need to generate new mapping files.")
        return False


# Look up mapping PDB <=> Uniprot
def PDB_Uniprot(search_key, pdb_to_uniprot_mapping_path="pdb_to_uni_mapping.json",
                uniprot_to_pdb_mapping_path="uni_to_pdb_mapping.json"):
    """Performs lookup in the .json-dumped dictionaries with mappings, both PDB => Uniprot and Uniprot => PDB.
    Search_key can be either a PDBid, a "PDBid Chain" pair or a Uniprot ID; no need to signal it to the function, it recognizes it automatically.
    If search_key == PDBid: return the list of all available mappings
    If search_key == PDBid + " " + Chain: return the only available mapping
    If search_key == UniprotID: return the list of all mappings"""
    if isfile(pdb_to_uniprot_mapping_path):
        try:
            # load PDB => Uniprot mapping
            pdb_to_uni = json_load(pdb_to_uniprot_mapping_path)
        except:
            print("Something went wrong loading the PDB => Uniprot .json mapping file.")
            return False

    if isfile(uniprot_to_pdb_mapping_path):
        try:
            # load Uniprot => PDB mapping
            uni_to_pdb = json_load(uniprot_to_pdb_mapping_path)
        except:
            print("Something went wrong loading the Uniprot => PDB .json mapping file.")

    # If none could be loaded, abort
    if 'pdb_to_uni' not in locals() and 'uni_to_pdb' not in locals():
        print("No mapping could be loaded. Aborting lookup.")
        return False

    # If search_key == PDBid: return the list of all available mappings
    if len(search_key.split()) == 1 and len(search_key) == 4:
        if 'pdb_to_uni' in locals():
            matching_chains = [x for x in pdb_to_uni.keys() if x.startswith(search_key)]
            if len(matching_chains) > 0:
                matches = [tuple([matching_key, pdb_to_uni[matching_key]]) for matching_key in matching_chains]
                return matches
            else:
                print("No mapping found for PDB ID:", search_key)
                return False
        else:
            print("Missing PDB => Uniprot mapping.")
            return False

    # If search_key == PDBid + " " + Chain: return the only available mapping
    elif len(search_key.split()) == 2 and len(search_key.split()[0]) == 4:
        if 'pdb_to_uni' in locals():
            if search_key in pdb_to_uni.keys():
                return pdb_to_uni[search_key]
            else:
                print("No mapping found for:", search_key)
                return False
        else:
            print("Missing PDB => Uniprot mapping.")
            return False

    # If search_key == UniprotID: return the list of all mappings
    elif len(search_key.split()) == 1 and len(search_key) > 4:
        if 'uni_to_pdb' in locals():
            if search_key in uni_to_pdb.keys():
                return uni_to_pdb[search_key]
            else:
                print("No mapping found for Uniprot ID:", search_key)
                return False
        else:
            print("Missing Uniprot => PDB mapping.")
            return False

    # Incorrect search_key
    else:
        print(
            "Incorrect search phrase. Are you sure you're trying to map PDB <=> Uniprot?\nLegal search phrases: 1s72, 1s72 A, P30615")
        return False
