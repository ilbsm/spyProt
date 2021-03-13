#!/usr/bin/env python3

#                                                                          |--------------|
#    .d8888b.                       8888888b.                  888         | v. alpha 0.3 |
#   d88P  Y88b                      888   Y88b                 888         |  19 Dec 2018 |
#   Y88b.                           888    888                 888         |--------------|
#    "Y888b.   88888b.  888  888    888   d88P 888d888 .d88b.  888888 
#      "Y88b. 888 "88b 888  888     8888888P"  888P"  d88""88b 888    
#         "888 888  888 888  888    888        888    888  888 888    
#   Y88b  d88P 888 d88P Y88b 888    888        888    Y88..88P Y88b.  
#    "Y8888P"  88888P"   "Y88888    888        888     "Y88P"   "Y888 
#              888           888                                   
#              888      Y8b d88P  |------------------------------------------                              
#              888       "Y88P"   | Rapid protein overview from multiple DBs |                              
#   ------------------------------------------------------------------------ |
#   |                            developed                                   |
#   |                               by                                       |
#   |                        Borys Jastrzębski                               |
#   ------------------------------------------------------------------------ |

import wget
from sys import argv
from datetime import datetime
import os
import re


# from prettytable import PrettyTable

###########################################################################################
#                                                                                         #
# To Think/Do: timestamps on result files?                                                #
#              maybe a graphical interface?                                               #
#              function for cleaning up old fastas?                                       #
#              option not to timestamp fastas?                                            #
#              reformat into a table                                                      #
#              collect clan info from the PFAM page                                       #
#              Solve dependency problems on Ubuntu (does the depInstall.sh work?)         #
#              Add option to accessPDB to include only the chain A in the fasta file      #
#              Add knot type distinction (more info needed from KnotProt)                 #
#              Check for multiple chains in PFAM and KnotProt (and also PDB if you don't) #
###########################################################################################

# Input to all: PDB ID

def accessPDB(protID, downloadFasta=True, fastaTimestamp=True):
    # Downloading a temp file to read from
    tempStamp = datetime.now().isoformat()
    try:
        wget.download("https://files.rcsb.org/header/" + protID + ".pdb", ".temps/" + protID + "_" + tempStamp + ".txt")
        print("\nPDB | Database succesfully accessed at ID:", protID)
    except:
        print("PDB | Can't access the database at ID:", protID)
        return None

    # Downloading the fasta sequence
    if downloadFasta == True:
        # Adding timestamp to the output file based on the preference
        if fastaTimestamp == True:
            fastaStamp = datetime.now().isoformat()
        else:
            fastaStamp = ""

        # Don't forget to add chain recognition and choice
        # Actually both are stored in one file
        # Ask Lina what to do about it
        try:
            wget.download(
                "https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=" + protID + "&compressionType=uncompressed",
                "results/fasta/" + protID + "_" + fastaStamp + ".fasta")
            print("\nPDB | The fasta file succesfully downloaded for ID:", protID)
        except:
            print("PDB | Couldn't download the fasta file for ID:", protID)

    # Processing header to extract the remaining data
    with open(".temps/" + protID + "_" + tempStamp + ".txt") as temp:
        # All data to a list of lines
        data = temp.readlines()

    srcOrganism = list(
        (line.split("ORGANISM_SCIENTIFIC:")[1].strip()[:-1] for line in data if "ORGANISM_SCIENTIFIC:" in line))
    protFunction = list((line.split("COMPND")[1].strip() for line in data if "COMPND" in line))

    # Clean up the temporary file
    os.remove(".temps/" + protID + "_" + tempStamp + ".txt")

    return (protID, srcOrganism, protFunction, fastaStamp)


# At the moment it returns a boolean: knotted (in the db) or unknotted (not present in the db)
# WARNING: At the moment it checks only for chain A of the protein.
def accessKnotProt(protID):
    tempStamp = datetime.now().isoformat()
    try:
        # Temporary file in the making
        wget.download("http://knotprot.cent.uw.edu.pl/browse/?pdbid=" + protID,
                      ".temps/" + protID + "_" + tempStamp + ".txt")
        print("\nKnotProt | Database succesfully accessed at ID:", protID)
    except:
        print("KnotProt | Can't access the database at ID:", protID)
        return None

    # Processing the file
    with open(".temps/" + protID + "_" + tempStamp + ".txt") as temp:
        # All data to a list of lines
        data = temp.read()

    # Check the topmost box to see whether knotted or unknotted
    resultBox = data.split('<div class="alert alert-info alert-dismissable">')[1].split("</div>")[0].replace("/",
                                                                                                             "").split(
        "<strong>")

    # Clean up the temporary file
    os.remove(".temps/" + protID + "_" + tempStamp + ".txt")

    if "Knotted" in resultBox:

        tempStamp = datetime.now().isoformat()
        try:
            # Temporary file in the making
            wget.download("http://knotprot.cent.uw.edu.pl/view/" + protID + "/A/",
                          ".temps/" + protID + "_" + tempStamp + ".txt")
            print("\nKnotProt | View succesfully accessed at ID:", protID)
        except:
            print("KnotProt | Can't access the view at ID:", protID)
            return None

        # Processing the file
        with open(".temps/" + protID + "_" + tempStamp + ".txt") as temp:
            # All data to a list of lines
            typeData = temp.read()

        # Check the tag span to check if it's slipknot or knot
        typeData = typeData.split('<span class="label label-primary">')[1].split("</span>")[0]
        # If knotType == True, it's a Slipknot. Else: it's a Knot.
        knotType = ["Knot", "Slipknot"]["Slipknot" in typeData]

        # Clean up the temporary file
        os.remove(".temps/" + protID + "_" + tempStamp + ".txt")

        print("The protein at PDB ID:", protID, "is knotted. It contains a", knotType)
        return knotType

    elif "Unknotted" in resultBox:
        print("The protein at PDB ID:", protID, "is unknotted.")
        return 0
    else:
        print("KnotProt | Error: Can't process the result box for ID:", protID)
        return None


# Takes info for chain A only, and from the PDB part of the PFAM table
# Rewrite to include other chains and then return them to the accessKnotProt function
# To get info on the knots for each chain
def accessPFAM(protID):
    tempStamp = datetime.now().isoformat()
    try:
        # Temporary file in the making
        wget.download("https://pfam.xfam.org/structure/" + protID + "#tabview=tab3",
                      ".temps/" + protID + "_" + tempStamp + ".txt")
        print("\nPFAM | Database succesfully accessed at ID:", protID)
    except:
        print("PFAM | Can't access the database at ID:", protID)
        return None

    with open(".temps/" + protID + "_" + tempStamp + ".txt", 'r') as temp:
        data = temp.read()

    # Clean up the temporary file
    os.remove(".temps/" + protID + "_" + tempStamp + ".txt")

    # Now let's do a simple string split to get the first <tr> of the tbody of the table out of the code.
    # (Only the first row because we're interested in chain A only)
    trOne = "".join(
        data.split('<table id="structuresTable"')[1].split("</table>")[0].split("<tbody>")[1].split("</tr>")[0].split()[
        2:])

    # Then a regex pattern for splitting to get the individual rows of the table
    rows = re.compile(r'<\/?td*>')
    ## Again, the first element will be the name of the chain so we can skip it as it is always A
    table = [x for x in rows.split(trOne) if x != ""][1:]

    # All the relevant information extracted and final values assigned to variables.
    ## The PDB prefix is here in case of potential extension to include seqs starts and ends from UniProt.
    PDBseqStart, PDBseqEnd, famID, famName = table[0], table[1], table[-1].split("family/")[1].split('"')[0], \
                                             table[-1].split(">")[1].split("<")[0]

    # Now we can get the clan info for the family, if it exists
    tempStamp = datetime.now().isoformat()
    try:
        # Temporary file in the making
        wget.download("https://pfam.xfam.org/family/" + famID + "#tabview=tab2",
                      ".temps/" + famID + "_" + tempStamp + ".txt")
        print("\nPFAM | Database succesfully accessed at Family ID:", famID)
    except:
        print("PFAM | Can't access the clan database at Family ID:", famID)
        return None

    with open(".temps/" + famID + "_" + tempStamp + ".txt", 'r') as temp:
        clanData = temp.read()

    # Clean up the temporary file
    os.remove(".temps/" + famID + "_" + tempStamp + ".txt")

    try:
        # Take first <p> of the blockContent div of the clanBlock div
        clanData = \
        clanData.split('<div class="block" id="clanBlock">')[1].split('<div class="blockContent">')[1].split("<p>")[
            1].split("</p>")[0]

        # Simple split to get the name and CL ID from between the <a> handles
        clanID = "CL" + clanData.split("/clan/CL")[1].split('"')[0]
    except:
        print("PFAM | No clan information found for protein family", famID)
        clanID = None

    return (clanID, famID, famName, PDBseqStart, PDBseqEnd)


#########################################################################################################
# SOLVED – THIS IS AN UNNECESSARY FUNCTION AT THE MOMENT. (after you change the display function)       #
# accessPFAM works properly now
# Could serve as a backup                                                                         #
#########################################################################################################
# The mapping file between PDB ID and PFAM comes from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/. #
# It is a file from 2015 – a bit old, I know, but it's just a temporary measure for the purpose         #
# of this proof of concept.                                                                             #
#########################################################################################################
def mapToPFAM(protID):
    with open(".resources/pdb_pfam_mapping.txt") as mapping:
        data = mapping.readlines()

    header = data[0]
    familyInfo = list((line for line in data if protID.upper() in line))

    if familyInfo == []:
        familyInfo = "No record"

    return (header, familyInfo)


#########################################################################################################
#########################################################################################################

# Deprecated / will be rewritten from scratch
def combineResults(fromPDB, fromKP, fromPFAM):
    # Write out results to a result file
    if fromPDB != None:

        # Combine everything
        # Write to a result file
        # fromPDB[0] = protID
        with open("results/result_" + fromPDB[0] + ".txt", 'w') as overallOutput:

            # New line in the result file
            def newLine():
                overallOutput.write("\n")

            overallOutput.write("PDB ID of the protein: " + fromPDB[0])
            overallOutput.write("\n--------------------------------------\n\n")

            overallOutput.write("[PDB DATA]")
            newLine()
            overallOutput.write("The protein comes from: " + ", ".join(fromPDB[1]))
            newLine()
            overallOutput.write("Short overview: " + ", ".join(fromPDB[2]))
            overallOutput.write("\n--------------------------------------\n\n")

            overallOutput.write("[KnotProt DATA]")
            newLine()
            overallOutput.write("The protein is: " + ["UNKNOTTED", "KNOTTED"][fromKP])
            overallOutput.write("\n--------------------------------------\n\n")

            overallOutput.write("[PFAM DATA]")
            newLine()
            overallOutput.write(fromPFAM[0])
            overallOutput.write("\n".join(fromPFAM[1]))
            overallOutput.write("\n--------------------------------------\n\n")

            try:
                with open("results/fasta/" + fromPDB[0] + "_" + fromPDB[3] + ".fasta", 'r') as fl:
                    fasta = fl.read()

                overallOutput.write(
                    "[FASTA Sequence from PDB] (results/fasta/" + fromPDB[0] + "_" + fromPDB[3] + ".fasta)")
                newLine()
                overallOutput.write(fasta)
                overallOutput.write("\n--------------------------------------\n\n")
            except:
                print("Trouble finding the fasta file when printing the results. Perhaps nothing was found.")

        print("Result file generated properly. Check the results folder.")
        return "results/result_" + fromPDB[0] + ".txt"

    else:
        print("No protein ID returned by the accessPDB function.")
        return False


def display(filename):
    # display the results file
    try:
        # For MacOS
        os.system("open " + filename)
    except:
        try:
            # For Windows
            os.system("start " + filename)
        except:
            print("Something went wrong displaying the file.")


# Function to clean up older fasta files for one protein
def cleanFastas():
    pass


# Safety measure
def cleanTemps():
    # Clean up all the remaining temporary files
    for temp in [f for f in os.listdir(".temps")]:
        os.remove(".temps/" + temp)


if __name__ == "__main__":
    # Creating necessary folders
    if not os.path.exists(".temps"):
        os.makedirs(".temps")

    if not os.path.exists("../results"):
        os.makedirs("results")

    if not os.path.exists("../results/fasta"):
        os.makedirs("results/fasta")

    # Main loop
    # Since not all of the main functions are ready, this is left commented out
    #
    # for proteinID in argv[1:]:
    #    results = [accessPDB(proteinID), accessKnotProt(proteinID), accessPFAM(proteinID)]
    #    #display(combineResults(*results))
    #    print(results)

    # Extra clean-up
    cleanTemps()
