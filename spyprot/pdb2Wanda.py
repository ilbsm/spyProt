import urllib

import numpy as np

#######
# Author: Pawel Dabrowski-Tumanski
# p.dabrowski@cent.uw.edu.pl
# Date 02.10.2015; version from 28.01.2016
# Script reads the pdb file and gets the chain information. Each chain forms individual class. In each class the coordinates of CA atoms, missing atoms, residues, SS-bonds and links are stored.
# Script converts pdb file to xyz file (one for each chain), in 5-column format (with residue as the last column), fills the gap with straight line and prints warning if the gap is longer than 6 residues.
# Finally it gives the commands for Wanda's program for finding the crossing of the surfaces.
#######

################################ Chain Class definition and possible amino acids ################################
amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CIT', 'CYS', 'GLU', 'GLN', 'GLY', 'HEY', 'HIS', 'HYP', 'ILE', 'LEU', 'LYS',
               'MET', 'ORN', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'SEC', 'PYL', 'ASX', 'GLX', 'XLE', 'XAA',
               'MSE', 'FGL', 'LLP', 'SAC', 'PCA', 'MEN', 'CSB', 'HTR', 'PTR', 'SCE', 'M3L', 'OCS', 'KCX', 'SEB', 'MLY',
               'CSW', 'TPO', 'SEP', 'AYA', 'TRN', 'DAL', 'DAR', 'DSG', 'DAS', 'DCY', 'DGN', 'DGL', 'DHI', 'DIL', 'DLE',
               'DLY', 'MED', 'DPN', 'DPR', 'DSN', 'DTH', 'DTR', 'DTY', 'DVA']


def find_index(number, arr):
    for k in range(len(arr)):
        if (arr[k][0] == number):
            return k


class Chain:
    def __init__(self, name):
        self.name = name
        self.residues = []
        self.coordinates = []
        self.bridges = []
        self.missing = []
        self.gaps = []
        self.helix = []
        self.sheet = []

    def find_length(self):
        self.length = max(len(self.residues), len(self.coordinates))

    def add_residue(self, line):
        line = line[19:70]
        for k in range(len(line.split())):
            if (line.split()[k] in amino_acids):
                self.residues.append(line.split()[k])
        self.find_length()

    def add_bridge(self, bridge):
        self.bridges.append(bridge)

    def add_helix(self, helix):
        self.helix.append(helix)

    def add_sheet(self, sheet):
        self.sheet.append(sheet)

    def add_missing(self, missing, residue):
        self.missing.append([missing, residue])

    def add_coordinate(self, index, coordinate, residue):
        if (len(self.coordinates) == 0):
            for k in range(len(self.missing)):
                if (index > self.missing[k][0]):
                    self.coordinates.append([self.missing[k][0], [], self.missing[k][1]])
            self.coordinates.append([index, coordinate, residue])
        else:
            diff = index - self.coordinates[len(self.coordinates) - 1][0]
            if (diff > 1):
                x = self.coordinates[len(self.coordinates) - 1][1][0]
                y = self.coordinates[len(self.coordinates) - 1][1][1]
                z = self.coordinates[len(self.coordinates) - 1][1][2]
                vec_diff = [(coordinate[0] - x) / diff, (coordinate[1] - y) / diff, (coordinate[2] - z) / diff]
                self.gaps.append([diff, self.coordinates[len(self.coordinates) - 1][0], index])
            for k in range(diff - 1):
                x = round(x + vec_diff[0], 3)
                y = round(y + vec_diff[1], 3)
                z = round(z + vec_diff[2], 3)
                #    res=self.missing[find_index(index-diff+k+1,self.missing)][1]
                #    if (res==None):
                res = "XXX"
                self.coordinates.append([index - diff + k + 1, [x, y, z], res])
            self.coordinates.append([index, coordinate, residue])

    ######### characterize the bond type
    def bond_type(self, res1, atom1, Nend, res2, atom2, Cend):
        if ((atom1[0], atom2[0]) == ("C", "N")):
            if ((((res1, atom1) == ("GLU", "CD")) or ((res1, atom1) == ("ASP", "CG"))) and (
                    (res2, atom2) == ("LYS", "NZ"))):
                return "Amide"
            else:
                # return "Amide-like"
                return "Amide"
        if ((atom1[0], atom2[0]) == ("N", "C")):
            if ((((res1, atom1) == ("LYS", "NZ")) or (Nend and atom1 == "N")) and (
                    ((res2, atom2) == ("GLU", "CD")) or ((res2, atom2) == ("ASP", "CG")) or (Cend and atom2 == "C"))):
                return "Amide"
            else:
                # return "Amide-like"
                return "Amide"
        if ((atom1[0], atom2[0]) == ("C", "O")):
            if ((((res1, atom1) == ("GLU", "CD")) or ((res1, atom1) == ("ASP", "CG"))) and (
                    ((res2, atom2) == ("SER", "OG")) or ((res2, atom2) == ("THR", "OG1")))):
                return "Ester"
            else:
                # return "Ester-like"
                return "Ester"
        if ((atom1[0], atom2[0]) == ("O", "C")):
            if ((((res1, atom1) == ("SER", "OG")) or ((res2, atom2) == ("THR", "OG1"))) and (
                    ((res2, atom2) == ("GLU", "CD")) or ((res2, atom2) == ("ASP", "CG")) or (Cend and atom2 == "C"))):
                return "Ester"
            else:
                # return "Ester-like"
                return "Ester"
        if ((atom1[0], atom2[0]) == ("C", "S")):
            if ((((res1, atom1) == ("GLU", "CD")) or ((res1, atom1) == ("ASP", "CG"))) and (
                    (res2, atom2) == ("CYS", "SG"))):
                return "Thioester"
            else:
                # return "Thioester-like"
                return "Thioester"
        if ((atom1[0], atom2[0]) == ("S", "C")):
            if (((res1, atom1) == ("CYS", "SG")) and (
                    ((res2, atom2) == ("GLU", "CD")) or ((res2, atom2) == ("ASP", "CG")) or (Cend and atom2 == "C"))):
                return "Thioester"
            else:
                # return "Thioester-like"
                return "Thioester"
        else:
            if atom1 == 'SG' and atom2 == 'SG':
                return 'SS'
            else:
                return 'Others'

    ######## check, whether it is N- or C-end
    def N_end(self, index):
        if (index == self.coordinates[0][0]):
            return True
        else:
            return False

    def C_end(self, index):
        if (index == self.coordinates[len(self.coordinates) - 1][0]):
            return True
        else:
            return False

    ######### getting taxonomy
    def get_taxonomy(self, name):
        url = "http://www.uniprot.org/uniprot/" + name
        source = urllib.urlopen(url)
        for line in source.decode('utf-8'):
            if ("/taxonomy/33208\"" in line):
                self.taxonomy = "Animal"
                break
            if ("/taxonomy/33090\"" in line):
                self.taxonomy = "Plant"
                break
            if ("/taxonomy/4751\"" in line):
                self.taxonomy = "Fungus"
                break
            if ("/taxonomy/10239\"" in line):
                self.taxonomy = "Virus"
                break
            if ("/taxonomy/2\"" in line):
                self.taxonomy = "Bacteria"
                break
            if ("/taxonomy/2157\"" in line):
                self.taxonomy = "Archea"
                break
            if ("/taxonomy/33630\"" in line):
                self.taxonomy = "Alveolata"
                break
            if ("/taxonomy/554915\"" in line):
                self.taxonomy = "Amoebozoa"
                break
            if ("/taxonomy/554296\"" in line):
                self.taxonomy = "Apusozoa"
                break
            if ("/taxonomy/1401294\"" in line):
                self.taxonomy = "Breviatea"
                break
            if ("/taxonomy/193537\"" in line):
                self.taxonomy = "Centroheliozoa"
                break
            if ("/taxonomy/3027\"" in line):
                self.taxonomy = "Cryptophyta"
                break
            if ("/taxonomy/33682\"" in line):
                self.taxonomy = "Euglenozoa"
                break
            if ("/taxonomy/207245\"" in line):
                self.taxonomy = "Fornicata"
                break
            if ("/taxonomy/38254\"" in line):
                self.taxonomy = "Glaucocystophyceae"
                break
            if ("/taxonomy/2830\"" in line):
                self.taxonomy = "Haptophyceae"
                break
            if ("/taxonomy/5752\"" in line):
                self.taxonomy = "Heterolobosea"
                break
            if ("/taxonomy/556282\"" in line):
                self.taxonomy = "Jakobida"
                break
            if ("/taxonomy/339960\"" in line):
                self.taxonomy = "Katablepharidophyta"
                break
            if ("/taxonomy/136087\"" in line):
                self.taxonomy = "Malawimonadidae"
                break
            if ("/taxonomy/66288\"" in line):
                self.taxonomy = "Oxymonadida"
                break
            if ("/taxonomy/5719\"" in line):
                self.taxonomy = "Parabasalia"
                break
            if ("/taxonomy/543769\"" in line):
                self.taxonomy = "Rhizaria"
                break
            if ("/taxonomy/2763\"" in line):
                self.taxonomy = "Rhodophyta"
                break
            if ("/taxonomy/33634\"" in line):
                self.taxonomy = "Stramenopiles"
                break
            if ("/taxonomy/33630\"" in line):
                self.taxonomy = "Alveolata"
                break
            else:
                self.taxonomy = "Other"

    ######### cleaning data
    def clean(self):
        ### adding C-end
        if (len(self.coordinates) != 0):
            for k in range(len(self.missing)):
                if (self.coordinates[len(self.coordinates) - 1][0] < self.missing[k]):
                    self.coordinates.append([self.missing[k][0], [], self.missing[k][1]])
        else:
            for k in range(len(self.missing)):
                self.coordinates.append([self.missing[k][0], [], self.missing[k][1]])
                ### clearing double CA atoms
        for k in range(len(self.coordinates) - 1, 0, -1):
            if (self.coordinates[k][0] == self.coordinates[k - 1][0]):
                self.coordinates.pop(k)
        #   if (self.coordinates_residue[k][0]==self.coordinates[k-1][0]):
        #    self.coordinates_residue.pop(k)
        ### clearing non-protein links
        for k in range(len(self.bridges) - 1, -1, -1):
            if ((self.bridges[k][1] not in amino_acids) or (self.bridges[k][4] not in amino_acids)):
                self.bridges.pop(k)
        ### clearing to small loops
        for k in range(len(self.bridges) - 1, -1, -1):
            if (abs(self.bridges[k][3] - self.bridges[k][6]) < 5):
                self.bridges.pop(k)
        ### defining bond type
        for k in range(len(self.bridges)):
            if (self.bridges[k][0] == "LINK"):
                self.bridges[k][0] = self.bond_type(self.bridges[k][1], self.bridges[k][2],
                                                    self.N_end(self.bridges[k][3]), self.bridges[k][4],
                                                    self.bridges[k][5], self.C_end(self.bridges[k][6]))
        ### removing bonds and links which do not exist!
        for k in range(len(self.bridges) - 1, -1, -1):
            if (find_index(self.bridges[k][3], self.coordinates) == None or find_index(self.bridges[k][6],
                                                                                       self.coordinates) == None):
                self.bridges.pop(k)
        ### removing amide bonds used to enlarge the backbone
        for k in range(len(self.bridges) - 1, -1, -1):
            # if ((self.bridges[k][0]=="Amide-like") and (abs(find_index(self.bridges[k][3],self.coordinates)-find_index(self.bridges[k][6],self.coordinates))==1)):
            if ((self.bridges[k][0] == "Amide") and (
                    abs(find_index(self.bridges[k][3], self.coordinates) - find_index(self.bridges[k][6],
                                                                                      self.coordinates)) == 1)):
                self.bridges.pop(k)
        ### dealing with different number of residues in different parts of PDB file
        k = len(self.residues) - len(self.coordinates)
        if (k > 0):
            check = 1  # check=1, if residues in N-terminus are same
            for i in range(min(k, len(self.residues), len(self.coordinates))):
                if (self.residues[i] != self.coordinates[i][2]):
                    check = 0
            if (check == 0):
                for i in range(k):
                    self.coordinates.insert(0, [self.coordinates[0][0] - 1, [], ""])
            else:
                if (len(self.residues) != 0 and len(self.coordinates) != 0):
                    for i in range(len(self.residues) - k, len(self.residues), 1):
                        self.coordinates.append([self.coordinates[len(self.coordinates) - 1][0] + 1, [], ""])
        if (k < 0):
            check = 1  # check=1, if residues in N-terminus are same
            for i in range(k):
                if (self.residues[i] != self.coordinates[i][2]):
                    check = 0
            if (check == 0):
                for i in range(k):
                    self.residues.insert(0, "UNK")
            else:
                for i in range(len(self.residues), len(self.residues) + k, 1):
                    self.residues.append("UNK")
                    ### length check, just for sure
        self.find_length()
        ######### checking gaps

    def check_gaps(self):
        communicate = "\n"
        for k in range(len(self.gaps)):
            if (self.gaps[k][0] > 7):
                communicate += "WARNING!!! In chain " + self.name + " there is a gap of length " + str(
                    self.gaps[k][0] - 1) + " between residues " + str(self.gaps[k][1]) + " and " + str(
                    self.gaps[k][2]) + "\n"
        communicate = communicate[:-1]
        return communicate

    ######### printing data
    def chain_print(self, PDB):
        output_file = open(PDB + "_" + self.name + ".xyz", 'w')
        print(self.check_gaps())
        if (len(self.residues) == self.length):
            for k in range(min(self.length, len(self.coordinates))):
                if (self.coordinates[k][1] != []):
                    output_file.write(str(self.coordinates[k][0]) + " " + str(self.coordinates[k][1][0]) + " " + str(
                        self.coordinates[k][1][1]) + " " + str(self.coordinates[k][1][2]) + " " + str(
                        self.residues[k]) + "\n")
        else:
            for k in range(self.length):
                if (self.coordinates[k][1] != []):
                    output_file.write(str(self.coordinates[k][0]) + " " + str(self.coordinates[k][1][0]) + " " + str(
                        self.coordinates[k][1][1]) + " " + str(self.coordinates[k][1][2]) + " " + str(
                        self.coordinates[k][2]) + "\n")
        output_file.close()
        output_file = open(PDB + "_" + self.name + ".pdb", 'w')
        for k in range(len(self.helix)):
            output_file.write(self.helix[k])
        for k in range(len(self.sheet)):
            output_file.write(self.sheet[k])
        if (len(self.residues) == self.length):
            for k in range(min(self.length, len(self.coordinates))):
                if (self.coordinates[k][1] != []):
                    output_file.write(
                        "ATOM  %(atom)5s  CA  %(resname)3s A%(res)4s    %(x)8s%(y)8s%(z)8s  1.00  1.00           C\n" % {
                            "atom": self.coordinates[k][0], "res": self.coordinates[k][0], "resname": self.residues[k],
                            "x": self.coordinates[k][1][0], "y": self.coordinates[k][1][1],
                            "z": self.coordinates[k][1][2]})
        else:
            for k in range(self.length):
                if (self.coordinates[k][1] != []):
                    output_file.write(
                        "ATOM  %(atom)5s  CA  %(resname)3s A%(res)4s    %(x)8s%(y)8s%(z)8s  1.00  1.00           C\n" % {
                            "atom": self.coordinates[k][0], "res": self.coordinates[k][0],
                            "resname": self.coordinates[k][2], "x": self.coordinates[k][1][0],
                            "y": self.coordinates[k][1][1], "z": self.coordinates[k][1][2]})
        output_file.write("END\n")
        output_file.close()

    ######### printing commands to Wanda's program
    def commands_print(self, PDB, flag=0):
        for k in range(len(self.bridges)):
            if (flag == 1):
                print(self.bridges[k][0] + " ./surfacesMyOrient " + PDB + "_" + self.name + ".xyz " + str(
                    self.bridges[k][3]) + " " + str(self.bridges[k][6]) + " 0 0")
            if (flag == 2):
                print(self.bridges[k][0] + " " + PDB + "_" + self.name + " " + str(self.bridges[k][3]) + " " + str(
                    self.bridges[k][6]))
            else:
                print("./surfacesMyOrient " + PDB + "_" + self.name + ".xyz " + str(self.bridges[k][3]) + " " + str(
                    self.bridges[k][6]) + " 0 0")

    def getBridges(self, filepath):
        result = []
        for k in range(len(self.bridges)):
            result.append((
                self.bridges[k][0],
                self.name.upper(),
                self.bridges[k][3],
                self.bridges[k][6],
                # "./surfacesMyOrient "+filepath+"_"+self.name+".xyz "+str(self.bridges[k][3])+" "+str(self.bridges[k][6])+" 0 0",
            ))
        return result

    def getGaps(self):
        first = True
        result = ''
        for g in self.gaps:
            if not first:
                result += ', '
            first = False
            result += str(g[1]) + '-' + str(g[2])
        return result

    ######### finding distance between residues
    def find_distance(self, res1, res2):
        for k in range(len(self.coordinates)):
            if (self.coordinates[k][0] == res1):
                vec1 = self.coordinates[k][1]
            if (self.coordinates[k][0] == res2):
                vec2 = self.coordinates[k][1]
        print(np.linalg.norm(np.asarray(vec1) - np.asarray(vec2)))

    ######### finding index of first residue with coordinates
    def find_first(self):
        for k in range(len(self.coordinates)):
            if (self.coordinates[k][1] != []):
                return self.coordinates[k][0]


################################ Main part ################################
### search for chains and build chain classes
def run_pdb2Wanda(infile, work_dir, filename):
    with open(infile, 'r') as f:
        lines = []
        for line in f.readlines():
            lines.append(line.strip('\n'))

    chains = []
    for line in lines:
        if (line[0:6] == "SEQRES"):
            if (line[11] not in chains):
                chains.append(line[11])
        if (line[0:4] == "ATOM"):
            if (line[21] not in chains):
                chains.append(line[21])

    ### fill the chain class information
    cross_chain_bridges = []
    result = []
    result_gaps = {}
    result_str_begin = {}
    chain_data = {}
    for k in range(len(chains)):
        chain_data[chains[k]] = Chain(chains[k])
        for line in lines:
            if ((line[0:10] == "REMARK 465") and (line[15:18] in amino_acids) and (
                    line[19] == chains[k]) and isinstance(line[21:26], int)):
                chain_data[chains[k]].add_missing(int(line[21:26]), line[15:18])
            #  if ((line[0:5]=="DBREF") and (line[26:32].strip()=="UNP") and (line[12]==chains[k])):	# gets the taxonomy, but slows down the program
            #   chain_data[chains[k]].get_taxonomy(line[33:41].strip())
            if ((line[0:6] == "SEQRES") and (line[11] == chains[k])):
                chain_data[chains[k]].add_residue(line)
            if ((line[0:5] == "HELIX") and (line[19] == chains[k]) and (line[31] == chains[k])):
                chain_data[chains[k]].add_helix(line)
            if ((line[0:5] == "SHEET") and (line[21] == chains[k]) and (line[32] == chains[k])):
                chain_data[chains[k]].add_sheet(line)
            if ((line[0:6] == "SSBOND") and (line[15] == chains[k]) and (line[29] == chains[k])):
                chain_data[chains[k]].add_bridge(
                    ["SS", line[11:14], "S", int(line[17:21]), line[25:28], "S", int(line[31:35])])
            if ((line[0:6] == "SSBOND") and (line[15] != line[29])):
                cross_chain_bridges.append(
                    ["SS", line[15], line[11:14], "S", int(line[17:21]), line[29], line[25:28], "S", int(line[31:35])])
            if ((line[0:4] == "LINK") and (line[21] == chains[k]) and (line[51] == chains[k])):
                chain_data[chains[k]].add_bridge(
                    ["LINK", line[17:20], line[12:16].strip(), int(line[22:26]), line[47:50], line[42:46].strip(),
                     int(line[52:56])])
            if ((line[0:4] == "LINK") and (line[21] != line[51])):
                cross_chain_bridges.append(
                    ["LINK", line[21], line[17:20], line[12:16].strip(), int(line[22:26]), line[51], line[47:50],
                     line[42:46].strip(), int(line[52:56])])
            if ((line[0:6] == "HETATM") and (line[13:15] == "CA") and (line[21] == chains[k])):
                chain_data[chains[k]].add_coordinate(int(line[22:26]),
                                                     [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                                                     line[17:20])
            if ((line[0:4] == "ATOM") and (line[13:15] == "CA") and (line[21] == chains[k])):
                chain_data[chains[k]].add_coordinate(int(line[22:26]),
                                                     [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                                                     line[17:20])
            if ((line[0:3] == "TER") and (len(chain_data[chains[k]].coordinates) != 0)):
                break
        chain_data[chains[k]].clean()  # clean bridges data - do not comment
        chain_data[chains[k]].chain_print(work_dir + filename)  # save coordinates to .xyz and .pdb file
        # chain_data[chains[k]].commands_print(sys.argv[1])	# print input commands to Wanda's program consistent with files saved line above. Addition of second argument=1 prints as first field the bond type
        result.extend(chain_data[chains[k]].getBridges(work_dir + filename))
        result_gaps[chain_data[chains[k]].name] = chain_data[chains[k]].getGaps()
        result_str_begin[chain_data[chains[k]].name] = chain_data[chains[k]].find_first()
    return (result, result_gaps, result_str_begin)
