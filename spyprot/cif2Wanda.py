import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

#######
# Author: Pawel Dabrowski-Tumanski
# p.dabrowski@cent.uw.edu.pl
# Date 02.10.2015; version from 28.01.2016
# Script reads the pdb file and gets the chain information. Each chain forms 
#   individual class. In each class the coordinates of CA atoms, missing atoms,
#   residues, SS-bonds and links are stored.
# Script converts pdb file to xyz file (one for each chain), in 5-column format
#   (with residue as the last column), fills the gap with straight line and 
#   prints warning if the gap is longer than 6 residues.
# Finally it gives the commands for Wanda's program for finding the crossing of
#   the surfaces.
#######

############ Chain Class definition and possible amino acids ##################
AMINOACIDS = ['ALA', 'ARG', 'ASN', 'ASP', 'CIT', 'CYS', 'GLU', 'GLN', 'GLY', 'HEY', 'HIS',
              'HYP', 'ILE', 'LEU', 'LYS', 'MET', 'ORN', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
              'SEC', 'PYL', 'ASX', 'GLX', 'XLE', 'XAA', 'MSE', 'FGL', 'LLP', 'SAC', 'PCA', 'MEN', 'CSB',
              'HTR', 'PTR', 'SCE', 'M3L', 'OCS', 'KCX', 'SEB', 'MLY', 'CSW', 'TPO', 'SEP', 'AYA', 'TRN',
              'DAL', 'DAR', 'DSG', 'DAS', 'DCY', 'DGN', 'DGL', 'DHI', 'DIL', 'DLE', 'DLY', 'MED', 'DPN',
              'DPR', 'DSN', 'DTH', 'DTR', 'DTY', 'DVA']


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
        self.length = None

    def find_length(self):
        self.length = max(len(self.residues), len(self.coordinates))

    def add_residue(self, resid):
        if resid in AMINOACIDS:
            self.residues.append(resid)
        self.find_length()  # relikt z pdb2wanda

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
                #                res=self.missing[find_index(index-diff+k+1,self.missing)][1]
                #                if (res==None):
                res = "XXX"
                self.coordinates.append([index - diff + k + 1, [x, y, z], res])
            self.coordinates.append([index, coordinate, residue])

    ######### characterize the bond type 
    def bond_type(self, res1, atom1, Nend, res2, atom2, Cend):
        a1, a2 = sorted((atom1[0], atom2[0]))
        if atom1 == atom2 == 'SG':
            return 'SS'
        elif a1 == "C":
            if a2 == "N":
                return "Amide"
            elif a2 == "O":
                return "Ester"
            elif a2 == "S":
                return "Thioester"
        return 'Others'

    ######## check, whether it is N- or C-end
    def N_end(self, index):
        return index == self.coordinates[0][0]

    def C_end(self, index):
        return index == self.coordinates[len(self.coordinates) - 1][0]

    ######### cleaning data
    def clean(self):
        ### adding C-end
        if len(self.coordinates) != 0:
            for missing in self.missing:
                if (self.coordinates[len(self.coordinates) - 1][0] < missing[0]):
                    self.coordinates.append([missing[0], [], missing[1]])
        else:
            for missing in self.missing:
                self.coordinates.append([missing[0], [], missing[1]])
                ### clearing double CA atoms
        for k in range(len(self.coordinates) - 1, 0, -1):
            if (self.coordinates[k][0] == self.coordinates[k - 1][0]):
                self.coordinates.pop(k)
        ### clearing non-protein links
        for k in range(len(self.bridges) - 1, -1, -1):
            if ((self.bridges[k][1] not in AMINOACIDS) or (self.bridges[k][4] not in AMINOACIDS)):
                self.bridges.pop(k)
        ### clearing to small loops
        for k in range(len(self.bridges) - 1, -1, -1):
            if (abs(self.bridges[k][3] - self.bridges[k][6]) < 5):
                self.bridges.pop(k)
        ### defining bond type
        for bridge in self.bridges:
            if (bridge[0] == "LINK"):
                bridge[0] = self.bond_type(bridge[1], bridge[2], self.N_end(bridge[3]),
                                           bridge[4], bridge[5], self.C_end(bridge[6]))
        ### removing bonds and links which do not exist!
        for k in range(len(self.bridges) - 1, -1, -1):
            if (find_index(self.bridges[k][3], self.coordinates) == None
                    or find_index(self.bridges[k][6], self.coordinates) == None):
                self.bridges.pop(k)
        ### removing amide bonds used to enlarge the backbone
        for k in range(len(self.bridges) - 1, -1, -1):
            index_diff = find_index(self.bridges[k][3], self.coordinates) \
                         - find_index(self.bridges[k][6], self.coordinates)
            if self.bridges[k][0] == "Amide" and abs(index_diff) == 1:
                self.bridges.pop(k)
        ### dealing with different number of residues in different parts of PDB file
        k = len(self.residues) - len(self.coordinates)
        if k > 0:
            check = 1  # check=1, if residues in N-terminus are same
            for i in range(min(k, len(self.residues), len(self.coordinates))):
                if (self.residues[i] != self.coordinates[i][2]):
                    check = 0
            if check == 0:
                for i in range(k):
                    self.coordinates.insert(0, [self.coordinates[0][0] - 1, [], ""])
            else:
                if len(self.residues) != 0 and len(self.coordinates) != 0:
                    for i in range(len(self.residues) - k, len(self.residues), 1):
                        self.coordinates.append([self.coordinates[len(self.coordinates) - 1][0] + 1, [], ""])
        if k < 0:
            check = 1  # check=1, if residues in N-terminus are same
            for i in range(k):
                if (self.residues[i] != self.coordinates[i][2]):
                    check = 0
            if check == 0:
                for i in range(k):
                    self.residues.insert(0, "UNK")
            else:
                for i in range(len(self.residues), len(self.residues) + k, 1):
                    self.residues.append("UNK")
                    ### length check, just for sure
        self.find_length()

        ######### checking gaps

    def check_gaps(self):
        communicate = ""
        for gap in self.gaps:
            if gap[0] > 7:
                communicate += "\nWARNING!!! In chain {} there is a gap of " \
                               "length {} between residues {} and {}\n".format( \
                    self.name, str(gap[0] - 1), str(gap[1]), gap[2])
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

    def getBridges(self):
        result = []
        for bridge in self.bridges:
            result.append((bridge[0], self.name.upper(), bridge[3], bridge[6]))
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
def run_cif2Wanda(infile, work_dir, outfile):
    # ignored search for secondary structures
    # infile -- file name or handle
    cifdict = MMCIF2Dict(infile)

    ### fill the chain class information
    chains = list(set(cifdict['_struct_asym.id']))
    chain_data = {chain: Chain(chain) for chain in chains}
    cross_chain_bridges = []

    # fill Chain objects with aminoacid sequence
    for chain, resid_seq, seqid, resid_struct in zip(cifdict['_pdbx_poly_seq_scheme.asym_id'],
                                                     cifdict['_pdbx_poly_seq_scheme.mon_id'],
                                                     cifdict['_pdbx_poly_seq_scheme.seq_id'],
                                                     cifdict['_pdbx_poly_seq_scheme.pdb_mon_id']):
        chain_data[chain].add_residue(resid_seq)
        if resid_struct == '?':
            chain_data[chain].add_missing(int(seqid), resid_seq)

    # fill Chain objects with ssbonds and other links
    if type(cifdict['_struct_conn.ptnr1_label_asym_id']) == list:
        for chain_beg, chain_end, seqid_beg, seqid_end, resid_beg, resid_end, \
            atom_beg, atom_end in zip(cifdict['_struct_conn.ptnr1_label_asym_id'],
                                      cifdict['_struct_conn.ptnr2_label_asym_id'],
                                      cifdict['_struct_conn.ptnr1_label_seq_id'],
                                      cifdict['_struct_conn.ptnr2_label_seq_id'],
                                      cifdict['_struct_conn.ptnr1_label_comp_id'],
                                      cifdict['_struct_conn.ptnr2_label_comp_id'],
                                      cifdict['_struct_conn.ptnr1_label_atom_id'],
                                      cifdict['_struct_conn.ptnr2_label_atom_id']):
            if chain_beg == chain_end:
                chain_data[chain_beg].add_bridge(["LINK", resid_beg, atom_beg, int(seqid_beg),
                                                  resid_end, atom_end, int(seqid_end)])
            else:
                cross_chain_bridges.append(["LINK", chain_beg, resid_beg, atom_beg, int(seqid_beg),
                                            chain_end, resid_end, atom_end, int(seqid_end)])
    elif type(cifdict['_struct_conn.ptnr1_label_asym_id']) == str:
        chain_beg, chain_end = cifdict['_struct_conn.ptnr1_label_asym_id'], cifdict['_struct_conn.ptnr2_label_asym_id']
        seqid_beg, seqid_end = cifdict['_struct_conn.ptnr1_label_seq_id'], cifdict['_struct_conn.ptnr2_label_seq_id']
        resid_beg, resid_end = cifdict['_struct_conn.ptnr1_label_comp_id'], cifdict['_struct_conn.ptnr2_label_comp_id']
        atom_beg, atom_end = cifdict['_struct_conn.ptnr1_label_atom_id'], cifdict['_struct_conn.ptnr2_label_atom_id']
        if chain_beg == chain_end:
            chain_data[chain_beg].add_bridge(["LINK", resid_beg, atom_beg, int(seqid_beg),
                                              resid_end, atom_end, int(seqid_end)])
        else:
            cross_chain_bridges.append(["LINK", chain_beg, resid_beg, atom_beg, int(seqid_beg),
                                        chain_end, resid_end, atom_end, int(seqid_end)])

    # fill Chain objects with coordinates
    for chain, seqid, resid, atom, x, y, z in zip(cifdict['_atom_site.label_asym_id'],
                                                  cifdict['_atom_site.label_seq_id'],
                                                  cifdict['_atom_site.label_comp_id'],
                                                  cifdict['_atom_site.label_atom_id'],
                                                  cifdict['_atom_site.Cartn_x'],
                                                  cifdict['_atom_site.Cartn_y'],
                                                  cifdict['_atom_site.Cartn_z']):
        if atom == 'CA':
            chain_data[chain].add_coordinate(int(seqid), [float(x), float(y), float(z)], resid)

    ### final collecting
    result = []
    result_gaps = {}
    result_str_begin = {}
    for chain in chains:
        chain_data[chain].clean()  # clean bridges data - do not comment
        chain_data[chain].chain_print(work_dir + outfile)  # save coordinates to .xyz and .pdb file
        # chain_data[chain].commands_print(sys.argv[1])	# print input commands to Wanda's program consistent with files saved line above. Addition of second argument=1 prints as first field the bond type
        result.extend(chain_data[chain].getBridges())
        result_gaps[chain_data[chain].name] = chain_data[chain].getGaps()
        result_str_begin[chain_data[chain].name] = chain_data[chain].find_first()
    return result, result_gaps, result_str_begin
