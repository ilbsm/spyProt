from re import compile
import gzip
import bz2


def convertXYZtoPDB(data_xyz, output_gz, start_idx='', stop_idx=''):
    '''Generate a simple PDB file from a XYZ format containing the 3D coordinates of residues
       '''
    converted = ""
    if start_idx != '': start_idx = int(start_idx)
    if stop_idx != '': stop_idx = int(stop_idx)
    if start_idx != '' and stop_idx != '':
        fragment = True
    else:
        fragment = False

    xyz_atm = compile(r"^\s*(?P<resid>[0-9]*)\s+(?P<x>-?\d+\.\d*)\s+(?P<y>-?\d+\.\d*)\s+(?P<z>-?\d+\.\d*).*$")
    for line in data_xyz:
        data = xyz_atm.match(line)
        if data:
            resid = int(data.groups()[0])
            x = float(data.groups()[1])
            y = float(data.groups()[2])
            z = float(data.groups()[3])
            if fragment:
                if resid >= start_idx and resid <= stop_idx:
                    converted += "ATOM%7d  CA  GLY A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n" % (
                        resid, resid, x, y, z)
            else:
                converted += "ATOM%7d  CA  GLY A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n" % (
                    resid, resid, x, y, z)

    fw = gzip.open(output_gz, 'wb')
    fw.write(converted)
    fw.close()

    return converted


def grepChain(filename, chain, bzip=True):
    atm = compile(r"^ATOM.{9}CA.{6}" + chain + ".*$|^END|^MODEL")
    if bzip:
        f = bz2.BZ2File(filename, "rb")
    else:
        f = open(filename, "r")
    lines = f.readlines()
    f.close()
    out_lines = []

    for line in lines:
        if atm.match(line):
            out_lines.append(line)

    if bzip:
        f = bz2.BZ2File(filename, "wb")
    else:
        f = open(filename, "w")

    f.write("".join(out_lines))
    f.close()


def getSubchain(filename, start_idx, stop_idx, bzip=True, fmt='pdb'):
    if bzip:
        f = bz2.BZ2File(filename, "rb")
        filetype = fmt
    else:
        f = open(filename, "r")
        filetype = fmt

    lines = f.readlines()
    f.close()
    out_lines = []

    if filetype == 'pdb':
        atm = compile(r"^ATOM.{9}CA.*$")
        mdl = compile(r"^END|^MODEL")

        for line in lines:
            if mdl.match(line):
                out_lines.append(line)
            elif atm.match(line):
                i = int(line[22:26])
                if start_idx <= i <= stop_idx:
                    out_lines.append(line)
    elif filetype == 'xyz':
        for line in lines:
            if "t" in line:
                out_lines.append(line)
            else:
                d = line.split()
                i = int(d[0])
                if start_idx <= i <= stop_idx:
                    out_lines.append(line)
    if bzip:
        f = bz2.BZ2File(filename, "wb")
    else:
        f = open(filename, "w")
    f.write("".join(out_lines))
    f.close()
