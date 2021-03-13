import tarfile
import pickle as _pickle
from os import path, makedirs, listdir

'''Utility methods used by other classes
   '''


def mkDirIfNotExist(dir):
    if not path.exists(dir):
        makedirs(dir)


def mkDirsIfNotExist(dirs):
    for dir in dirs:
        mkDirIfNotExist(dir)


def bz2Files(files, compressed_file, logger=None):
    tar = tarfile.open(compressed_file, "w:bz2")
    for name in files:
        try:
            tar.add(name, path.basename(name))
        except OSError:
            if logger:
                logger.error("Missing file while bz2: " + name)
            else:
                print("Missing file while bz2: " + name)
    tar.close()


def bz2Folder(folder, compressed_file, logger=None):
    mkDirIfNotExist(path.dirname(compressed_file))
    tar = tarfile.open(compressed_file, "w:bz2")
    for name in listdir(folder):
        try:
            tar.add(path.join(folder, name), path.basename(name))
        except OSError:
            if logger:
                logger.error("Missing file while bz2: " + name)
            else:
                print("Missing file while bz2: " + name)
    tar.close()


def arraytostring(ar):
    i = []
    for e in ar:
        i.append("%d-%d" % (e))
    return ", ".join(i)


# Gzip - necessary for JSMOL viewer
def gzipPDB(code, chain, path, suffix=''):
    return _gzip(path + '/' + code + '_' + chain + suffix + '.pdb')


# Gzip
def _gzip(file):
    import gzip
    import shutil
    with open(file, 'rb') as f_in, gzip.open(file + '.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    return


def read_from_bin_file(filename):
    # Convert digital data to binary format
    with open(filename, 'rb') as file:
        return _pickle.load(file)


def write_file(data, filename):
    # Convert binary data to proper format and write it on Hard Disk
    with open(filename, 'wb') as file:
        file.write(data)
