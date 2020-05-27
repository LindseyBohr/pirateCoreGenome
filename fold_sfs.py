#!/usr/bin/python

import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

################################################################################
# this script takes output sfs_dadi_in.txt (from using sfs_to_dadi.py) and
# folds it for input into dadiAnalysis.py and the other modeling scripts
################################################################################

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='unfolded dadi input to folded')
    parser.add_argument("sfs_in", help="sfs_dadi_in.txt",
        action=FullPaths,
        type=is_file),
    return parser.parse_args()

def read_sfs(sfs_in):
    infile = open(sfs_in,"r")
    for i, line in enumerate(infile):
        line = line.strip().split()
        if i == 0:
            numiso = line[0]
        if i > 0:
            sfs = line[1:] # remove the zero from the beginning since we don't fold that
    return sfs, numiso

def fold_sfs(sfs): # takes list of sfs values starting with singletons
    new_sfs = []
    for i, value in enumerate(sfs,1):
        if i <= (len(sfs)/2):
            tempSum = int(value) + int(sfs[-i])
            new_sfs.append(tempSum)
    if len(sfs)%2 == 1:
        new_sfs.append(sfs[len(sfs)/2+1])
    return new_sfs
   
def write_sfs(sfs_ds,n):
    sfs_out = open("folded_sfs_dadi_in.txt","w")
    line_start_int = n
    sfs_out.write(str(line_start_int)+" folded\n0 ")
    for item in sfs_ds:
        sfs_out.write(str(item)+" ")
    sfs_out.write("\n")
    sfs_out.close()    

args = get_args()

sfs, numiso = read_sfs(args.sfs_in)
sfs = fold_sfs(sfs)
write_sfs(sfs, numiso)

