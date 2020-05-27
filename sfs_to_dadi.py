#!/usr/bin/python

import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

################################################################################
# this script takes output sfs.txt (from using SiteFrequencySpectra.py) and
# reformats it for input into dadiAnalysis.py and the other modeling scripts
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
    parser = argparse.ArgumentParser(description='format converter to dadi input')
    parser.add_argument("sfs_in", help="sfs_input.txt",
        action=FullPaths,
        type=is_file),
    parser.add_argument("num_isolates", help="<number_isolates(witout_outgroup)>")
    return parser.parse_args()

def read_sfs(sfs_in):
    sfs_list = []
    infile = open(sfs_in,"r")
    for i, line in enumerate(infile):
        line = line.strip().split("\t")
        if i > 0:
            syn = line[1]
            sfs_list.append(syn)
    return sfs_list
    
def write_sfs(sfs_ds,n):
    sfs_out = open("sfs_dadi_in.txt","w")
    line_start_int = int(n) + 1
    sfs_out.write(str(line_start_int)+" unfolded\n0 ")
    for item in sfs_ds:
        sfs_out.write(item+" ")
    sfs_out.write("\n")
    sfs_out.close()    

args = get_args()

sfs = read_sfs(args.sfs_in)
write_sfs(sfs,args.num_isolates)

