#!/usr/bin/python

import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

################################################################################
# this script takes output sfs.txt (from using SiteFrequencySpectra.py) and
# folds all of it for plotting sSNPs and nsSNPS in R 
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
    return parser.parse_args()

def read_sfs(sfs_in):
    syn_sfs_list = []
    non_syn_sfs_list = []
    infile = open(sfs_in,"r")
    for i, line in enumerate(infile):
        line = line.strip().split("\t")
        if i > 0:
            syn = line[1]
            non_syn = line [2]
            syn_sfs_list.append(syn)
            non_syn_sfs_list.append(non_syn)
    return syn_sfs_list, non_syn_sfs_list

def fold_sfs(sfs): # takes list of sfs values starting with singletons
    new_sfs = []
    for i, value in enumerate(sfs,1):
        if i <= (len(sfs)/2):
            tempSum = int(value) + int(sfs[-i])
            new_sfs.append(tempSum)
    if len(sfs)%2 == 1:
        new_sfs.append(sfs[len(sfs)/2+1])
    return new_sfs
    
def write_sfs(syn, non):
    sfs_out = open("sfs_folded.txt","w")
    sfs_out.write("Frequency\tSynonymous\tNonsynonymous\n")
    freq = 1
    for i,item in enumerate(syn):
        sfs_out.write(str(freq)+"\t"+str(item)+"\t"+str(non[i])+"\n")
        freq += 1
    sfs_out.write("\n")
    sfs_out.close()    

args = get_args()

syn_sfs, non_syn_sfs = read_sfs(args.sfs_in)
folded_syn = fold_sfs(syn_sfs)
folded_non = fold_sfs(non_syn_sfs)
write_sfs(folded_syn,folded_non)

