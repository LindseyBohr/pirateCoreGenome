#!/usr/bin/python

import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

################################################################################
# This script creates an aligned core genome from a specified list of isolates
# so you don't have to rerun roary or pirate
# Have to use subsample_output.pl to get subsampled gene_output file 
# Output includes gff file
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
    parser = argparse.ArgumentParser(description='Create core gene aln')
    parser.add_argument("gene_info", help="subsampled_pirate_output", 
	action=FullPaths,
        type=is_file)
    parser.add_argument("sequences",
        help="Directory with .fasta files of gene sequences in nucleotides",
        action=FullPaths, type=is_dir)
    parser.add_argument("threshold",
        help="gene frequency to be called in core (ie 95 or 100)")
    return parser.parse_args()

def getCoreGenesList(infile,threshold):
    """Build core genome list"""
    infile = open(infile,"r")
    groupsTuples = []
    idDict = {}
    for i, line in enumerate(infile):
        if i == 0:
            idLine = line.strip().split("\t")
            idList = idLine[22:]
            for item in idList:
                idDict[item] = idList.index(item)
            ### grab isolate index info
        if i > 0:
            line = line.strip()
            line = line.split("\t")
            order = line[20]
            gene_product = line[2] ### this is the gene_product to put in .gff3
            len_of_gene = line[18] ### this is the max len
            num_genomes = line[6]
            num_fission_loci = int(line[10])
            allele_name = line[0]
            gene_family = line[1]
            num_threshold = int(threshold) * .01 * len(idList)
            total_num_iso = len(idList)
            num_threshold_rounded = int(round(num_threshold))
            if num_fission_loci == 0:   ### if theres any genomes that contain a fission version of this gene, skip the gene entirely to avoid accumulating gaps
                if int(num_genomes) >= num_threshold_rounded: ### in 100% of clade 1
                    groupProteins = line[22:]
                    groupProteins=[make_protein_fancy(protein,len_of_gene) for protein in groupProteins]
                    groupsTuples.append((gene_family,groupProteins,gene_product))
    return groupsTuples, idDict, num_threshold_rounded, total_num_iso

def make_protein_fancy(protein,length):
    """Deals with fission, duplicates, etc in gene family index"""
    if not protein:
        return "none_"+length
    elif ":" and ";" in protein: ### this and the other colon one are fusions and don't actually get used with the above num_fission_loci if statement
        return make_protein_fancy(protein.split(";")[1],length)
    elif ";" in protein:
        return protein.split(";")[0]
    elif ":" in protein:
        protein=protein[1:]
        return protein.split(":")[0]
    return protein


def actually_concatenate_aln(groupsTuple, idDict, directory, total_num_iso):
    ### go through groupsTuples
    seqlist = [None]*total_num_iso
    gffpos = 1
    gffdata = []
    for i,j in idDict.items():
        seqlist[j] = SeqRecord(Seq(""),id=i)
    for geneTuple in groupsTuple:
        gffdone = 0
        genelength = 0 ##for isolates where it isn't found
        gene,isolates,gene_product = geneTuple
        foundlist = [None]*total_num_iso ##To make sure we add dashes for 'none' isolates
        for record in SeqIO.parse(directory+"/"+gene+".nucleotide.fasta","fasta"):
            if record.id in isolates:
                isoid="_".join(record.id.split("_")[:-1])
                foundlist[idDict[isoid]] = record.id
                seqlist[idDict[isoid]].seq += record.seq
                if not gffdone:  #gff data
                    gffdone = 1
                    genelength = len(record)
#                    gffpos += 1
                    temppos = gffpos + genelength
                    gffdata.append([gene_product,gffpos,temppos])
                    gffpos = temppos
        for num,iso in enumerate(isolates):
            if iso not in foundlist:
                print(num,iso, gene)
                seqlist[num].seq += "-"*genelength #more accurate than using the gene_len_est from the gene_info.tsv
    write_out_gff(gffdata,gffpos)### write out line to gff file? or store values in gff data structure and write out in a dif function
    SeqIO.write(seqlist,"newCoreGenomeAln.fasta","fasta")
    pass

def write_out_gff(gene_ds,aln_length):
    """writes out core genome gff file"""
    outGff = open("new_core.gff","w")
    outGff.write("##gff-version 3\n##sequence-region core_genome 1 "+str(aln_length)+"\n")
    for gene in gene_ds: # list of lists [[gene_product,start,stop],[g2,s2,s2], [etc]]
        product = gene[0]
        start = gene[1]
        stop = gene[2]
        outGff.write("core_genome\tRoary\tCDS\t{0}\t{1}\t.\t+\t0\tID={2};gene={2}:locus_tag={2}\n".format(start,stop,product))
    outGff.write("##FASTA\n")
    outGff.close()
    return
        

args = get_args()

groupsTuples, idDict, num_iso, total_num_iso  = getCoreGenesList(args.gene_info, args.threshold)
seqIDs = actually_concatenate_aln(groupsTuples, idDict, args.sequences, total_num_iso)

# input: list of isolates, gene presence absence matrix, dir to 
# pangenome aln-ed genes
# output: core gene aln of concatenated core genes by isolate, gff file
# output: number of genes in core, length of core

