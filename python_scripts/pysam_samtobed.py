#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:15:28 2020

@author: charlotteb
"""

import pysam
import logging as L
L.basicConfig(filename='example.log', level=L.DEBUG)
#logging.basicConfig(filename='example.log', filemode='w', level=logging.DEBUG)

import argparse
parser = argparse.ArgumentParser(description='pysam - bam to bed file')

parser.add_argument('-i', '--bam', type=str, help='input file (bam format)')
parser.add_argument('-o', '--bed', type=str, help='output file (bed format)')

args = parser.parse_args()




#define output and input input full path, rb means its a bam file, with open always has colon afterwards
with open (args.bed, "w") as output:
    infile = pysam.AlignmentFile(args.bam, "rb")
    #the iter will grab all alignments. fetch belongs to the alignment file object
    #the class is AlignmentFile 
    iter = infile.fetch()
=======
print(args)
print(args.bam)

#This is to prompt the user to feed in the input file directly from command line
if args.bam == "-" :
    input_bam = sys.stdin
    infile = pysam.AlignmentFile(input_bam, "rb")
    iter = infile.fetch(until_eof = True)
else:
    input_bam = args.bam
    infile = pysam.AlignmentFile(input_bam, "rb")
    iter = infile.fetch()

#define output and input input full path, rb means its a bam file, with open always has colon afterwards
with open(args.bed, "w") as output:
    #the iter will grab all alignments. fetch belongs to the alignment file object
    #the class is AlignmentFile    

    #for each of the alignments do this
    for aln in iter:
        if aln.is_proper_pair:
            chrom = aln.reference_name
            bed_start = aln.reference_start
            bed_end = aln.reference_end
            output.write(f'{chrom}\t{bed_start}\t{bed_end}\n')
            L.info(chrom)
            L.info(bed_start)
            L.info(bed_end) 
        else:
            L.warn('not a proper pair')
            continue
  

      
 #are we doing beginning of the read and end of read for each pair
# or the beginning and end of the fragment  
# where both reads are properly mapped we will write the whole sequence

#Each iter returns a AlignedSegment object which represents a single read along with its fields and optional tags:
#we now have an object
#we want to consider bed file format:
#chrom/chromstart/stop
#how will we deal with paired reads? 
#

#example code
# pairedreads = pysam.AlignmentFile("allpaired.bam", "wb", template=samfile)
# for read in samfile.fetch():
#      if read.is_paired:
#              pairedreads.write(read)


#BE AWARE Coordinates in pysam are always 0-based (following the python convention). 
#SAM text files use 1-based coordinates.