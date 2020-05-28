#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:27:41 2020

@author: rhodgson
"""
import pysam
import sys 
import logging as L
L.basicConfig(filename='/Users/rhodgson/example.log', filemode= 'w', level= L.DEBUG)

#Import argparse
import argparse
#define a parser that takes an input and output file
#First defining parser - initialising
#next need to add a description
parser = argparse.ArgumentParser(description='Converting bam to bed with pysam')
#add the argument - for bar is name/flag then need type - so string as it's just the name of your file. help is what you will add, ie inputfile
parser.add_argument('-i', '--bam', type=str, help ='inputfile')
parser.add_argument('-o', '--bed', type=str, help ='outputfile')
args = parser.parse_args()



#stdin and stdout - had to add the end of file here (until_eof=True) 
#Could also do the same here for stdout - change args.bed to "-" 
#gzip too - writing to a compressed file 
#Add option (do you want to truntcate to 1st bit)
if args.bam == "-" :
    input_bam = sys.stdin
    infile = pysam.AlignmentFile(input_bam, "rb")
    iter = infile.fetch(until_eof=True)
else:
    input_bam = args.bam
    infile = pysam.AlignmentFile(input_bam, "rb")
    iter = infile.fetch()


#Figuring out the stdout bit - surely we'll have to write the out 
if args.bed == "-" :
    output_bed = sys.stdout
else:
    output_bed = open (args.bed, "w")
    
    
with output_bed as output:
    for aln in iter:
        if aln.is_proper_pair:
            chr = aln.reference_name
            start = aln.reference_start
            end = aln.reference_end
            output.write(f'{chr}\t{start}\t{end}\n')
            L.info(chr)
            L.info(start)
            L.info(end)
        else:
            L.warn('not a proper pair')
            continue
    
        