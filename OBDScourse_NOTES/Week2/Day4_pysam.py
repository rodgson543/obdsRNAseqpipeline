#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:24:29 2020

@author: rhodgson
"""


#Day4: Using pysam to convert a bam to bed using pysam 
#First going to convert sam file in CGAT cluster week 2 obds folder to BAM
#Use samtools in cgat cluster to convert 
#samtools view -S -b ERR1755082.test.sam > ERR1755082.test.bam #converts to bam
#samtools sort ERR1755082.test.bam -o ERR1755082.test.sorted.bam #sorts bam
#samtools index ERR1755082.test.sorted.bam #indexes bam
#copy file locally (need the index file too) ~/OBDStestdata

#Going to do all this work on a new branch
#git checkout master #switch to master and git pull
#git checkout -b Rose #create/switch to new branch
#git push -u origin Rose

#Bed files have 3 compulsory fields: chromosome, chromosome start and chromosome end

#Now want to convert from bam to bed using pysam
import pysam 

#Define input file #rb is used when the input file is not text
#Define output file as a bedfile, want to use the input file as a template 

with open("/Users/rhodgson/OBDStestdata/ERR1755082.new.bed", "w") as output:
    #Making an object of the class AlignmentFile below
    infile = pysam.AlignmentFile("/Users/rhodgson/OBDStestdata/ERR1755082.test.sorted.bam", "rb")
    #The iter will grab all data from the file, fetch is a method that belongs to the AlignmentFile class
    iter = infile.fetch() 
    #For each of the alignments, please do this(an alignment is each line of the file)
    #We now want to look at each of the lines of bam and bed and what each has to convert it
    #We want to only take lines that are properly mapped (both of the reads)
    for aln in iter:
        if aln.is_proper_pair:
            chr = aln.reference_name
            start = aln.reference_start
            end = aln.reference_end
        else:
            continue
#Now want to get rid of unmapped region        
        
        print(chr)
        print(start)
        print(end)
        print(type(aln))
        output.write(f'{chr}\t{start}\t{end}\n')



#This code is to run from command line


import pysam
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

with open(args.bed, "w") as output:
    #Making an object of the class AlignmentFile below
    infile = pysam.AlignmentFile(args.bam, "rb")
    #The iter will grab all data from the file, fetch is a method that belongs to the AlignmentFile class
    iter = infile.fetch() 
    #For each of the alignments, please do this(an alignment is each line of the file)
    #We now want to look at each of the lines of bam and bed and what each has to convert it
    #We want to only take lines that are properly mapped (both of the reads)
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

 
        
        
#We are now going to try and run using stdin stdout



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
    



with open(args.bed, "w") as output:
    #For each of the alignments, please do this(an alignment is each line of the file)
    #We now want to look at each of the lines of bam and bed and what each has to convert it
    #We want to only take lines that are properly mapped (both of the reads)
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

    

#Now want to 
#Fixed the stdout
#cat /Users/rhodgson/OBDStestdata/ERR1755082.test.sorted.bam | python /Users/rhodgson/GitHub/OBDS_Training_Apr_2020/BAMtoBEDpysam.py -i - -o - | wc -l > counting.txt
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
    

        
