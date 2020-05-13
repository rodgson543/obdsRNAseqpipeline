#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:30:07 2020

@author: rhodgson
"""


#Copy the SAM file from shared/week2/ERRetc

#Now write script to convert the SAM file to a BED file
#comparing SAM file format to BED file format
#For each line of my sam file, read it and convert to output line (use for loop and return)
import argparse

parser = argparse.ArgumentParser(description='sam to bed file')
parser.add_argument('-i', '--sam', type=str, help='input file (sam format)')
parser.add_argument('-o', '--bed', type=str, help='output file (bed format)')

args = parser.parse_args()

with open(args.bed, "w") as bed:
    #We open a file in write mode and connect it to the loop that follows
    
    with open(args.sam, "r") as sam:
        lastid = -1
        for line in sam:
    #Now need to look within the file and if it contains the @ symbol then continue
            if line.startswith('@'):
                continue
            splitline = line.split('\t')
            id = splitline[0] #Use ID (first column of sam file) as unique variable for each row  
            if lastid == id: #to remove duplicates, drop iteration if ID = previousID 
                continue
            #chr needs to go in first column of bed file
            chrom = str(splitline[2])
            if chrom == "*":
                 continue #this if command ensures that if chrom contains * the 
                 #rest of the loop is skipped and you move to the next iteration
            start = int(splitline[3])-1
            start_mate = int(splitline[7])-1
            if start < start_mate:
                start_bed = start
            else:
                start_bed = start_mate
            
            stop_bed = start_bed+abs(int(splitline[8])) 
            #print(chrom,start_bed, stop_bed)
            
            bed.write(f'''{chrom}\t{start_bed}\t{stop_bed}\n''')
            #the line above writes to your output file. 
            #\n to ensure it prints each new sequence to a new line 
            #\T to ensure its delimited 
            
            lastid = id # assigns current id as new value to lastid 
            

            
#if you take sequence and you have paired reads
#the pair will be somewhere else in the file
#all you need to do is take the start and add the tlength
#only doing paired ones find the end position of the chromosome: so for paired alignment need to count the characters and add to the staert sit
  
# convert every SAM line into a BED line
#If PNEXT >POS(4) then its a forward read -TLEN
# Is POS>PNEXT then its reverse read
# Then identify duplicates / multimapping downstream 

#OK we're going to do it the hard way: Process the sam file format to 
#Supply the SAM file name on the command line using –i or --input


#Supply the BED file name on the command line using –o or –output
#Format your output using F-strings
#Provide a command line argument to pad the intervals in the bed file
#Output a file with the coordinates of the sequenced fragments