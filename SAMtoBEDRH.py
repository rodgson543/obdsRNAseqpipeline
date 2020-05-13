#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:30:07 2020

@author: rhodgson

"""
#This script will run the following pipeline for everyone, you just need to type
#python /Users/rhodgson/GitHub/OBDS_Training_Apr_2020/SAMtoBEDRH.py -i ERR1755082.test.sam -o ERR1755082NEW.test.bed


#Import argparse
import argparse
#define a parser that takes an input and output file
#First defining parser - initialising
#next need to add a description
parser = argparse.ArgumentParser(description='Converting sam to bed')
#add the argument - for bar is name/flag then need type - so string as it's just the name of your file. help is what you will add, ie inputfile
parser.add_argument('-i', '--sam', type=str, help ='inputfile')
parser.add_argument('-o', '--bed', type=str, help ='outputfile')
args = parser.parse_args()


with open(args.bed, "w") as output:
    with open(args.sam, "r") as sam:
        lastid = -1
        for line in sam:
            if line.startswith('@'):
                continue 
            splitline= (line.split("\t"))
            id = splitline[0]
            if lastid == id:
                continue
            chr = str(splitline[2])  
            if chr == "*":
                continue
            start = int(splitline[3])-1
            start_mate = int(splitline[7])-1
            if start < start_mate:
                start_bed = start
            else: 
                start_bed = start_mate
            stop_bed = start_bed+ abs(int(splitline[8]))
            
            print(chr, start_bed, stop_bed)
            output.write(f'''{chr}\t{start_bed}\t{stop_bed}\n''')

            lastid = id

           
        


#OK we're going to do it the hard way: Process the sam file format to 
#Supply the SAM file name on the command line using –i or --input


#Supply the BED file name on the command line using –o or –output
#Format your output using F-strings
#Provide a command line argument to pad the intervals in the bed file
#Output a file with the coordinates of the sequenced fragments