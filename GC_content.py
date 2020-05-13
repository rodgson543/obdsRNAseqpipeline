# -*- coding: utf-8 -*-
"""
Created on Tue May  5 15:11:24 2020

@author: spandey
"""

# import module
# initialise the command line parser
# add arguments to the command line parser
# parse the command line arguments
# check the validity of the arguments parsed
# use those arguments in your program

import argparse

parser = argparse.ArgumentParser(description='Calculating GC content.')
parser.add_argument('-f', '--file', type=str, help='input file (fasta format)')

args = parser.parse_args()

#open test file
gc = 0
len_seq = 0
with open(args.file, "r") as f:
    for line in f:
        for character in line:
            if character == "A" or character == "G" or character == "C" or character == "T":
                len_seq += 1
                if character == "G" or character == "C":
                    gc += 1
gc_frac = float(gc/len_seq)
result = gc_frac * 100
print(result)