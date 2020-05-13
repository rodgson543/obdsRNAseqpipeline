#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 10:21:58 2020

@author: rhodgson
"""

#This pipeline is to run rnaseq pipelines from the fastq file from all 
#I think I will create my github repository here 


#rsync -a /Users/rhodgson/GitHub/OBDS_Training_Apr_2020/rnaseq_pipeline.py rose@cgatui.imm.ox.ac.uk:/ifs/obds-training/apr20/rose/pipelines/rnaseqpipeline
#Now import section
import sys
from ruffus import *
from cgatcore import pipeline as P
import gzip


#write parameters
Params = P.get_parameters("pipeline.yml")

#First going to do the fastqc on the fastq files 
#This will create fastqc.html files and fastqc.zip
#It also feeds evertyhing into a fastqc.zip file

@transform('*.fastq.gz', suffix('.fastq.gz'), '_fastqc.html')
def fastqqc(infile, outfile):
    statement = '''fastqc %(infile)s> %(outfile)s.log'''
    P.run(statement)










#main - will allow to run from cgat core - always goes to the end of the file

if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )