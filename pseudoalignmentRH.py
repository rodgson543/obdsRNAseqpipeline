#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 11:00:36 2020

@author: rhodgson
"""
#Import everything we need

from ruffus import *
from cgatcore import pipeline as P
import sys



#Import parameters - sure this will change
Params = P.get_parameters ("pipeline_rna_seq.yml")

@transform('*.fastq.gz', suffix('.fastq.gz'),'_fastqc.html')
def fastqc (infile , outfile):
    statement = '''fastqc %(infile)s > %(outfile)s.log'''
    P.run(statement)     