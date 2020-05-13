#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 11:00:36 2020

@author: rhodgson
"""
#To run this:
    #python <script name> make -v 5

#Import everything we need
#Working directory: /ifs/obds-training/apr20/rose/pipelines/pseudoalignment
#Make a directory with where we want the files- made a symbolic link to fastq files from/ifs/obds-training/apr20/exercises/rnaseq/
#will be rsyncing script:
    #rsync -a /Users/rhodgson/GitHub/obdsRNAseqpipeline/pseudoalignmentRH* rose@cgatui.imm.ox.ac.uk:/ifs/obds-training/apr20/rose/pipelines/pseudoalignment
#We will be using kallisto

from ruffus import *
from cgatcore import pipeline as P
import sys

#So the first thing we still want to do is the QC of the fastq files (including multiQC)
#see file pipeline_rna_seq.py for notes on these functions



#Import parameters - sure this will change
Params = P.get_parameters ("pseudoalignmentRH.yml")

#Fastqc - just added an output folder here
#as this is a transform thing, we needed to use a regular expression to put these reports into a new folder

@follows(mkdir('fastqc_reports'))
@transform('*.fastq.gz', regex(r'(.*).fastq.gz'),r'fastqc_reports/\1_fastqc.html')
def fastqc (infile , outfile):
    statement = '''fastqc --outdir fastqc_reports  %(infile)s > %(outfile)s.log'''
    P.run(statement)     

#Multiqc -hmm
@follows(mkdir('multiqc_reports'))
@merge(fastqc, 'multiqc_reports/multiqc.html')
def multiqc (infiles , outfile):
    statement = ''' export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 
    && multiqc . -f -n %(outfile)s '''
    P.run(statement)  
    
#Next we will do the pseudo alignment, ie the kallisto
#Build a Kallisto index - first look at kallisto help
#kallisto index wants a string (-i filename) FASTA files
#We want the fasta files of all the transcripts of the genome - go to ensembl and download a file
#ensembl > mouse > cdna> README gives info on these files
#wget ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz

#This step is literally taking the fasta genome file and converting it to a debruijn graph
#At this point, we are not using the fastq files
@transform(Params["kallisto_index"], regex(r'(.*).fa.gz'), 'transcripts.idx')
def kallisto_index(infile, outfile):
    statement = ''' kallisto index -i %(outfile)s %(infile)s'''
    P.run(statement)  
    
#Next step is to do kallisto quant - transcript.idx 
#use quant to get a tsv file
#We are inputting the fastq files to do the quantification
#to put the output of each pair of fastq files into a newfolder: 
#use regex - name folder with \1 as the name of the sample
    #using the r protects the backslashes
#need to sort out the output_dir - says find this bit, and replace with an empty string

#We want to capture kallisto output (from normally stdout to a log file)
#> redirects stdout, 2 > redirects stderr, < redirects stdin

@follows(mkdir('kallisto_quant_output'),kallisto_index)
@collate('*.fastq.gz', regex(r'(.+)_[12].fastq.gz$'), r'kallisto_quant_output/\1/abundance.tsv')
def kallisto_quant(infiles, outfile):
    input_list = " ".join(infiles)
    output_dir = outfile.replace('/abundance.tsv' , '' )
    statement = '''kallisto quant  -t %(kallisto_threads)s %(kallisto_options)s 
    -i transcripts.idx -o %(output_dir)s %(input_list)s > %(output_dir)s/stdout.log 
    2> %(output_dir)s/stderr.log'''
    P.run(statement, job_threads = Params["kallisto_threads"])

#Now want to do multiqc after
#tells you how many reads it pseduo aligned and how many it cant 
@follows(mkdir('multiqc_reports'))
@merge(kallisto_quant, 'multiqc_reports/multiqc_kallisto.html')
def multiqc_kallisto (infiles , outfile):
    statement = ''' export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 
    && multiqc kallisto_quant_output/ -f -n %(outfile)s '''
    P.run(statement)  
 
#The main bit to let it run from cgat core
if __name__ == "__main__":
    sys.exit( P.main (sys.argv))