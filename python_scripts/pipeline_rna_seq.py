#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 10:16:38 2020

@author: rose
"""
#RNAseq pipeline
#Whole point of this is to enable you to run all of these commmands one after another using ruffus

#rsync -a /Users/rhodgson/GitHub/OBDS_Training_Apr_2020/pipeline_rna_seq.py rose@cgatui.imm.ox.ac.uk:/ifs/obds-training/apr20/rose/pipelines/rnaseqpipeline/

#import packages 
from ruffus import *
from cgatcore import pipeline as P
import sys



#Import parameters: contains info for probelms like hisat, you would adjust this file for different data, ie diff genome
#Make sure this yml file is in the same folder as script and data
Params = P.get_parameters ("pipeline_rna_seq.yml")

#This part is the fastq to generate fastqc html files
#Input the fastq.qz anything ending in this, defining wiht suffix. 
#Then tell it output being fastqc.html etc
#Decorator transform changes one files to another, need to define infile/outfile: variable names
# ruffus statement is what the command line takes, variables have to be in %()s - 
#def = define function 
#P.run takes the statement - ad local stuff, makes statement and submits it to the cluster
@transform('*.fastq.gz', suffix('.fastq.gz'),'_fastqc.html')
def fastqc (infile , outfile):
    statement = '''fastqc %(infile)s > %(outfile)s.log'''
    P.run(statement)                
    
    
#Next want to do multiqc to make a nice report containing all the fastqc from each fastq file
#We want to merge the input from fastqc files into multiqc report
#First use decorator @follows to make a directory for output
#@follows says dont do antyhing until you've done this, ie dont do merge until made folder

#use the input from fastqc in the merge function, output will be multiqc.html report
#define multiqc function, need multiple infiles - get these from fastqc function
#run this like  . to look in cd for all fastqc files in that directory
#Name the output file usin -n and specify output directory using -o 
#the -f overwrites previous runs from other directories
#Would have to look at the multiqc file on our own comp using rysnc 

@follows(mkdir('multiqc_reports'))
@merge(fastqc, 'multiqc_reports/multiqc.html')
def multiqc (infiles , outfile):
    statement = '''export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && 
    multiqc . -f -n %(outfile)s'''
    P.run(statement)  


#Next function is mapping to genome - Hisat
#collate - use this to give 2 input files and 1 output file using a regex
#Use a regex to allow input for two sets of file, ie for paired end sequence mapping
#The regex 
#To define the output name: \1 takes matched group between () as same, ie refers back to that group of files, reuses it for name of the bam
#tuple here fastq1 and fastq2 to define infiles
#parameterise using a parameters yml file that we supply in the same folder  
#define this file in import parameters file above (yml file)
#statement gives you the command line use thing
#we piped the sam file to a sorted bam file to cut out steps, and indexed it to get the bai file
#the dash after samtools - will act as stdout and pipe it through
#&& will do "only if no previous errors from previous command, then do this"
@follows(mkdir('bamfiles'))
@collate('*.fastq.gz', regex(r'(.+)_[12].fastq.gz$'),r'bamfiles/\1.bam')
def hisat (infiles , outfile):
    fastq1,fastq2 = infiles
    statement = '''hisat2 %(hisat_options)s --threads %(hisat_threads)s -x %(hisat_genome)s 
    -1 %(fastq1)s -2 %(fastq2)s --summary-file %(outfile)s.txt | 
    samtools sort - -o %(outfile)s 
    && samtools index %(outfile)s '''
    P.run(statement, job_threads = Params["hisat_threads"], job_queue = "all.q")    


#Run in the output files from hisat to transform into idxstats thing 
#this will send files to the bamfiles folder as we prviously specificied this in hisat command above
#idxstats: reports alignment summary stats to a text file
    
@transform(hisat, suffix('.bam'), '.idxstats')
def idxstats(infile, outfile):
    statement = '''samtools idxstats %(infile)s > %(outfile)s'''
    P.run(statement)

#Next bit will be to do the flag stat thing 
#The input file will be 

@transform(hisat, suffix('.bam'), '.flagstat')
def flagstat(infile, outfile):
    statement = '''samtools flagstat %(infile)s > %(outfile)s'''
    P.run(statement)

#Then to do feature counts - use merge to go from all the bam files to one count table
#how to put in input files as a list, takes list, create a string where each thing is 
#Again use the yml parameters file to give some options from the yml file
#-T is threads, annotation and options 
#in p.run statement we want to give it 12 threads from cgat core rather than only requesting 1 
@merge(hisat, 'count.table')
def featureCounts(infiles, outfile):
    input_list = " ".join(infiles)
    statement = ''' featureCounts %(featureCounts_options)s -T %(featureCounts_threads)s -a %(featureCounts_annotation)s 
    -o %(outfile)s %(input_list)s'''
    P.run(statement, job_threads = Params["featureCounts_threads"], job_queue = "all.q")


#The main bit to let it run from cgat core/server. Gives control to pipeline software
#when you run this script, run the main argument from cgat core - pass it all of the arguments and run it all

if __name__ == "__main__":
    sys.exit( P.main (sys.argv))