#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ChIP-seq pipeline
"""
#copy files from shared to my directory. They are already symbolic links so have to copy
#cp -d /ifs/obds-training/apr20/shared/week3/chipseq/* .
#will be rsyncing script
#rsync -a /Users/rhodgson/GitHub/OBDS_Training_Apr_2020/chipseq_pipeline.* rose@cgatui.imm.ox.ac.uk:/ifs/obds-training/apr20/rose/pipelines/chipseqpipeline

import sys
from ruffus import *
from cgatcore import pipeline as P

params = P.get_parameters("chipseq_pipeline.yml")

#fastq input
#Do fastqc of the the files and put them in a new directory. name file using the first part of the regex
#no group is an option for fastqc - where you have plots, it doesn't average over bases. give you bases over evry base. 
#params["q"] in param file - put the name of the q to allow you to put the pipeline from one to another q

@follows(mkdir("fastqc"))
@transform("*.fastq.gz", regex(r"(.*).fastq.gz"), r"fastqc/\1_fastqc.html")
def fastqc(infile, outfile):
    statement = "fastqc --nogroup -o fastqc %(infile)s > %(outfile)s.log"
    P.run(statement,job_queue = params["q"])

#Next bit -multiqc

@follows(mkdir("reports"))
@merge(fastqc, "reports/fastqc_report.html")
def multiqc(infiles, outfile):
    statement = """export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8;
    multiqc fastqc/ -f -n %(outfile)s"""
    P.run(statement,job_queue = params["q"])

#bam - collate - input from fastq, regex to take files of the same prefix with either 1 or 2
#output using regex, to put files in bam folder and gice them a specific name based on input and regex
#bowtie, specific options. 
#index is genome (params)
#stderr pieped to output and pipe stdout to samtools and then index the file
#-@ means threads in samtools
#
@follows(mkdir("bam"))
@collate("*.fastq.gz", regex(r"(.+)_[12].fastq.gz"), r"bam/\1.bam")
def bowtie2(infiles, outfile):
    read1, read2 = infiles
    statement = """bowtie2 %(bowtie2_options)s 
                   -p %(bowtie2_threads)s 
                   -x %(bowtie2_index)s 
                   -1 %(read1)s -2 %(read2)s 
                   2> %(outfile)s.log | 
                   samtools sort -o %(outfile)s -@ %(bowtie2_threads)s - 
                   && samtools index %(outfile)s"""
    P.run(statement,job_queue = params["q"],
          job_threads = params["bowtie2_threads"],
          job_memory = params["bowtie2_memory"])
    
#input from bowtie 2 
#regex gets files ending in bam. output will end in picardmetrics and starts with filename. 
#picard manual. runs using java - ok because we are usign this in a statement for command line 
#as long as it's on the command line - ie java installed in environment. normally we have to provide java jar before the tool
#within obds-py3 - conda provides a convinience wrapper, so when you type in picard, conda behind the scenes adds java jar so you can just run picard

#alignment summary metrics is the function we want. - option is collectalignmentsummarymetrics
#R I O #ref - fasta#input#output - txt file
#overheads of memory as an option here - in params we're adding 2gs onto final memory as picard is memory intensive


@transform (bowtie2, regex(r'(.*).bam'), r'\1.picardmetrics')
def picard (infile, outfile):
     final_memory = str(int(params['picard_memory'])+ 2)+'g'
     statement = """picard -Xmx%(picard_memory)sg CollectAlignmentSummaryMetrics 
                    R=%(picard_genome)s I= %(infile)s O=%(outfile)s"""
     P.run(statement, job_queue = params['q'], job_memory= final_memory)


#Doing idxstats on bam files (mapped files right?)
@transform (bowtie2, regex(r'(.*).bam'), r'\1.idxstats')
def idxstats (infile, outfile):
     statement = """samtools idxstats %(infile)s > %(outfile)s"""
     P.run(statement, job_queue = params['q'])    
     
#flagstats - looking at reads per flag
@transform (bowtie2, regex(r'(.*).bam'), r'\1.flagstats')
def flagstats (infile, outfile):
     statement = """samtools flagstats %(infile)s > %(outfile)s"""
     P.run(statement, job_queue = params['q'])    

#Adding step to remove duplicates
#input output logfile
#removed duplicates here using picard
#also had to index the bam file
@transform (bowtie2, regex(r'(.*).bam'), r'\1.rmdup.bam')
def picard_rmdup (infile, outfile):
     final_memory = str(int(params['picard_memory'])+ 2)+'g'
     statement = """picard -Xmx%(picard_memory)sg MarkDuplicates REMOVE_DUPLICATES = true
                    I= %(infile)s O=%(outfile)s M= %(outfile)s.metrics
                    && samtools index %(outfile)s"""
     P.run(statement, job_queue = params['q'], job_memory= final_memory)

#We're now going to filter unmapped reads using samtools 
#Look at flags and think abotu values for the filters flags: http://broadinstitute.github.io/picard/explain-flags.html
#-f 1 keeps the paired reads, -f 2 keeps reads mapped in a proper pair
#-F gets rid of unmapped pairs
#-F 0x100 gets rid of not primary alignment
#-b makes it a bam file, -h keeps the header
#Dave Tang page to explain https://davetang.org/wiki/tiki-index.php?page=SAMTools
#Index at the end
#So so far at this step - we want to filter the output of picard removed duplicates so input will be picardrmdup
#sort out regex so the rmdup will be removed 
#To filter out mitochondrial reads too, use after input file chr{1..19} - this will select only those files
@transform(picard_rmdup, regex(r'(.*).rmdup.bam'), r'\1.filter_read.bam')
def filter_read(infile, outfile):
    statement = """samtools view -h -b -f1 -f2 -F 4 -F 0x100 %(infile)s chr{1..19} > %(outfile)s
    && samtools index %(outfile)s"""
    P.run(statement,job_queue = params['q'])


     
#Now to filter blacklisted peaks/regions
#This function will remove the regions from the params file using bedtools, also index to 
#We give the file in the yml parameters file to filter on. abam says write a bam
#
@transform (filter_read, regex(r'(.*).filter_read.bam'), r'\1.whitelist.bam')
def blacklist_filter(infile, outfile):
    statement = """ bedtools intersect -v -abam %(infile)s -b %(blacklist)s >
    %(outfile)s && samtools index %(outfile)s """
    P.run(statement, job_queue = params['q'])

#merging these QC params into a multuQC file - givinga list of functions that it's merging 
#will take stats from each of these 
#reports folder is already made - but chances are mapping will take an age so we know the first multiqc and folder creation
@merge ([picard, idxstats, flagstats], 'reports/bamqc_report.html') 
def multiqc_bam (infiles, outfile):
    statement= """export LC_ALL=en_US.UTF-8;export LANG=en_US.UTF-8; 
                  multiqc bam/ -n %(outfile)s"""
    P.run (statement, job_queue = params["q"])
    
    
#macs2 does the peak calling: has 2 modes, compares chip sample to an input. OR can look at treatment alone, looks at background to decide whether something is a peak or not
#so you can either give just treatment, or treatment and control. 
#regex ip is input and gr is glucoR - not very generic pipeline and based on naming convention of samples. We should generalise it
#You could give a txt file that contains a list of your files if your naming system is more complex
#to name - replacing peaks narrows - pulling out just the \1 to recover name of sample and name the output.
#tuple to define input
#callpeak is function - give it the treatment and the control.   
#format is BAM - probably should have put in BEDPE - as this may only take the first read
#genome in yml
#redircts to the outdir peaks

@follows (mkdir ('peaks'))
@collate(blacklist_filter, regex(r"bam/(.*)_(ip|gr).bam"), r"peaks/\1_peaks.narrowPeak")
def macs2(infiles, outfile):
    treatment, control = infiles
    name = outfile.replace('_peaks.narrowPeak', "").replace('peaks/', '')
    statement = """macs2 callpeak -t %(treatment)s -c %(control)s
                   -n %(name)s -f BAM -g %(macs2_genome)s --outdir peaks"""
    P.run(statement,job_queue = params["q"]) 

#Take output of macs2 and use bedtools to merge into a single peak -  
#Creates an allpeaks bed file - this will quantify all peaks in all samples - gives you exisiting peaks in genome 

@merge(macs2,"allpeaks.bed")
def merge_peaks(infiles, outfile):
    inlist = " ".join(infiles)
    statement = """cat %(inlist)s | sort -k1, 1 -k2, 2n | mergeBed -i stdin  > %(outfile)s"""
    P.run(statement, job_queue = params['q'])

#Next we want to count the reads under the peaks using bedtools - report counts in alignments of each inoput that overlap 
#Output will be table of counts
#bring in all bam files, add input - hardcode allpeaks.bed 
#redirect output to outfile - this is the counts 
@merge(blacklist_filter, "peak_counts")
def coverage(infiles, outfile):
    inlist = " ".join(infiles)
    statement = """bedtools multicov -bams %(inlist)s -bed allpeaks.bed > %(outfile)s"""
    P.run(statement,job_queue = params["q"]) 

#Next going to use homer to annotate the peaks
#We can download reference genomes 
#we didnt do this as homer is an absolute nightmate
#Basically homer uses perl - it's not well written and you need to install the reference genomes
#seems like you have to find a path for where homer is installed etc. 
#Would run this basically
#annotatePeaks.pl peaks/rep1_dex_peaks.narrowPeak mm10

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
