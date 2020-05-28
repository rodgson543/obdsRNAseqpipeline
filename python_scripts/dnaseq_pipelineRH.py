#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
====================
DNAseq pipeline
====================
   1. Quality control reads with fastqc and multiqc
   2. Align to genome using gapped alignment (BWA-MEM)
   3. Check alignment quality and target region coverage (Picard)
   4. Call varinats with Strelka2
   5. Annotate variants with Ensembl VEP

"""
#strelka2 - mutation caller. bench

import sys
from ruffus import *
from cgatcore import pipeline as P

params = P.get_parameters("dnaseq_pipeline.yml")

#fastq input
#Do fastqc of the the files and put them in a new directory. name file using the first part of the regex
#no group is an option for fastqc - where you have plots, it doesn't average over bases. give you bases over evry base. 
#params["q"] in param file - put the name of the q to allow you to put the pipeline from one to another q


@follows(mkdir("fastqc"))
@transform("*.fastq.gz", regex(r"(.*).fastq.gz"), r"fastqc/\1_fastqc.html")
def fastqc(infile, outfile):
    statement = "fastqc --nogroup -o fastqc %(infile)s > %(outfile)s.log"
    P.run(statement,job_queue = params["q"])

#MultiQC 
@follows(mkdir("reports"))
@merge(fastqc, "reports/fastqc_report.html")
def multiqc(infiles, outfile):
    statement = """export LC_ALL=en_US.UTF-8; export LANG=en_US.UTF-8;
                            multiqc fastqc/ -f -n %(outfile)s"""
    P.run(statement,job_queue = params["q"])

#BWA #need to install this in conda environmentn 
#Mapping reads - index is the genome reference to align to
#paired end read - tuple, collate gives a tuple as a list - need to supply these separately on CL, here's read one and read 2 
#pipe to samtools sort to index - will generate a bam per pair of fastq files
#Have to give threads to both the program and to SGE in run statemtn

@follows(mkdir("bam"))
@collate("*.fastq.gz", regex(r"(.+)_[12].fastq.gz"), r"bam/\1.bam")
def bwa(infiles, outfile):
    read1, read2 = infiles
    statement = """bwa mem
                           -t %(bwa_threads)s
                           %(bwa_index)s
                           %(read1)s %(read2)s
                           2> %(outfile)s.log |
                           samtools sort -o %(outfile)s -@ %(samtools_threads)s -
                           && samtools index %(outfile)s"""
    P.run(statement,job_queue = params["q"],
          job_threads = params["bwa_threads"],
          job_memory = params["bwa_memory"])


#Take bwa bams - transform them to remove duplicates 
#Picard is dodgy - needs a buffer for final memory 

@transform (bwa, regex(r'(.*).bam'), r'\1.rmdup.bam')
def rmdup (infile, outfile):
     final_memory = str(int(params['picard_memory'])+ 2)+'g'
     statement = """picard -Xmx%(picard_memory)sg
                                MarkDuplicates
                                I= %(infile)s O=%(outfile)s M=%(outfile)s.metrics
                                && samtools index %(outfile)s"""
     P.run(statement, job_queue = params['q'], job_memory= final_memory)

#This step is look at the alignment summary metrics
#R= reference gneome
#Gives you a metircs file

@transform (bwa, regex(r'(.*).bam'), r'\1.picardmetrics')
def picard (infile, outfile):
     final_memory = str(int(params['picard_memory'])+ 2)+'g'
     statement = """picard -Xmx%(picard_memory)sg
                                CollectAlignmentSummaryMetrics
                                R=%(picard_genome)s I= %(infile)s O=%(outfile)s"""
     P.run(statement, job_queue = params['q'], job_memory= final_memory)

#IDXstats
@transform (bwa, regex(r'(.*).bam'), r'\1.idxstats')
def idxstats (infile, outfile):
     statement = """samtools idxstats %(infile)s > %(outfile)s"""
     P.run(statement, job_queue = params['q'])

#flagstats - quality check
@transform (bwa, regex(r'(.*).bam'), r'\1.flagstats')
def flagstats (infile, outfile):
     statement = """samtools flagstats %(infile)s > %(outfile)s"""
     P.run(statement, job_queue = params['q'])

#merge the output of the aligned bwa files - QC
@merge ([picard, idxstats, flagstats, rmdup], 'reports/bamqc_report.html')
def multiqc_bam (infiles, outfile):
    statement= """export LC_ALL=en_US.UTF-8;export LANG=en_US.UTF-8;
                            multiqc bam/ -n %(outfile)s -f"""
    P.run (statement, job_queue = params["q"])


#Take marked duplicates file - has it removed the duplicates?!
#strelka - joint allele caller variant caller
#From internet: You should always have PCR duplicates either marked or filtered before running Strelka. If PCR duplicates are marked, these reads will be ignored by the method.
#Strelka - uses python2 so in P run function - new conda environment - pls run in strelka environment (but only this job!)
#-options: -m ? local refers to wherever job is submitted
#conda create -n strelka strelka (environemnt based on one package)#conda will sort this out and wrap it 

@follows(mkdir("strelka"))
@transform (rmdup, regex(r'bam/(.*).rmdup.bam'), r'strelka/\1/results/variants/variants.vcf.gz')
def strelka (infile, outfile):
     outdir = outfile.replace("/results/variants/variants.vcf.gz","")
     statement = """configureStrelkaGermlineWorkflow.py
                                --bam %(infile)s
                                --referenceFasta %(picard_genome)s
                                --runDir %(outdir)s  > %(outfile)s.log &&
                                %(outdir)s/runWorkflow.py -m local -j 12"""
     P.run(statement, job_queue = params['q'], job_condaenv="strelka", job_threads=12)

#VEP - variant annotation tool using vep
#transform strelka to captue vcf.gz files, output will be a referral back to parantheses are coming back
#conda create -n vep ensembl-vep #do all this in vep environment
#Do we want to zip this? --vcf makes it go to vcf - seems like it and index theoutput fille with tabix

@follows(mkdir("vep"))
@transform (strelka, regex(r'strelka/(.*)/results/variants/variants.vcf.gz'), r'vep/\1.vcf.gz')
def vep (infile, outfile):
     statement = """vep --cache --dir %(vep_dir)s --vcf --compress output gzip
                            -i %(infile)s  -o %(outfile)s  > %(outfile)s.log 
                            && tabix -p vcf %(outfile)s"""
     P.run(statement, job_queue = params['q'], job_condaenv="vep")
     
#Next function we want to use transform to filter the variant filtering from vep based on depth
#have to use bcftools, FORMAT/DP is a column  
#-O gives the type of output 
#had to create a conda environment wiith bcftools

@transform(vep, regex(r'vep/(.*).vcf.gz'), r'vep/\1.filtered.vcf.gz')
def bcftools(infile, outfile):
    statement = """bcftools filter -i 'FORMAT/DP[0]>10' -O z%(infile)s > %(outfile)s """
    P.run(statement, job_queue = params['q'], job_condaenv="bcftools")
    
    
#filter based on region - provided this in a file using -R. Got file off here and created it as a bed file by specifying 
#https://m.ensembl.org/biomart/martview/29e2740001b97f8ea744104c6a53f1c8
@transform (bcftools, regex(r'vep/(.*).filtered.vcf.gz'), r'vep/\1.filteredregion.vcf.gz')
def bcftools_region (infile, outfile):
    statement = """bcftools filter -R %(bcftools_region)s -O z %(infile)s > %(outfile)s"""
    P.run(statement, job_queue = params['q'], job_condaenv="bcftools")    
    
#Run on the chat server
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
