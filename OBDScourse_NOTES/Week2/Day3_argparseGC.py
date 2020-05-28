# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Day 3 in OBDS course
#This will be just for notes

#Supply input FASTA and output TXT file names on the command line
#Open GC program from yesrerday 
#https://docs.python.org/2/howto/argparse.html

#argparse is a way of running scripts from the command line
#will be able to give the file name and script in the command line
#so a way of creating a pipeline etc - define input and output and just feed in files
#to run this script from command line: python <pathtoscript> -f <pathtofile>
#Need to make sure you'e in python 3 environment otherwise it wont work

#Import argparse
import argparse
#define a parser that takes an input and output file
#First defining parser - initialising
#next need to add a description
parser = argparse.ArgumentParser(description='GC content')
#add the argument - for bar is name/flag then need type - so string as it's just the name of your file. help is what you will add, ie inputfile
parser.add_argument('-f', '--file', type=str, help ='inputfile')
args =parser.parse_args()
print(args.file)


#open test file using with open
#We want to count GC content and creat an algorith
#For every character, we want to count it 
#We want to work out it's a G or a C
#If it's a GC we count it as a second feature
#Then we do GC count/total count
gc = 0
len_seq = 0
#here we are pointing to a text - we've now defined the file as args.file so we can run from command line
with open(args.file, "r") as f:
    for line in f:
        for character in line:
            if character == "A" or character == "G" or character == "C" or character == "T":
                len_seq += 1
                if character == "G" or character == "C":
                    gc += 1
#These print statements need to be after the end of the for loops so they print out after it's finished
print(gc & len_seq)
#Calculate the percentage, a // gives you 0 as an integer         
print(gc/len_seq*100)

gc_frac =float(gc/len_se
result =gc_frac*100
print(result)

#Copy the SAM file from shared/week2/ERRetc

#Now write script to convert the SAM file to a BED file
import pysam
#Define input file #rb is used when the input file is not text
infile = pysam.AlignmentFile("-", "rb")
#Define output file to write to "w" 
outfile = pysam.AlignmentFile("-", "w", template=infile)
for s in infile:
    outfile.write(s)
    
#Or try and do as in command line linux: samtools view -S -b example.sam > example.bam
pysam.view("-S", "-b" "input", "output")

#Copy the SAM file from shared/week2/ERRetc

#Now write script to convert the SAM file to a BED file
#comparing SAM file format to BED file format
#Now need to look within the file and if it contains the @ symbol then continue to next lie
#For each line of my sam file, read it and convert to output line (use for loop and return)
#Going to index the file and then rewrite it to a new file by specifying 
#Looking at bed file, want chr then the chr start and end. First Chr in column
#Then the start need to be in -1 
#To work out the end of that region for the bed file, need to do a calculation
#pnext and pos - field number 4 and field number . Problem this is paired end reads so two lines for every
#WE dont know if we need pnext or pos - we just need the start so as long as we take the lowest one
#would need to add tlen to whatever start is to get the stop
#problem is that we have paired reads so a forward and a reverse so p next and pos are sometimes start and end and vice versa
#Now want to separate by tab to convert to BED
#so start of length of 1st fragment to end in next one:

    
#Now want to write a file - create this bed file (line with open "w")
#This will all print to one line unless we include \n at the end
#The fstrings uses the file.write function and f-strings - the f' is essential at start

#Next in the new bed file, there are * for the unmapped regions - inserted if chr == "*" then contiue
#How to remove the duplicates, in the sam file, each pair of reads has an id
#had to define the id of each line, thats the first column (0) 
#We then said if the id of this line is the same as id of last, then skip it (cintinue)
#The splitline #we had to define "\t" to split into each line - otherwise the sam is just continuous






with open("/Users/rhodgson/ERR1755082.test.bed", "w") as output:
    with open("/Users/rhodgson/ERR1755082.test.sam", "r") as sam:
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



##Next slide




#Correct code:
with open("/Users/rhodgson/ERR1755082.test.bed", "w") as output: #This line says this will be the output file to write to
    with open("/Users/rhodgson/ERR1755082.test.sam", "r") as sam: #This is the file to read from
        lastid = -1 #Defines the lastid thing for removing duplicated lines
        for line in sam: #Start the for loop
            if line.startswith('@'): #ignore the lines that start with @, they are headers
                continue #skip
            splitline= (line.split("\t")) #variable splitline splits the files into lines when a space is inserted 
            id = splitline[0] #Assigns id as being the value in the first column of the line
            if lastid == id: #Says that if the lastid (defined above) is the same as id, skip
                continue
            chr = str(splitline[2])  #defines chromosome column as being the 3rd column
            if chr == "*": #syas if chromosome value is * then ignore it, this is because these reads are unmapped
                continue
            start = int(splitline[3])-1 #Defining the start as rname - problem as this is paired ends so may be the reverse stand which will read backwars
            start_mate = int(splitline[7])-1 #Define rnext
            if start < start_mate: #If the start mate is bigger then take that as the start next
                start_bed = start
            else: 
                start_bed = start_mate #otherwise take the first one (which ever value is smaller basically will become the start)
            stop_bed = start_bed+ abs(int(splitline[8])) #Define the stop as being the startbed+ the tlen (length of the temp)
            
            print(chr, start_bed, stop_bed)
            output.write(f'''{chr}\t{start_bed}\t{stop_bed}\n''') #write this line to file

            lastid = id #resets the id to the last id so the loop can continue
        
#Download 2 BED files from ENCODE
#Different tracks for the same cell line
#Same track for different cell lines
#Write a Python script to count the number of intervals that overlap
#between the two  les
#Both  les should be provided as command line arguments
#Advanced: Calculate the number of overlapping bases

#Setting up the arg code: 
#    https://docs.python.org/3/library/argparse.html#dest
#dest 
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

#Then change the file path to args.sam and args.bed



