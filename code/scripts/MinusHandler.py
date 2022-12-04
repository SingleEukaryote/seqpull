from Bio.Seq import Seq
# Used in one line for reverse complementing minus sequences
import os
import sys
from os.path import exists
import time
startTime = time.time()
'''
This script reverse complments sequences that blasted against a subject in the minus direction.
The output file is a fasta file with all plus sequences in the blast result, and all reverse complemented minus sequences (tagged in seq name too).
REQUIRES a conda environment with biopython

run:

python minus_handler.py blastresult.tsv sequencesfile.fas outputfilename.fas

* If a sequence's blast results include both plus and minus, it'll output both a normal and a reverse complemented sequence. 
'''

########## INPUT ###########
blastresult = sys.argv[1]
############   blastresults.tsv (output format: -outfmt '6 qseqid sstrand')
############   qseqid must the be first column and sstrand must be the last.
sequencesfile = sys.argv[2]
#############   sequencesfile.fas (fasta file you queried against the database)
outputfilename = sys.argv[3]
##############   outputfilename.fas (the name of the file that contains the plus contigs and reverse complmented minus contigs)
############################

## Log
fastaname = sequencesfile.split('/')[-1]
listname = blastresult.split('/')[-1]
print('\n============================================== MinusHandler.py ==============================================')
print(f'> Reverse complementing sequences minus sequences from {fastaname} \n> Using these blast results: {listname}\n')
## 


### Blast results check



###
#################################################################################################################
###
###### Parse blast results tsv for minus and plus subject sequence names and put them into separate files #######
###
#################################################################################################################

blastresults = open(f'{blastresult}', 'r') # Read blastresults.tsv
minseqnames = open(f'{blastresult}.minus', 'w') # Write temporary file containing sequence names that blasted minus
plusseqnames = open(f'{blastresult}.plus', 'w') # Write temporary file containing sequence names that blasted plus

extractcount = []
linenum = 0
for line in blastresults:
    if line.endswith('minus\n'):
        seqname = line.split('\t')[0] + "\n"
        # Takes the first column (qseqid) string and writes it to a file.
        minseqnames.write(seqname)
        extractcount.append(seqname)
        linenum +=1
    elif line.endswith('plus\n'):
        seqname = line.split('\t')[0] + "\n"
        plusseqnames.write(seqname)
        extractcount.append(seqname)
        linenum +=1
    else:
        linenum +=1
        print("\nYOUR BLAST TABLE IS WRONG\n Make sure output format: -outfmt '6 qseqid sstrand'\n Remember, sstrand must be at the end and qseqid must be at the beginning. \n Error at line#: " +  str(linenum))

plusseqnames.close()
minseqnames.close()

extcountset = set(extractcount)
extcount = len(extcountset)
######################################
###
###### Extract minus sequences #######
###
######################################

# A script for extracting certain sequences from within a FASTA file.
 
# First, convert FASTA file into file with one line per sequence.
# Make sure the name of your FASTA file doesn't contain any dots 
# besides the one before the extension!

contig_list = f'{blastresult}.minus'
output_file = f'{blastresult}.minusseqs'

linenum = 0
onelinefile = f'{blastresult}.minus.tmp'
with open(onelinefile, 'w') as outfile:
    for line in open(sequencesfile, 'r'):
        line = line.strip()
        if len(line) > 0:
            if line[0] == '>':
                if linenum == 0:
                    outfile.write(line + '\t')
                    linenum += 1
                else:
                    outfile.write('\n' + line + '\t')
            else:
                outfile.write(line)
 
# Populate a dictionary using the 'allononeline' file
all_seqs = {}
with open(onelinefile, 'r') as allseqsfile:
    for line in allseqsfile:
        line = line.strip().split('\t')
        all_seqs[line[0][1:]] = line[1]
 
# Generate a set of the sequences you wish to retrieve
desired_seqs = set()
with open(contig_list, 'r') as listfile:
    for line in listfile:
        #print('heyheyhey ' + line)
        desired_seqs.add(line.strip())

# Find the overlap between the total sequences and the desired ones
all_seqs_names = set(all_seqs.keys())
#print(desired_seqs)
#print(all_seqs_names)
toextract = all_seqs_names.intersection(desired_seqs)
#print(toextract)
# Use 'toextract' set to generate desired file
with open(output_file, 'w') as extractfile:
    for name in toextract:
        extractfile.write('>' + name + '\n' + all_seqs[name] + '\n')
 
if exists(onelinefile):
    os.remove(onelinefile)
else:
    print('File empty')

#################################################
###
###### Reverse complement minus sequences #######
###
#################################################
print('Reverse complemented sequences (subject strand minus):')
print('---------------------------------------------------------')

minusseqs = f'{blastresult}.minusseqs'
revcomped = f'{blastresult}.minusseqs.revcomp'

minussequences = open(f'{minusseqs}', 'r')
revcompedsequences = open(f'{revcomped}', 'w')
for line in minussequences:
    if line.startswith('>'):
        sequence = minussequences.readline()
        revcomp = Seq(sequence).reverse_complement()
        header = line.split('\n')[0] + '.revcomp'
        revcompedsequences.write(header)
        revcompedsequences.write(str(revcomp) + '\n')
        print(header)
        if sequence.split('\n')[0] == revcomp:
            print('Reverse complement failed for: ' + header)
        else:
            pass
        #print('Original sequence:')
        #print(sequence + '\n')
        #print('Reverse complemented sequence:')
        #print(revcomp + '\n\n')
minussequences.close()
revcompedsequences.close()
print('---------------------------------------------------------\n')
#####################################
###
###### Extract plus sequences #######
###
#####################################
print('Sequences kept the same (subject strand plus):')
print('---------------------------------------------------------')

contig_list = f'{blastresult}.plus'
output_file = f'{blastresult}.plusseqs'

linenum = 0
onelinefile = f'{blastresult}.plus.tmp'
with open(onelinefile, 'w') as outfile:
    for line in open(sequencesfile, 'r'):
        line = line.strip()
        if len(line) > 0:
            if line[0] == '>':
                if linenum == 0:
                    outfile.write(line + '\t')
                    linenum += 1
                else:
                    outfile.write('\n' + line + '\t')
            else:
                outfile.write(line)
 
# Populate a dictionary using the 'allononeline' file
all_seqs = {}
with open(onelinefile, 'r') as allseqsfile:
    for line in allseqsfile:
        line = line.strip().split('\t')
        all_seqs[line[0][1:]] = line[1]
 
# Generate a set of the sequences you wish to retrieve
desired_seqs = set()
with open(contig_list, 'r') as listfile:
    for line in listfile:
        desired_seqs.add(line.strip())
 
# Find the overlap between the total sequences and the desired ones
all_seqs_names = set(all_seqs.keys())
toextract = all_seqs_names.intersection(desired_seqs)

# Use 'toextract' set to generate desired file
with open(output_file, 'w') as extractfile:
    for name in toextract:
        extractfile.write('>' + name + '\n' + all_seqs[name] + '\n')
 
if exists(onelinefile):
    os.remove(onelinefile)
else:
    print('File empty')

with open(output_file, 'r') as seqs:
    for line in seqs:
        if line.startswith('>'):
            print(line.split('\n')[0])
print('---------------------------------------------------------')

#########################################
###
###### Write output sequence file #######
###
#########################################

os.system(f'cat {blastresult}.plusseqs > {outputfilename}')
os.system(f'cat {blastresult}.minusseqs.revcomp >> {outputfilename}')

filesize = os.path.getsize(outputfilename)
if filesize == 0:
    print("\n\n======================== MINUS HANDLER ================================\n")
    print(f"{outputfilename} is empty: " + str(filesize) + "\n Deleting empty file.")
    print("\n=======================================================================\n\n")
    #os.remove(outputfilename)
else:
    pass

#### LOGS

# Count sequences
## Input file
output = open(f'{sequencesfile}', "r")
data = output.read()
originalseqs = data.count('>')
output.close()
## Plus seqs
output = open(f'{blastresult}.plusseqs', "r")
data = output.read()
plusseqs = data.count('>')
output.close()
## Reverse complemented seqs
output = open(f'{blastresult}.minusseqs.revcomp', "r")
data = output.read()
minusseqs = data.count('>')
output.close()
## Amount of seqs in final file
output = open(f'{outputfilename}', "r")
data = output.read()
realnumseqs = data.count('>')
output.close()


print('Sequences extracted from blast results:       ' + str(extcount))
print('Sequences only kept the same:                 ' + str(plusseqs-(realnumseqs-extcount)))
print('Sequences only reverse complemented:          ' + str(minusseqs-(realnumseqs-extcount)))
print('Sequences both kept and reverse complemented: '+ str(realnumseqs-extcount) + ' times two = ' + str(2*(realnumseqs-extcount)))
print('Final number of sequences:                    ' + str(realnumseqs))

########################################
###
###### Remove intermediate files #######
###            (optional)
########################################


#os.remove(f'{blastresult}.plus')
#os.remove(f'{blastresult}.plusseqs')
#os.remove(f'{blastresult}.minus')
#os.remove(f'{blastresult}.minusseqs')
#os.remove(f'{blastresult}.minusseqs.revcomp')


executionTime = (time.time() - startTime)
print('Script run time:                              ' + str(executionTime) + ' seconds')
print('============================================== MinusHandler.py ==============================================\n')
