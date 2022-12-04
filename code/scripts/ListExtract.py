#! /usr/bin/env bash
#! /usr/bin/env python

'''
Made by Tristan ^^

python ListExtract.py [fasta file] [list of search terms] [output file]

This script:
- Extracts sequences from a FASTA file using a list of search terms.

The list may be super specific (exact sequence names), or super general (the letter 'g').
The more general the search terms are, the longer the script will take.

Ex. The list file below:
--------------------
g
t
s
--------------------
   will extract all sequences with 'g', 't', or 's' letters in the sequence name.

Real ex.
- Extract sequences from a datafile using blast results
- If you blastn with -outfmt '6 sseqid', you'll get a file with only sequence names. Use it and the original fasta file to extract sequences based on blast results

Ignores duplicate lines and outputs 1 copy of duplicate sequences. 

If the terms you wanted to search for resulted in no sequences, then nothing will happen.
(Empty files are deleted)

* Accepts multiple line sequences too: convert FASTA file into file with one line per sequence.

!!! NEW FEATURE
---------------
Don't extract a sequence.
For example if you list 

'proteobacteria'
'(:not:)Gammaproteobacteria'

Then it'll extract all sequences with 'proteobacteria' in the name, unless it also has 'Gammaproteobacteria' in the name.

If

'g'
'(:not:)g'

is used, then it will do nothing.
---------------
Change
'''
#################### IMPORTS
import sys
import os
from os.path import exists
#### To measure the run time
import time
startTime = time.time()
############################

######### INPUTS #########
fasta_file = sys.argv[1]  # File full of sequences
searchlist = sys.argv[2]    # File with search terms in every line (ex. 'potato')
output_file = sys.argv[3] # Desired name and path of the output file 
##########################

## Log
fastaname = fasta_file.split('/')[-1]
listname = searchlist.split('/')[-1]
print('\n========================================= ListExtract.py =========================================')
print(f'> Extracting sequences from: {fastaname}')
print(f'> Using the search list:     {listname}')

## 


#########################################################################################
# Deals with blast result names -outfmt '6 sseqid' and spaces in sequence name headers
# Finds full sequence name based on list of names

## remove duplicate lines (makes the script faster)
lines_seen = set() # holds lines already seen
nots_seen = set()
terms_seen = set()
no_dupes = open(f'{searchlist}_nodupes.temp', "w")


notnumber = 0
notlist = []
for line in open(f'{searchlist}', "r"):
    if line not in lines_seen: # not a duplicate
        no_dupes.write(line)
        lines_seen.add(line)
        terms_seen.add(line.split('\n')[0])
    if line.startswith('(:not:)') and line.split('\n')[0] not in nots_seen: # Searches for sequence terms you do not want extracted
        notthis = line.split('(:not:)')[1].split('\n')[0]
        notlist.append(notthis)
        nots_seen.add(line.split('\n')[0])
        notnumber += 1
no_dupes.close()


fullseqlist = open(f'{searchlist}_fullnames', 'w')

## List all sequence names in searchlist file
## Ex. TRINITY_DN22402_c0_g1_i1

with open(f'{searchlist}_nodupes.temp','r') as f:
    queries = [line.strip() for line in f]

## Query your list against sequence headers in fasta file, and write the matches to a new file
## Ex. TRINITY_DN22402_c0_g1_i1 len=300 path=[0:0-299] will be found
seqsfound = 0

termSet = set(terms_seen)
termDict = dict.fromkeys(termSet, 0) # Get metrics on how many sequences a term found to extraction

notSet = set(nots_seen)
notSetlog = 0
if len(notSet) == 0:         # this prevents the empty set bug below. If there are no nope in notSet, then the sequence list won't be written. I know it isn't elegant...
    notSet.add('(:not:)fuckdootledothisstringshouldnotbeinanysequenceknowntoman')
    notSetlog = 1
notDict = dict.fromkeys(notSet, 0) # Get metrics on how many sequences prevented from extraction


### Meat of the script
redundancies = 0
realseqsnoped = 0
with open(f'{fasta_file}','r') as f:
    for line in f:
        if line.startswith('>'):
            redundant = 0
            for query in queries:
                no = 0
                if query in line:
                    seqsfound += 1
                    termDict[query] += 1
                    redundant += 1
                    if redundant >= 2:
                        redundancies += 1
                    for nope in notSet:
                        notthis = nope.split('(:not:)')[1].split('\n')[0]
                        if notthis in line: # Counts occurences of above
                            no += 1
                            notDict[nope] += 1
                            if redundant == 1:
                                realseqsnoped += 1
                            break
                if query in line and no == 0:
                    nocarrot = line[1:]
                    fullseqlist.write(nocarrot) # Write sequence name that has a search term in it but not a (:not:)term in it

fullseqlist.close()
os.remove(f'{searchlist}_nodupes.temp')
#########################################################################################
##
#time.sleep(2)
##

print(f'> Putting sequences into:    {output_file}')
## Extract sequences based on full sequence name list
linenum = 0
onelinefile = f'{output_file}.tmp'
with open(onelinefile, 'w') as outfile:
    for line in open(fasta_file, 'r'):
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
with open(f'{searchlist}_fullnames', 'r') as listfile:
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
    print('Sequence file is empty')


###############################################################################
### Logs and quality assurance, remove below if you think it's unnecessary or time consuming
### Making these logs was really a fun coding challenge for myself
###############################################################################
### Check which search terms did not bring up results ###
#logTime = time.time()

searchterms = open(f'{searchlist}', 'r')

termlist = []
listcount = 0
for line in searchterms:
    term = line.split('\n')[0]
    if term in notlist:
        pass
    else:
        termlist.append(term)
        listcount += 1

searchterms.close()

termset = set(termlist)

searchfail = []
searchpass = []
terms = 0
passterms = 0
failterms = 0
totalseqs = 0


totalnopes = 0
print('\n----------- extracted sequences -----------')
for term in termset:
    seqcount = 0
    if term not in searchpass and term not in searchfail:
        sequences = open(f'{output_file}', 'r')
        for line in sequences:
            if line.startswith('>') and term in line:
                seqcount += 1
                searchpass.append(term)
            else:
                searchfail.append(term)

    sequences.close()

    passset = set(searchpass)
    failset = set(searchfail)
    if term in passset and term not in notSet:
        print(str(termDict[term]) + " (+) sequences found using '" + term + "'")
        totalseqs += seqcount
        terms += 1
        passterms += 1
    if term in failset and term not in passset and term not in notSet:
        print(str(termDict[term]) + " (+) sequences found using '" + term + "', but they were noped.")
        terms += 1
        failterms += 1
    if term in notSet:
        print(str(notDict[term]) + " (-) successful searches noped using '" + term + "'")
        totalnopes += notDict[term]

output = open(f'{output_file}', "r")
data = output.read()
realnumseqs = data.count('>')
output.close()
print('-------------------------------------------')
print(str(seqsfound) + ' successful searches.')
print(str(redundancies) + " redundant searches.")
print(str(seqsfound-redundancies) + " sequences found in searches.")
if notSetlog == 0:
    print(str(realseqsnoped) + ' sequences noped using (:not:) terms.')
#print(str(totalseqs) + ' total extractions.')
print(str(realnumseqs) + " total sequences extracted.")

print('-------------------------------------------')

print('Total search terms:    ' + str(listcount))
print('Unique search terms:   ' + str(terms))
if notSetlog == 0:
    print('(:not:) search terms:  ' + str(notnumber))
print('Successful searches:   ' + str(passterms))
print('Failed searches:       ' + str(failterms))
executionTime = (time.time() - startTime)
#reallogtime = (time.time() - logTime)
print('Script run time:       ' + str(executionTime) + ' seconds')
#print(str(reallogtime))
filesize = os.path.getsize(output_file)
if filesize == 0:
    print('\n No sequences extracted. Try making a better search list.')
    os.remove(output_file)
    os.remove(f'{searchlist}_fullnames')                            
                            
    print('  ___  ___  _ __ _ __ _   _ ')
    print(' / __|/ _ \|  __|  __| | | |')
    print(' \__ \ (_) | |  | |  | |_| |')
    print(' |___/\___/|_|  |_|   \__, |')
    print('                       __/ |')
    print('                      |___/ ')
    print('Generated files have been deleted.')
    executionTime = (time.time() - startTime)
    print('Script run time: ' + str(executionTime) + ' seconds')
    print('========================================= ListExtract.py ==========================================\n')
    sys.exit()

else:
    pass
print('========================================= ListExtract.py ==========================================\n')

#os.remove(f'{searchlist}_fullnames')   