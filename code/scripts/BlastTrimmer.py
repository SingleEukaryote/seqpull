import sys
import os
from os.path import exists
import time
startTime = time.time()
'''
Uses qstart and qend from blast result tsvs (-outfmt '6 qseqid qstart qend')
Trims sequences based on the lowest blast query start site and highest blast query end site

Run:
python blast_trimmer.py blastresults.tsv sequencesfile.fas outputfilename.fas
'''

########## INPUT ###########
blastresult = sys.argv[1]
############   blastresults.tsv (output format: -outfmt '6 qseqid qstart qend')
############   qseqid must the be first column.
sequencesfile = sys.argv[2]
#############   sequencesfile.fas (fasta file you queried against the database)
speciesnames = sys.argv[3] 
# This file is PR2 18S pipeline specific, remove if you are use this script for other purposes
# Make 'null' if not in use
#############
outputfile = sys.argv[4]
##############   outputfile.fas (the name of the file that contains the trimmed sequences.)
############################

## Log
fastaname = sequencesfile.split('/')[-1]
listname = blastresult.split('/')[-1]
print('\n============================================== BlastTrimmer.py ==============================================')
print(f'> Trimming sequences in {fastaname} \n> Using these blast results: {listname}\n')
## 
################################################################################################
###
###### Parse blast results tsv and make a dictionary of start and end sites per sequence #######
###
################################################################################################


# Creating an empty dictionary
SEDict = {}


# Fill dictionary with sequence names
blastresults = open(f'{blastresult}', 'r') # Read blastresults.tsv
for line in blastresults:
    seqname = line.split('\t')[0]
    SEDict[f"{seqname}"] = []
    
blastresults.close()

# Define sequences with start and end numbers

blastresults = open(f'{blastresult}', 'r') 
for line in blastresults:
    
    seqname = line.split('\t')[0]
    qstart = line.split('\t')[1]
    qend = line.split('\t')[2].split('\n')[0]

    SEDict[f"{seqname}"].append(int(qstart))
    SEDict[f"{seqname}"].append(int(qend))


#print('\n\n--- Trim Dictionary ---\n')
#pprint.pprint(SEDict)
#print('\n-----------------------\n\n')

## Make a set of contigs


#####################################################################################################
###
###### Use dictionary to fish out and trim sequences based on min and max numbers in the list #######
###
#####################################################################################################



linenum = 0
onelinefile = f'{sequencesfile}.tmp'
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
 
# Generate a set of the sequences you wish to retrieve from the dictionary keys
desired_seqs = set(SEDict.keys())
 
# Find the overlap between the total sequences and the desired ones
all_seqs_names = set(all_seqs.keys())
toextract = all_seqs_names.intersection(desired_seqs)

# Use 'toextract' set to generate desired file

if speciesnames == 'null':
    pass
else:
    #############################################################################
    ######## PR2 SPECIES APPENDING, REMOVE IF NOT USING PR2 18S PIPELINE ########
    #############################################################################
    # Creating an empty dictionary
    specieDict = {}

    # Fill dictionary with sequence names
    speciename = open(f'{speciesnames}', 'r') # Read speciename.tsv
    for line in speciename:
        seqname = line.split('\t')[0]
        specieDict[f"{seqname}"] = []
    speciename.close()
    
    # Define sequences with specie names

    speciename = open(f'{speciesnames}', 'r') 
    for line in speciename:
        
        seqname = line.split('\t')[0]
        specie = line.split('\t')[1].split('\n')[0]

        specieDict[f"{seqname}"].append(specie)

    #print('\n\n--- Specie Dictionary ---\n')
    #print(specieDict)
    #print('\n--------------------------\n\n')

    #############################################################################
    #############################################################################
    #############################################################################
bpscut = 0
with open(outputfile, 'w') as extractfile:
    for name in toextract:

        ############################### Trim module ##################################
        ### USE Start End DICTIONARY TO GET START AND END SITES FOR THE SEQUENCE NAME
        seqlist = SEDict[name]
        start = min(seqlist) - 1
        end = max(seqlist)
        ### SEQUENCE DICTIONARY, chooses the sequence correlating to the contig name
        seq = all_seqs[name]
        ### TRIM SEQUENCES USING A RANGE
        trim = seq[start:end]
        ### SEQUENCE LENGTHS
        olen = len(seq)
        nlen = len(seq[start:end])

        ##############################################################################
        if speciesnames == 'null':
            print('--------------------------------------------------') 
        else:
            ############################# PR2 specie module ##############################
            # Defines specie for contig using specie name dictionary from above
            # >bf222b3.pool25.18S.Contig1.len5886.revcomp
            if name.endswith('revcomp'):
                norevcomp = name.split('.revcomp')[0]
                specie = specieDict[str(norevcomp)]
            else:
                specie = specieDict[str(name)]

            for line in open(f'{speciesnames}', 'r'):
                contig = line.split('\t')[0]

                ## See if the contig name is equal
                if contig.split('.')[3] == name.split('.')[3]:
                    test1 = contig.split('.')[3]
                    test2 = name.split('.')[3]

                    if test1 == test2:
                        pass
                        #print('--------------------------------------------------')       
                        #print('test1 = ' + test1)
                        #print('test2 = ' + test2) 
                        #print('GOOD! ' + test1 + ' = ' + test2)
                        
                    else:
                        print(f'===========ERROR============= {test2} ============ERRROR=============')
                        print('test1 = ' + test1)
                        print('test2 = ' + test2) 
                        print('WARNING!!!!!!!!!!!!!!!!' + '.\n!\n!\n!\n!\n!\n!\n! TELL TRISTAN TO FIX THE DAMN SCRIPT\n!\n!\n!\n!\n!\n!\n!\n!\n')



            ################### REMOVE IF NOT USING PR2 18S PIPELINE #####################

        ################# Write new file with trimmed sequences ######################
        if speciesnames == 'null':
            extractfile.write('>' + name + '.nlen' + str(nlen) + '\n' + trim + '\n')
            print(name + "\n")
        else:
            extractfile.write('>' + name + '.nlen' + str(nlen) + str(specie).split("'")[1] + '\n' + trim + '\n')
            print(name + "\n" + str(specie).split("'")[1].split('_')[1] + "% " + str(specie).split("'")[1].split('_')[2] + ' ' + str(specie).split("'")[1].split('_')[3])
        ##############################################################################

        ##############################################################################

        
        print('Old length:     ' + str(olen))
        print(f'Min start:      {start} (based on blast results)\nMax end:        {end} (based on blast results)')
        print('New length:     ' + str(nlen))
        print(f'Base pairs cut: {olen-nlen}')
        print('--------------------------------------------------') 
        ##############################################################################
        bpscut += (olen-nlen)
extractfile.close()

if exists(onelinefile):
    os.remove(onelinefile)
else:
    print('File empty')


### For pipeline ###
#print('============================================================')
#print('EMPTY FILES?')
#print('============================================================')
#os.system('find ./ -type f -empty -print -delete')
#print('============================================================')
####################



## Input file
output = open(f'{sequencesfile}', "r")
data = output.read()
originalseqs = data.count('>')
output.close()
## Amount of seqs in final file
output = open(f'{outputfile}', "r")
data = output.read()
realnumseqs = data.count('>')
output.close()

print('Input sequences:      ' + str(originalseqs))
print('Output sequences:     ' + str(realnumseqs))
print('Total bp trimmed off: ' + str(bpscut))
executionTime = (time.time() - startTime)
print('Script run time:      ' + str(executionTime) + ' seconds')
print('============================================== BlastTrimmer.py ==============================================\n')
