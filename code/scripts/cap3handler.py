#!/usr/bin/env python 
# for CAP3
# Descript: Converts multiline FASTAs to single line FASTAs
# Also checks if files are 
# Usage: FastaMLtoSL.py <sequences.faa> 
# Example: FastaMLtoSL.py mySeqs.faa 
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports:
import os
import sys
import re
import time
#===========================================================================================================

# Stores file one for input checking.
inFile  = sys.argv[1]
outFile = sys.argv[2]
stall = sys.argv[3] # Input a number of seconds to stall the script for. Use 0 seconds if manual. 
# Note there is a bug with CAP3 where I need this script to stall for a couple seconds until the contig file is full


filename = inFile.split('/')[-1]
print('\n===============================================================================================')
print(f'> Collapsing contigs of {filename}')


## Stall time ###
print(f'>> Stalling for {stall} second(s), please remove if you do not think this step is necessary.')
time.sleep(int(stall))
#################


# Reads sequence file list and stores it as a string object. Safely closes file:
try:	
	with open(inFile,"r") as newFile:
		sequences = newFile.read()
		sequences = re.split("^>", sequences, flags=re.MULTILINE) # Only splits string at the start of a line.
		del sequences[0] # The first fasta in the file is split into an empty empty element and and the first fasta
						 # Del removes this empty element.
		newFile.close()
except IOError:
	print("Failed to open " + inFile)
	exit(1)
# Conversts multiline fasta to single line. Writes new fasta to file.
try:	
	with open(outFile,"w") as newFasta:
		for fasta in sequences:
			try:
				header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
			except ValueError:
				print(fasta)
			header = ">" + header + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
			sequence = sequence.replace("\n","") + "\n" # Replace newlines in sequence, remember to add one to the end.
			newFasta.write(header + sequence)
		newFasta.close()
except IOError:
	print("Failed to open " + inFile)
	exit(1)

# check if size of file is 0
if os.stat(outFile).st_size == 0:
    outpath = outFile.split('1')[0]
    # ../output/{gene}/{filename}/1.{filename}.{gene}.fas.cap.contigs.oneline.fas"
    print('>>ERROR, CAP3 assembly did not work. This is likely because the input sequences are too long or too few to require assembly.')
    print('\n')
    os.system(f'rm {outpath}/*cap')
    os.system(f'rm {outpath}/*cap*')
    #open(f"{outpath}/_CAP3.assembly.failed.txt", "w")
    exit()
else:
    print(">>> Success!")
print('===============================================================================================\n')