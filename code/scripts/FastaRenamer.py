import sys
import os
'''
1. Inserts sys.argv[3] into the sequence header. It also appends the sequence length to the title.
2. python3 fastarenamer.py {input file} {output file name} {organism name to add to sequence header}

- Appends the sequence name and length to the end of the filename provided.
'''
inFile = sys.argv[1]
outFile = sys.argv[2]
filename = sys.argv[3]

filesize = os.path.getsize(inFile)
if filesize == 0:
    print("\n\n======================== fastarenamer ================================\n")
    print(f"{inFile} is empty: " + str(filesize) + "\n Deleting empty file.")
    print("\n========================================================================\n\n")
    os.remove(inFile)
    sys.exit()
else:
    pass

name = inFile.split('/')[-1]
print('\n===============================================================================================')
print(f'> Renaming sequences of {name}\n')

print('New sequence names:')

####################################### Rename "Contig" or "Trinity" ###############################################
renameFile = open(f'{inFile}', 'r')
renamedFile = open(f'{outFile}', 'w')
for line in renameFile:
    if line.startswith('>'):
        seq = renameFile.readline()
        name = line.replace(" ", "|")
        ################### MOST IMPORTANT LINE
        ################### Modify newheader to name new sequences according to the input file convention
        newheader = '>' + filename + '.' + name.split('>')[1].split("\n")[0] + '.len' + str(len(seq)) + "\n"
        
        #######################################\
        print('\nFrom ' + line.split('\n')[0])
        print('To   ' + newheader.split('\n')[0])
        renamedFile.write(newheader)
        renamedFile.write(seq)  
renameFile.close()
renamedFile.close()
##########################################################################################################

print('===============================================================================================\n')
