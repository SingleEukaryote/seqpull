import sys
from Bio import SeqIO

# Run: "python SequenceCleaner.py {input file} {minimum sequence length} {percent Ns allowed} {output file}"

fasta_file = sys.argv[1]
min_length = sys.argv[2]  # The minimum sequence length to keep (if 500, it'll remove all sequences below 500bp)
por_n = sys.argv[3]       # Default should be 100, to keep all sequences
outputfile = sys.argv[4]

#def sequence_cleaner(fasta_file, min_length=0, por_n=100):
    # Create our hash table to add the sequences
    
sequences = {}

# Using the Biopython fasta parse we can read our fasta input
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    # Take the current sequence
    sequence = str(seq_record.seq).upper()
    # Check if the current sequence is according to the user parameters
    if (
        len(sequence) >= int(min_length)
        and (float(sequence.count("N")) / float(len(sequence))) * 100 <= int(por_n)
    ):
        # If the sequence passed in the test "is it clean?" and it isn't in the
        # hash table, the sequence and its id are going to be in the hash
        if sequence not in sequences:
            sequences[sequence] = seq_record.id
        # If it is already in the hash table, we're just gonna concatenate the ID
        # of the current sequence to another one that is already in the hash table
        else:
            sequences[sequence] += "_" + seq_record.id

# Write the clean sequences

# Create a file in the same directory where you ran this script
with open(f'{outputfile}', "w+") as cleanedFile:
    for sequence in sequences:
        cleanedFile.write(">" + sequences[sequence] + "\n" + sequence + "\n")




#sequence_cleaner(sys.argv[1],sys.argv[2],sys.argv[3])
name = fasta_file.split('/')[-1]
output = open(f'{outputfile}', "r")
data = output.read()
seqsafter = data.count('>')
output.close()
name = fasta_file.split('/')[-1]
output = open(f'{fasta_file}', "r")
data = output.read()
realnumseqs = data.count('>')
output.close()
# log
print('\n=======================================================================================')
print(f'> Removing sequences < {min_length} from {name}')
print(f'> Amount of sequences before: {realnumseqs}')
print(f'> Amount of sequences after:  {seqsafter}')
print('=======================================================================================\n')
