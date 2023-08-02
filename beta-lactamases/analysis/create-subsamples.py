# Make subsamples of CTX-M-65 dataset

from Bio import SeqIO
from random import sample

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

def read_fasta(fasta_file):
    '''Simply uses SeqIO to read in fasta as dict.'''
    return(to_dict_remove_dups(SeqIO.parse(fasta_file, 'fasta')))

fasta_file = "/Users/Liam/Downloads/test-bl/input/contigs/CTX-M-65_contigs.fa"

seqs = read_fasta(fasta_file)

subsamplings = [10, 50, 100, 200, 300, 400]
always_include = ["NZ_CP047194.1", "NZ_CP066256.1"]

for s in subsamplings:
    subsampled = sample(sorted(seqs), s)
    subsampled += always_include
    subsampled = set(subsampled)
    with open("CTX-M-65_subsampled_"+str(s)+".fa", "w") as f:
        for sub in subsampled:
            f.write(">%s\n%s\n" % (sub, str(seqs[sub].seq)))
    

