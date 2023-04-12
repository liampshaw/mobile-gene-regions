from Bio import SeqIO
import sys
import argparse
import os

def get_options():
    parser = argparse.ArgumentParser(description='Splits a fasta file into separate fasta files, one per contig')
    parser.add_argument('--fasta', help='Input multifasta file', type=str)
    parser.add_argument('--outdir', help='Output directory', type=str)
    return parser.parse_args()

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

def read_fasta(fasta_file):
    '''Simply uses SeqIO to read in fasta as dict.'''
    return(to_dict_remove_dups(SeqIO.parse(fasta_file, 'fasta')))

def main():
    args = get_options()
    seqs = read_fasta(args.fasta)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    for seq in seqs:
        seq_name = seq.split('.')[0]
        seq_seq = str(seqs[seq].seq)
        with open(args.outdir+"/"+seq_name+'.fa', 'w') as f:
            f.write('>%s\n%s\n' % (seq_name, seq_seq))

if __name__=="__main__":
    main()
