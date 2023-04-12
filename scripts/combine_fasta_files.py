# For a list of accessions
# combine the fasta files into a single fasta 
# so that they can be blasted for the gene

fasta_files = [line.strip("\n") for line in open(snakemake.input, "r").readlines()]
