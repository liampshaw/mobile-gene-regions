from Bio import SeqIO
import sys
import re
seqs = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))

# We handle the gene name as just everything before "-" and make uppercase for CARD
# This may not work depending on the gene
gene = re.sub("-.*", "", str(sys.argv[2]).upper())
db = sys.argv[3]

genes_seen = []
with open(sys.argv[4], 'w') as f:
	for seqid in seqs:
		if db=="CARD":
			gene_id = seqid.split("|")[-1]
		elif db=="NCBI":
			gene_id = re.sub("bla", "", seqid.split("|")[5])
		if gene+"-" in gene_id:  # addition is for the numbering, to avoid e.g. CMY2 being a match for "CMY" search string 
			if gene_id not in genes_seen: # Need to check for duplicates - if we have already found a sequence, we don't add another (because of duplicate names)
				genes_seen.append(gene_id)
				f.write(">"+gene_id+"\n")
				f.write(str(seqs[seqid].seq)+"\n")
