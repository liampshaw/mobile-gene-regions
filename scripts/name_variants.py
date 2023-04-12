# 

from Bio import SeqIO
import re
import argparse

def get_options():
    parser = argparse.ArgumentParser(description='Assign variants for gene')
   # parser.add_argument('--gene', help='Gene', type=str) # assumes stuff about directories, that pipeline has been run
    parser.add_argument("--variant_fasta", help="Fasta file of variants to use", type=str, required=True)
    parser.add_argument('--output_file', help='output file of assignments', type=str, required=False, default='gene-variants')
    parser.add_argument('--input_fasta',help='Input fasta', required=True)
    parser.add_argument('--outputdir', help='output dir', type=str, default='./', required=False)
    return parser.parse_args()


def main():
	args = get_options()
	variant_seqs = SeqIO.to_dict(SeqIO.parse(args.variant_fasta, 'fasta'))

	#gene = args.gene
	#gene_family = re.sub("-.*", "", gene) # technically doesn't handled CTX-M properly, but ok

	fasta_file = args.input_fasta
	seqs_to_assign = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

	#gene_seqs = {k:v for k, v in variant_seqs.items() if gene_family in k}


	seqs_aa = {str(v.seq.translate()): re.sub('.*\\|', '', k)  for k, v in variant_seqs.items()}

	seqs_nt = {str(v.seq):re.sub('.*\\|', '', k) for k, v in variant_seqs.items()}

	seq_names_dict = {}
	seq_counts_dict = {}
	for seqid, seq_obj in seqs_to_assign.items():
		seq_nt = str(seq_obj.seq)
		seq_aa = str(seq_obj.seq.translate())

		# Check for truncation
		if seq_aa.count('*')>1:
			seq_names_dict[seqid] = 'truncated'
			#print(seqid, 'truncated')
		else:
			if seq_aa in seqs_aa.keys():
				#print(seq_aa)
				seq_names_dict[seqid] = seqs_aa[seq_aa]
			else:
				seq_names_dict[seqid] = 'unnamed'
				#print(seqid, 'unnamed')
		if seq_nt in seq_counts_dict.keys():
			seq_counts_dict[seq_nt] +=1
		else:
			seq_counts_dict[seq_nt] = 1

	#print(seq_names_dict)
	#print(seq_counts_dict)

	with open(args.output_file, "w") as f:
		for seqid, name in seq_names_dict.items():
			f.write("%s,%s\n" % (seqid, name))

	#alignment_file = args.inputprefix+'_focal_gene.dedup.aln'
	#metadata_file = args.inputprefix+'_focal_gene.dedup.txt' 
	#seq_groups = [line.split('\t')[1]  for line in open(metadata_file, 'r').readlines()]

	

	#aligned_seqs_dedup = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))

	# counts_seqs = {}
	# for seqid in aligned_seqs_dedup.keys():
	# 	n_identical_seqs = [len(x.split(',')) for x in seq_groups if seqid in x]
	# 	if len(n_identical_seqs)==0:
	# 		n_identical_seqs = 1
	# 	else:
	# 		n_identical_seqs = n_identical_seqs[0]
	# 	counts_seqs[seqid] = n_identical_seqs

	# Write a nucleotide alignment, with protein names where applicable (different nt seq can have same protein name because)
	# of synonymous mutations
	# with open(args.outputdir+args.outputprefix+'.aln','w') as f:
	# 	for seqid, seq_obj in aligned_seqs_dedup.items():
	# 		seq = re.sub('\n', '', str(seq_obj.seq))
	# 		f.write('>%s | %s %d\n%s\n' % (seqid, seq_names_dict[seqid], counts_seqs[seqid], seq))

	# # Write the sequence variants
	# with open(args.outputdir+args.outputprefix+'.csv','w') as f:
	# 	for seq, name in seq_names_dict.items():
	# 		f.write('%s,%s\n' % (seq,name) )


if __name__=="__main__":
	main()
