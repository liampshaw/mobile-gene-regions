import treetime
from Bio import Phylo
from Bio import SeqIO
import argparse
import numpy as np

# Heavily indebted to code by Marco Molari - see notes on source throughout
# https://github.com/mmolari/pangenome-evo/blob/dev/local/

def get_options():
    parser = argparse.ArgumentParser(description='Convert tree for plotting')
    parser.add_argument('--tree_in', help='input tree file', type=str) 
    parser.add_argument('--aln', help='alignment that was used to generate tree', type=str) 
    parser.add_argument('--tree_out', help='output tree file', type=str) 
    return parser.parse_args()



def rescale_branch_length(node, factor):
	# Source: https://github.com/mmolari/pangenome-evo/blob/dev/local/scripts/pangraph/refine_coretree_treetime.py
    if node.branch_length is not None:
        node.branch_length *= factor
    for child in node:
        rescale_branch_length(child, factor)


def count_snps(aln_dict):
	# Given a SeqIO dict of MSA, returns number of SNPs
	# Adapted from reduce_aln function 
	# Source: https://github.com/mmolari/pangenome-evo/blob/dev/local/scripts/pangraph/reduced_core_alignment.py
	L = len(aln_dict[list(aln_dict.keys())[0]])
	N = len(aln_dict.keys())
	M = np.empty((N, L), dtype=str)
	for n, s in enumerate(list(aln_dict.keys())):
		M[n] = list(aln_dict[s])
	# remove gaps
	has_gap = np.any(M == "-", axis=0)
	n_gaps = has_gap.sum()
	M = M[:, ~has_gap]
	# remove consensus
	is_consensus = np.all(M == M[0], axis=0)
	n_consensus = is_consensus.sum()
	M = M[:, ~is_consensus]
	# remove extra characters
	is_extra = ~np.all(np.isin(M, list("acgtACGT")), axis=0)
	n_extra = is_extra.sum()
	M = M[:, ~is_extra]
	return int(M.shape[1])


def main():
	# Source: https://github.com/mmolari/pangenome-evo/blob/dev/local/scripts/pangraph/refine_coretree_treetime.py
	args = get_options()

	tree_in = args.tree_in
	aln = args.aln
	tree_out = args.tree_out	

	seqs = SeqIO.to_dict(SeqIO.parse(aln, 'fasta'))

	aln_L = max([len(v.seq) for k, v in seqs.items()]) # get maximum length of aln
	aln_snps = count_snps(seqs)

	tree = Phylo.read(tree_in, format="newick")
	rescale_branch_length(tree.root, aln_snps)
	tree.root_at_midpoint() # By default root to midpoint; can reroot on preferred when plotting
	tree.ladderize()


	# instantiate treetime
	myTree = treetime.TreeAnc(
	    gtr="Jukes-Cantor",
	    tree=tree,
	    aln=aln,
	    verbose=0,
	    seq_len=aln_L
	)

	myTree.tree.root.branch_length = 0.0
	# optimize branch length
	print("--------- Treetime ---------")
	print("    < before optimizing >")
	print("tot. branch length:", myTree.tree.total_branch_length())
	print("n. nonterminals:", len(myTree.tree.get_nonterminals()))
	myTree.optimize_tree(prune_short=True)
	print("    < after optimizing >")
	print("tot. branch length:", myTree.tree.total_branch_length())
	print("n. nonterminals:", len(myTree.tree.get_nonterminals()))


	# Scale so that total branch length corresponds to total SNPs
	# i.e. branch length of 1 -> 1 SNP
	tree = myTree.tree 
	rescale_branch_length(tree.root, aln_snps/tree.total_branch_length())


	Phylo.write(
	    myTree.tree,
	    tree_out,
	    format="newick",
	    # format_branch_length="%1.10f",
	    format_branch_length="%.5e",
	)

if __name__=="__main__":
	main()


