import argparse
import subprocess
import os
import re

def get_options():
    parser = argparse.ArgumentParser(description='Make a multi fasta file for a given gene.',
                                     prog='make_multi_fasta')
    parser.add_argument('--gene', help='Name of gene (e.g. CTX-M-15, VIM-1, NDM-1).', required=True)
    parser.add_argument('--threshold', help='SNP threshold for distance from gene variant.', required=False, default=25)
    parser.add_argument('--fastadir', help='Directory with fasta files.', required=True)
    parser.add_argument('--outputdir', help='Output directory', required=False, default='./')
    parser.add_argument('--multiplehits', help='Whether to allow multiple hits', required=False, default=False, action='store_true')
    return parser.parse_args()

def main():
    args = get_options()
    if not os.path.exists(args.outputdir):
        os.mkdir(args.outputdir)

    gene_family = "-".join(args.gene.split("-")[:-1])
    print(gene_family)

    # I don't think the next 5 lines are needed...?
    get_seq = "grep -A 1 '"+args.gene+" ' ../data/nucleotide_fasta_protein_homolog_model.fasta | tail -n 1"
    get_seq_out = subprocess.check_output(get_seq, shell=True)
    gene_seq = get_seq_out.decode().strip('\n')
    with open(args.outputdir+'/'+args.gene+'.fa', 'w') as f:
        f.write('>%s\n%s' % (args.gene, gene_seq))

    get_gene_clusters = 'python get_accessions_for_gene.py '+gene_family+' '+\
                                            args.gene+' '+\
                                            str(args.threshold)+\
                                             " --singlehits | cut -d ',' -f 1 | sed -e 's/$/.fa/g' | "+\
                                             "sed -e 's~^~'"+args.fastadir+"'/~g' > "+\
                                             args.outputdir+"/"+args.gene+"_fastas.txt"
    get_gene_clusters_out = subprocess.check_output(get_gene_clusters, shell=True)

    cat = "xargs cat < "+args.outputdir+"/"+args.gene+"_fastas.txt > "+args.outputdir+"/"+args.gene+"_seqs.fa"
    cat_out = subprocess.check_output(cat, shell=True)


if __name__=='__main__':
    main()
