import subprocess
import argparse
import re
import sys
import shutil
import os

def get_options():
    parser = argparse.ArgumentParser(description="Run pipeline to analyse flanking regions of a focal gene.")
    parser.add_argument("--contigs", required=True,
                        help="fasta with contigs containing central/focal gene")
    parser.add_argument("--gene_fasta", required=True,
                        help="fasta with nucleotide sequence of focal gene")
    parser.add_argument("--focal_gene_name", required=True,
                        help="name of focal gene\n(NOTE: if using gffs, must match name of protein product e.g. IMP-4 not blaIMP-4)")
    parser.add_argument("--flanking_region", required=False,
                        help="size of flanking region (N.B. currently symmetrical upstream/downstream). Default: 5000",
                        default="5000")
    parser.add_argument("--output_dir", required=False,
                        help="output directory. Default: output",
                        default="output")
    parser.add_argument("--force", required=False,
                        help="whether to overwrite existing input files. Default: False",
                        action='store_true')
    parser.add_argument("--gff", required=False,
                        help="file with gff annotations for contigs. Default: none",
                        default="")
    parser.add_argument("--panx_export", required=False,
                        help="whether to export panX output from pangraph. Default: False",
                        action="store_true")
    parser.add_argument("--bandage", required=False,
                        help="whether to run Bandage on pangraph. Default: False",
                        action="store_true")
    parser.add_argument("--snv_threshold", required=False,
                        help="SNV threshold for focal gene. Default: 25",
                        default="25")
    parser.add_argument("--pangraph_polish", required=False,
                        help="whether to polish the pangraph. Default: False",
                        action="store_true")
    parser.add_argument("--pangraph_aligner", required=False,
                        help="aligner to use for building pangraph. Default: minimap2",
                        default="minimap2",
                        choices=["minimap2", "mmseqs"])
    parser.add_argument("--pangraph_seed", required=False,
                        help="random seed for pangraph (for reproducibility). Default: 0",
                        default="0")
    parser.add_argument("--pangraph_alpha", required=False,
                        help="value of alpha parameter for pangraph. Default: 100",
                        default="100")
    parser.add_argument("--pangraph_beta", required=False,
                        help="value of beta parameter for pangraph. Default: 10",
                        default="10")
    parser.add_argument("--pangraph_dist_backend", required=False,
                        help="distance backend for calculation of pangraph. Default: native",
                        default="native",
                        choices=["native", "mash"])
    parser.add_argument("--pangraph_minblocklength", required=False,
                        help="minimum block length for pangraph. Default: 100",
                        default="100")
    parser.add_argument("--pangraph_edgeminlength", required=False,
                        help="minimum edge length for pangraph when exporting gfa. Default: 0",
                        default="0")
    parser.add_argument("--breakpoint_minimap2", required=False,
                        help="whether to also calculate breakpoint distances using minimap2. Default: False",
                        action="store_true")
    return parser.parse_args()

def run_snakemake_pipeline(config_file):
    snakemake_command = f"snakemake --cores 1 --configfile {config_file} -r prepare_DB run_pangraph calculate_distances make_plots"
    print(snakemake_command)
    subprocess.run(snakemake_command, shell=True)

default_options =  {"version":"default",
                    "focal_genes": "[]",
                    "complete": "True",
                    "version" : "default",
                    "panx_export": "False",
                    "bandage": "True",
                    "snv_threshold": "25",
                    "region_upstream": "5000",
                    "region_downstream": "5000",
                    "pangraph_polish": "False",
                    "pangraph_aligner": "minimap2",
                    "pangraph_seed":"0",
                    "pangraph_alpha": "100",
                    "pangraph_beta": "10",
                    "pangraph_dist_backend": "native",
                    "pangraph_minblocklength": "100",
                    "pangraph_edgeminlength": "0",
                    "include_gff": "False",
                    "output_dir": "output",
                    "breakpoint_minimap2":"False"}
boolean_options = ["complete", "panx_export", "pangraph_polish",
                    "include_gff", "breakpoint_minimap2", "bandage"]


def write_config_file(user_options, config_file):
    with open(config_file, 'w') as f:
        for key in default_options.keys():
            if key in user_options.keys():
                if key=="focal_genes":
                    f.write("%s: [\"%s\"]\n" % (key, user_options[key]))
                elif key in boolean_options:
                    f.write("%s: %s\n" % (key, user_options[key]))
                else:
                    f.write("%s: \"%s\"\n" % (key, user_options[key]))
            else:
                if key in boolean_options:
                    f.write("%s: %s\n" % (key, default_options[key]))
                else:
                    f.write("%s: \"%s\"\n" % (key, default_options[key]))


def rewrite_gene_fasta(gene_fasta, new_prefix, new_fasta):
    with open(gene_fasta, 'r') as f:
        gene_fasta_string = f.read()
    if gene_fasta_string.count(">")>1:
        print("Warning: gene fasta "+gene_fasta+" seems to have more than one fasta header")
        sys.exit(1)
    elif gene_fasta_string.count(">")==1:
        re.sub(">", ">"+new_prefix+" ", gene_fasta_string)
        with open(new_fasta, "w") as f:
            f.write(gene_fasta_string)

if __name__ == "__main__":

    # Parse command line arguments
    args = get_options()

    user_options = {"output_dir": args.output_dir, 
                    "focal_genes": args.focal_gene_name,
                    "panx_export": args.panx_export,
                    "bandage": args.bandage,
                    "snv_threshold": args.snv_threshold,
                    "region_upstream": args.flanking_region, 
                    "region_downstream": args.flanking_region,
                    "pangraph_polish": args.pangraph_polish,
                    "pangraph_aligner": args.pangraph_aligner,
                    "pangraph_seed":args.pangraph_seed,
                    "pangraph_alpha":args.pangraph_alpha,
                    "pangraph_beta": args.pangraph_beta,
                    "pangraph_dist_backend": args.pangraph_dist_backend,
                    "pangraph_minblocklength": args.pangraph_minblocklength,
                    "pangraph_edgeminlength": args.pangraph_edgeminlength,
                    "breakpoint_minimap2": args.breakpoint_minimap2}

    if not os.path.exists(args.output_dir):
        print("doesn't exist")
        os.makedirs(args.output_dir)

    if os.path.exists(args.output_dir+"/"+args.focal_gene_name):
        print("Warning: output directory "+args.output_dir+"/"+args.focal_gene_name+" exists.")
        print("If you encounter snakemake errors, rerun the snakemake command below and add --unlock.")
        if args.force:
            print("You used --force! Removing previous results and rerunning pipeline.")
            shutil.rmtree(args.output_dir+"/"+args.focal_gene_name)


    config_file = "configs/"+"config_"+args.focal_gene_name+".yaml"

    # copy files into input directory if needed 
    # (wasteful, but because of how I've written snakemake pipeline to function with multiple genes)
    contigs_file = "input/contigs/"+args.focal_gene_name+"_contigs.fa"
    gene_file = "input/focal_genes/"+args.focal_gene_name+".fa"

    if args.gff!="":
        user_options["include_gff"] = "True"
        gff_file = "input/gffs/"+args.focal_gene_name+"_annotations.gff"
        if args.gff==gff_file:
            pass
        else:
            shutil.copy(args.gff, gff_file)

    if args.contigs==contigs_file:
        pass
    else:
        shutil.copy(args.contigs, contigs_file)

    if args.gene_fasta==gene_file:
        with open(args.gene_fasta, "r") as f:
            header = f.readline().split(" ")[0]
            if header!=args.focal_gene_name:
                print("Warning: your focal gene fasta does not have the same name as the focal gene. Please amend.")
                sys.exit(1)
    else:
        rewrite_gene_fasta(args.gene_fasta, args.focal_gene_name, gene_file)
    
    # Write the config file
    write_config_file(user_options, config_file)

    # Call the Snakemake pipeline
    run_snakemake_pipeline(config_file)
