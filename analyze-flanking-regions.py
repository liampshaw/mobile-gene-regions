import subprocess
import os
import argparse


def get_options():
    parser = argparse.ArgumentParser(description="Run Snakemake pipeline")
    parser.add_argument("--contigs", required=True,
                        help="fasta with contigs containing central/focal gene")
    parser.add_argument("--gene", required=True,
                        help="central/focal gene fasta")
    parser.add_argument("--prefix", required=True,
                        help="prefix to use for gene in config file")
    parser.add_argument("--flanking_region", required=False,
                        help="size of flanking region")
    parser.add_argument("--force", required=False,
                        help="whether to overwrite existing input files",
                        action='store_true')
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
                    "DB": "CARD",
                    "include_gff": "False",
                    "output_prefix": "output",
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

        # f.write("version: \"default\"\n")
        # f.write("focal_genes : [\""+prefix+"\"]\n")
        # f.write("complete: True\n")
        # f.write("panx_export : False\n")
        # f.write("bandage : True\n")
        # f.write("snv_threshold : \"25\"\n")
        # f.write("region_upstream : \"5000\"\n")
        # f.write("region_downstream : \"5000\"\n")
        # f.write("pangraph_polish : False\n")
        # f.write("pangraph_aligner : \"minimap2\"\n")
        # f.write("pangraph_seed : 0\n")
        # f.write("pangraph_alpha : 100\n")
        # f.write("pangraph_beta : 10\n")
        # f.write("pangraph_dist_backend: \"native\"\n")
        # f.write("pangraph_minblocklength : \"100\"\n") # minimum length of blocks 
        # f.write("pangraph_edgeminlength : \"0\"\n")
        # f.write("DB: [\"CARD\"]\n") 
        # f.write("include_gff: False\n")
        # f.write("output_prefix : \"output_"+prefix+"\"\n")
        # f.write("breakpoint_minimap2: True") 

if __name__ == "__main__":

    # Parse command line arguments
    args = get_options()

    config_file = "configs/"+"config_"+args.prefix+".yaml"
    contigs_file = "input/contigs/"+args.prefix+"_contigs.fa"
    gene_file = "input/focal_genes/"+args.prefix+".fa"
    if os.path.isfile(contigs_file):
        if not args.force:
            print("Warning: contigs file "+contigs_file+" exists already. Use --force to overwrite")
        else:
            os.remove(contigs_file)
            os.symlink(args.contigs, contigs_file)
    else:
        os.symlink(args.contigs, contigs_file)
    if os.path.isfile(gene_file):
        if not args.force:
            print("Warning: gene file "+gene_file+" exists already. Use --force to overwrite")
        else:
            os.remove(gene_file)
            os.symlink(args.gene, gene_file)
    else:
        os.symlink(args.gene, gene_file)
    #else:
    #    print("(File with prefix exists in input/contigs, not overwriting)")
    #if not os.path.isfile(gene_file):
    #os.symlink(args.gene, gene_file)
    #else:
    #print("(File with prefix exists in input/focal_gene, not overwriting)")
    #shutil.copyfile(args.contigs, "input/contigs/"+args.prefix+"_contigs.fa")
    #shutil.copyfile(args.gene, "input/focal_genes/"+args.prefix+".fa")
    
    write_config_file({"prefix": args.prefix, "focal_genes": args.prefix,
        "region_upstream": args.flanking_region, "region_downstream": args.flanking_region}, config_file)

    # Call the Snakemake pipeline
    run_snakemake_pipeline(config_file)
