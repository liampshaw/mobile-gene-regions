import subprocess
import shutil



def run_snakemake_pipeline(config_file):
    snakemake_command = f"snakemake --cores 1 --configfile {config_file} -r prepare_DB run_pangraph calculate_distances make_plots"
    subprocess.run(snakemake_command, shell=True)

def write_config_file(prefix, config_file):
    with open(config_file, 'w') as f:
        f.write("version: \"default\"\n")
        f.write("focal_genes : [\""+prefix+"\"]\n")
        f.write("complete: True\n")
        f.write("panx_export : False\n")
        f.write("bandage : True\n")
        f.write("snv_threshold : \"25\"\n")
        f.write("region_upstream : \"5000\"\n")
        f.write("region_downstream : \"5000\"\n")
        f.write("pangraph_polish : False\n")
        f.write("pangraph_aligner : \"minimap2\"\n")
        f.write("pangraph_seed : 0\n")
        f.write("pangraph_alpha : 100\n")
        f.write("pangraph_beta : 10\n")
        f.write("pangraph_dist_backend: \"native\"\n")
        f.write("pangraph_minblocklength : \"100\"\n") # minimum length of blocks 
        f.write("pangraph_edgeminlength : \"0\"\n")
        f.write("DB: [\"CARD\"]\n") 
        f.write("include_gff: False\n")
        f.write("output_prefix : \"output_"+prefix+"\"\n")
        f.write("breakpoint_minimap2: True") 

if __name__ == "__main__":
    import argparse

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run Snakemake pipeline")
    parser.add_argument("--contigs", required=True,
                        help="fasta with contigs containing central/focal gene")
    parser.add_argument("--gene", required=True,
                        help="central/focal gene fasta")
    parser.add_argument("--prefix", required=True,
                        help="prefix to use for gene in config file")
    args = parser.parse_args()

    config_file = "configs/"+"config_"+args.prefix+".yaml"

    shutil.copyfile(args.contigs, "input/contigs/"+args.prefix+"_contigs.fa")
    shutil.copyfile(args.gene, "input/focal_genes/"+args.prefix+".fa")
    
    write_config_file(args.prefix, config_file)

    # Call the Snakemake pipeline
    run_snakemake_pipeline(config_file)
