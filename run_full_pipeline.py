import subprocess

def run_snakemake_pipeline(config_file):
    snakemake_command = f"snakemake --cores 1 --configfile {config_file} -r prepare_DB run_pangraph calculate_distances make_plots"
    subprocess.run(snakemake_command, shell=True)

if __name__ == "__main__":
    import argparse

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run Snakemake pipeline")
    parser.add_argument("--config-file", dest="config_file", required=True,
                        help="Path to the Snakemake config file")
    args = parser.parse_args()

    # Call the Snakemake pipeline
    run_snakemake_pipeline(args.config_file)