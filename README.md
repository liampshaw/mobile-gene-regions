# Visualizing and quantifying structural diversity around mobile AMR genes

__Liam P. Shaw__

This repository contains a pipeline to analyse the flanking regions of genes using [pangraph](https://github.com/neherlab/pangraph). 

It was developed for the analysis of mobile AMR genes, but can be used for any input gene. 

Basic usage (run from the current directory):

```
python analyze-flanking-regions.py \
        --contigs contig_fasta.fa \
        --gene focal_gene.fa 
        --prefix gene_name 
```

For more details, see the linked paper:

Visualizing and quantifying structural diversity around mobile AMR genes  
Liam P. Shaw and Richard A. Neher, *biorxiv* (2023)  
doi: []()

### Dependencies

The pipeline requires [pangraph](https://github.com/neherlab/pangraph) to be installed. 

The pipeline uses [snakemake](https://snakemake.readthedocs.io/en/stable/index.html). 

Dependencies can be installed with conda or mamba (recommended: mamba):

```
mamba env create -f flanking-regions.yaml 
conda activate flanking-regions
```

### Usage

The pipeline has flexible parameters such as the size of the flanking region, the number of single nucleotide variants (SNVs) allowed in the focal gene, and pangraph parameters:


```
options:
  -h, --help            show this help message and exit
  --contigs CONTIGS     fasta with contigs containing central/focal gene
  --gene_fasta GENE_FASTA
                        fasta with nucleotide sequence of focal gene
  --focal_gene_name FOCAL_GENE_NAME
                        name of focal gene (NOTE: if using gffs, must match name of protein
                        product e.g. IMP-4 not blaIMP-4)
# Optional parameters
  --flanking_region FLANKING_REGION
                        size of flanking region (N.B. currently symmetrical
                        upstream/downstream). Default: 5000
  --output_dir OUTPUT_DIR
                        output directory. Default: output
  --force               whether to overwrite existing input files. Default: False
  --gff GFF             file with gff annotations for contigs. Default: none
  --panx_export         whether to export panX output from pangraph. Default: False
  --bandage             whether to run Bandage on pangraph. Default: False
  --snv_threshold SNV_THRESHOLD
                        SNV threshold for focal gene. Default: 25
  --pangraph_polish     whether to polish the pangraph. Default: False
  --pangraph_aligner {minimap2,mmseqs}
                        aligner to use for building pangraph. Default: minimap2
  --pangraph_seed PANGRAPH_SEED
                        random seed for pangraph (for reproducibility). Default: 0
  --pangraph_alpha PANGRAPH_ALPHA
                        value of alpha parameter for pangraph. Default: 100
  --pangraph_beta PANGRAPH_BETA
                        value of beta parameter for pangraph. Default: 10
  --pangraph_dist_backend {native,mash}
                        distance backend for calculation of pangraph. Default: native
  --pangraph_minblocklength PANGRAPH_MINBLOCKLENGTH
                        minimum block length for pangraph. Default: 100
  --pangraph_edgeminlength PANGRAPH_EDGEMINLENGTH
                        minimum edge length for pangraph when exporting gfa. Default: 0
  --breakpoint_minimap2
                        whether to also calculate breakpoint distances using minimap2.
                        Default: False

```

A tutorial on how to prepare data and use the pipeline is available in `tutorial/Tutorial.md`.





### Repository structure



`analyze-flanking-regions.py` is a helper script that creates a config file and calls the snakemake pipeline. 


For an example, the snakemake pipeline can be run on some example data (mcr-1) with the following command:

```
snakemake --cores 1 \
          --configfile configs/default_config.yaml \
          -r prepare_DB run_pangraph calculate_distances make_plots
```

Equivalently you can call the full pipeline on a config file with:

```
python run_full_pipeline.py --config-file configs/default_config.yaml
```

Outputs for a given gene are saved in `{output_prefix}/{focal_gene}` and include output data files as well as pdf plots (generated with R) and html plots (generated with [altair](https://altair-viz.github.io/) in Python).



