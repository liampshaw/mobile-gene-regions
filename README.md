# Visualizing and quantifying structural diversity around mobile AMR genes

__Liam P. Shaw__

This repository contains a pipeline to analyse the flanking regions of genes using [pangraph](https://github.com/neherlab/pangraph). 

It was developed for the analysis of mobile AMR genes, but in principle can be used for any input gene. 

For more details, see the paper:

Visualizing and quantifying structural diversity around mobile AMR genes  
Liam P. Shaw and Richard A. Neher, *biorxiv* (2023)  
doi: [10.1101/2023.08.07.551646](https://doi.org/10.1101/2023.08.07.551646)

### Dependencies

The pipeline requires [pangraph](https://github.com/neherlab/pangraph) to be installed first. Currently it uses the `feat/unbalanced-tree` branch which can be installed with:

```
git clone -b feat/unbalanced-tree https://github.com/neherlab/pangraph.git && cd pangraph
```

And then follow the installation instructions from the pangraph github as normal.

Other dependencies can be installed with conda or mamba (recommended: mamba):

```
mamba env create -f flanking-regions.yaml 
conda activate flanking-regions
```

There are also R packages required for some plots. These can be installed with:

```
Rscript scripts/install_libraries.R
```


### Usage

Basic usage (run from the current directory):

```
python analyze-flanking-regions.py \
        --contigs contig_fasta.fa \
        --gene_fasta focal_gene.fa \
        --focal_gene_name gene_name \
        --output_dir output
```

The underlying pipeline has flexible parameters such as the size of the flanking region, the number of single nucleotide variants (SNVs) allowed in the focal gene, and pangraph parameters:


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

`analyze-flanking-regions.py` is just a helper script that creates a config file and calls the pipeline for a single gene. However, the pipeline is written in [snakemake](https://snakemake.readthedocs.io/en/stable/index.html), so as to potentially enable analysis of multiple genes in a single command once the input files have been formatted correctly and put in the `input` directory.

Outputs for a given gene are saved in `{output_dir}/{focal_gene_name}`. They include output data files as well as interactive html plots (generated with [altair](https://altair-viz.github.io/) in Python). 

A tutorial on how to prepare data and use the pipeline to explore flanking regions is available in `tutorial/Tutorial.md`.

### Beta-lactamases example

The analysis presented in the paper for twelve beta-lactamases uses the version of the pipeline in `beta-lactamases`. The plots generated can be seen at [https://liampshaw.github.io/flanking-regions/](https://liampshaw.github.io/flanking-regions/).





