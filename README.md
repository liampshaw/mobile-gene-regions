# Understanding the evolution of the regions surrounding beta-lactamase genes 

__Liam P. Shaw__

Repository of scripts and data for analysis of the genomic contexts of a given focal gene - presents a pipeline to analyse the surrounding regions of genes using [pangraph](https://github.com/neherlab/pangraph) and generate output plots. The original intention was to analyse beta-lactamase genes, but the pipeline is now generalised for any input gene. 

![](tutorial/images/ecdf_screenshot.png)

![](tutorial/images/linear_blocks_screenshot.png)



### Repository structure

The pipeline uses [snakemake](https://snakemake.readthedocs.io/en/stable/index.html). 

Dependencies can be installed with conda or mamba (recommended: mamba):

```
mamba env create -f bl-region-env.yaml
conda activate bl-region-env
```

For a given 'focal' gene, there are two input files:

* `input/focal_genes/{gene}.fa` - the focal gene of interest 
* `input/focal_genes/{gene}_contigs.fa` - multi-fasta with the contigs which have been identified as containing the focal gene 

A small dataset is provided for the *mcr-1.1* gene with 5kb flanking regions. Parameters for the analysis can be specified in the config file: 

Within the `Snakefile` genes are specified in a list of `FOCAL_GENES`:

```
FOCAL_GENE_DICT = {"mcr-1.1"}
```

The full snakemake pipeline can be run with the following command:

```
snakemake --cores 1 \
          --configfile configs/default_config.yaml \
          -r run_pangraph calculate_distances make_plots
```

Outputs are written to `output` and include output data files as well as pdf plots (generated with R) and html plots (generated with altair in Python):

![](tutorial/images/example_screenshot.png)

A more detailed tutorial on how to prepare data and use the pipeline is available in `tutorial/Tutorial.md`.
