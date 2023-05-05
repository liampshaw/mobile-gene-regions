# Understanding the evolution of the regions surrounding beta-lactamase genes 

__Liam P. Shaw__

Repository of scripts and data for analysis of the genomic contexts of a given focal gene - with a focus on beta-lactamase genes. 

For the beta-lactamases, using CARD prevalence data, we obtain sequences containing beta-lactamases from clinically important beta-lactamase families and analyse the surrounding regions of genes using [pangraph](https://github.com/neherlab/pangraph).

### Repository structure

NEW:

For a given 'focal' gene, there are two input files required:

* `input/focal_genes/{gene}.fa` - the focal gene of interest 
* `input/focal_genes/{gene}_combined_contigs.fa` - multi-fasta with the contigs putatively containing the focal gene 

Example data is provided for MCR-1.1 as downloaded from [MicroBIGG-E](https://www.ncbi.nlm.nih.gov/pathogens/microbigge/#mcr-1.1).

Within the `Snakefile` genes are specified in a `FOCAL_GENE_DICT`:

```
FOCAL_GENE_DICT = {"MCR-1.1":"MCR"} # this syntax is due to e.g. CTX-M-15 being within CTX-M family
```

The full snakemake pipeline can be run with the following command:

```
snakemake --cores 1 --configfile configs/laptop_config.yaml -r prepare_DB run_pangraph calculate_distances make_plots
```


OLD:
Scripts for the analysis are provided in `scripts`. Some of the necessary data files are provided in `data`, but for reasons of size complete fastas are not. To rerun the analysis these need to be downloaded with (e.g.) `ncbi-acc-download` (see `data` folder).

Assuming fasta files have been downloaded and put into `fastadir` in the config file, the complete analysis can be run using `snakemake` (see [here](https://snakemake.readthedocs.io/en/stable/index.html)). This requires a number of dependencies listed in `betalactamase-regions.yaml`. 

There are two possible config files in `configs`. Recommended for reanalysis is `configs/laptop_config.yaml` (the other is for analysis on a slurm cluster). (to do: make pipeline more coherent and modular)

The directory `manuscript` contains manuscript files (to be added).

Outputs are written to `output`. 
