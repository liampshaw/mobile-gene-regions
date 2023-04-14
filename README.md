# Understanding the evolution of the regions surrounding beta-lactamase genes 

__Liam P. Shaw__

Repository of scripts and data for analysis of the genomic contexts of different clinically important beta-lactamases genes.

Using CARD prevalence data, we obtain sequences containing beta-lactamases from clinically important beta-lactamase families and analyse the surrounding regions of genes using [pangraph](https://github.com/neherlab/pangraph).

### Repository structure

Scripts for the analysis are provided in `scripts`. Some of the necessary data files are provided in `data`, but for reasons of size complete fastas are not. To rerun the analysis these need to be downloaded with (e.g.) `ncbi-acc-download` (see `data` folder).

Assuming fasta files have been downloaded and put into `fastadir` in the config file, the complete analysis can be run using `snakemake` (see [here](https://snakemake.readthedocs.io/en/stable/index.html)). This requires a number of dependencies listed in `betalactamase-regions.yaml`. 

There are two possible config files in `configs`. Recommended for reanalysis is `configs/laptop_config.yaml` (the other is for analysis on a slurm cluster). (to do: make pipeline more coherent and modular)

The directory `manuscript` contains manuscript files (to be added).

Outputs are written to `output`. 
