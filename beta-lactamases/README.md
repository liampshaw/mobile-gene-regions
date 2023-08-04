This folder contains a version of the pipeline specifically for analysing the twelve-betalactamase focal genes, as shown in the paper. 

### Getting the input data

The input data files are large and so must first be downloaded from Zenodo: [https://doi.org/10.5281/zenodo.8208376](https://doi.org/10.5281/zenodo.8208376).

The compressed archive should then be extracted within `input` in this folder:

```
cd input
tar -xzvf betalactamase_dataset.tar.gz
```

### Running the pipeline

```
snakemake --cores 1 --configfile betalactamases_config.yaml -r prepare_DB run_pangraph calculate_distances make_plots
```

This takes about 45 minutes to run (on my laptop: 2 GHz Quad-Core Intel Core i5 processors, 16 GB RAM). It should be possible to specify more cores, but currently this isn't working due to some error with calling blast (I hope to fix!).

### Further analysis

The folder `analysis` contains downstream analysis to make the figures in the paper from the pipeline outputs.

