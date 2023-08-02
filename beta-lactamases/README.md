This folder contains a version of the pipeline specifically for analysing the twelve-betalactamase focal genes, as shown in the paper. 

### Getting the input data

The input data files are large and so must first be downloaded from Zenodo: doi-link-to-be-added 


### Running the pipeline

```
snakemake --cores 1 --configfile betalactamases_config.yaml -r prepare_DB run_pangraph calculate_distances make_plots
```

(may take some time)

### Further analysis

The folder `analysis` contains downstream analysis to make the figures in the paper from the pipeline outputs.

